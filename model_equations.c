// Implementation of the age-structured VSEIR2 model with preschool booster
// Parameter model_vac controls whether aP and wP, or wP and natural infection are assumed similar
// Vaccine coverage is assumed to have a linear between time t0 (start of vaccination) and t1

#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <pomp.h>

// Expit transform for parameters constrained in interval [a,b]
static double expitCons(double x, double a, double b) {
	double out = (a + b * exp(x)) / (1.0 + exp(x));
	if(ISNAN(out)) out = (b + a * exp(-x)) / (1.0 + exp(-x)); // If x=+Inf, must return b
	return out;
}

// Logit transform for parameters constrained in interval [a,b]
static double logitCons(double x, double a, double b) {
	x = (x <= a) ? a : (x >= b ? b : x);
	double out = log((x - a) / (b - x));
	return out;
}

// State variables 
#define V(J) (x[stateindex[0] + (J)]) // Vaccinated
#define S1(J) (x[stateindex[1] + (J)]) // Susceptibles with no previous infection
#define E1(J) (x[stateindex[2] + (J)]) // Susceptibles with at least one previous infection
#define I1(J) (x[stateindex[3] + (J)]) // Exposed with no previous infection
#define S2(J) (x[stateindex[4] + (J)]) // Exposed with at least one previous infection
#define E2(J) (x[stateindex[5] + (J)]) // Individuals with primary infection
#define I2(J) (x[stateindex[6] + (J)]) // Individuals with repeat infection
#define R(J) (x[stateindex[7] + (J)]) // Recovered
#define C1(J) (x[stateindex[8] + (J)]) // Total case count for primary infections
#define C2(J) (x[stateindex[9] + (J)]) // Total case count for repeat infections

// Variations
#define DV(J) (f[stateindex[0] + (J)]) // Vaccinated
#define DS1(J) (f[stateindex[1] + (J)]) // Susceptibles with no previous infection
#define DE1(J) (f[stateindex[2] + (J)]) // Susceptibles with at least one previous infection
#define DI1(J) (f[stateindex[3] + (J)]) // Exposed with no previous infection
#define DS2(J) (f[stateindex[4] + (J)]) // Exposed with at least one previous infection
#define DE2(J) (f[stateindex[5] + (J)]) // Individuals with primary infection
#define DI2(J) (f[stateindex[6] + (J)]) // Individuals with repeat infection
#define DR(J) (f[stateindex[7] + (J)]) // Recovered
#define DC1(J) (f[stateindex[8] + (J)]) // Total case count for primary infections
#define DC2(J) (f[stateindex[9] + (J)]) // Total case count for repeat infections

// Parameters on the natural scale (all rates are per year)
#define iota (p[parindex[0]]) // Background force of infection
#define q1(J) (p[parindex[1] + (J)]) // Probability of infection given contact
#define q2(J) (p[parindex[2] + (J)]) // Age-specific infectiousness
#define omegaC(J) (p[parindex[3] + (J)]) // Seasonality coefficients in 5-10 y
#define omegaT(J) (p[parindex[4] + (J)]) // Seasonality coefficients in 10-20 y
#define CR(I, J) (p[parindex[5] + (I + 15 * J)]) // Contact rates between age groups
#define theta (p[parindex[6]]) // Relative infectiousness of repeat vs. primary infections
#define beta_sd (p[parindex[7]]) // White noise intensity in infection rate
#define sigma (p[parindex[8]]) // Inverse of latency period
#define gamma (p[parindex[9]]) // Inverse of infectious period
#define epsilonV (p[parindex[10]]) // Leakiness of vaccine-induced immunity
#define epsilonI (p[parindex[11]]) // Leakiness of natural infection-induced immunity
#define alphaV (p[parindex[12]]) // Waning rate of vaccine-induced immunity
#define alphaI (p[parindex[13]]) // Waning rate of natural infection-induced immunity
#define rho1(J) (p[parindex[14] + (J)]) // Reporting probability for primary infections, age-specific
#define rho2 (p[parindex[15]]) // Reporting probability for repeat infections
#define tau (p[parindex[16]]) // Reporting overdispersion
#define epsA (p[parindex[17]]) // Primary vaccine failure
#define v1 (p[parindex[18]]) // Vaccine coverage for primary course
#define v2 (p[parindex[19]]) // Vaccine coverage for preschool booster
#define t0 (p[parindex[20]]) // Time at which vaccination is started
#define t1 (p[parindex[21]]) // Time at which the maximal vaccine coverage is reached
#define t2 (p[parindex[22]]) // Date of first data point for vaccine coverage
#define delta(J) (p[parindex[23] + (J)]) // Vector of aging rates in each age group
#define N(J) (p[parindex[24] + (J)]) // Vector of population sizes in each age group
#define model_rho (p[parindex[25]]) // Model code, determines which hypothesis is made for the age-specific reporting probabilities
#define model_q (p[parindex[26]]) // Model code, determines which hypothesis is made for the age-specific q
#define model_vac (p[parindex[27]]) // Model code, determines the hypothesis for vaccinal immunity (0: aP=wP; 1:wP=natural infection)
#define IC(J) (p[parindex[28] + (J)]) // Initial conditions

// Parameters on the transformed scale (all rates are per year)
#define Tiota (pt[parindex[0]]) // Background force of infection
#define Tq1(J) (pt[parindex[1] + (J)]) // Probability of infection given contact
#define Tq2(J) (pt[parindex[2] + (J)]) // Age-specific infectiousness
#define TomegaC(J) (pt[parindex[3] + (J)]) // Seasonality coefficients in 5-10 y
#define TomegaT(J) (pt[parindex[4] + (J)]) // Seasonality coefficients in 10-20 y
#define TCR(I, J) (pt[parindex[5] + (I + 15 * J)]) // Contact rates between age groups
#define Ttheta (pt[parindex[6]]) // Relative infectiousness of repeat vs. primary infections
#define Tbeta_sd (pt[parindex[7]]) // White noise intensity in infection rate
#define Tsigma (pt[parindex[8]]) // Inverse of latency period
#define Tgamma (pt[parindex[9]]) // Inverse of infectious period
#define TepsilonV (pt[parindex[10]]) // Leakiness of vaccine-induced immunity
#define TepsilonI (pt[parindex[11]]) // Leakiness of natural infection-induced immunity
#define TalphaV (pt[parindex[12]]) // Waning rate of vaccine-induced immunity
#define TalphaI (pt[parindex[13]]) // Waning rate of natural infection-induced immunity
#define Trho1(J) (pt[parindex[14] + (J)]) // Reporting probability for primary infections, age-specific
#define Trho2 (pt[parindex[15]]) // Reporting probability for repeat infections
#define Ttau (pt[parindex[16]]) // Reporting overdispersion
#define TepsA (pt[parindex[17]]) // Primary vaccine failure
#define Tv1 (pt[parindex[18]]) // Vaccine coverage for primary course
#define Tv2 (pt[parindex[19]]) // Vaccine coverage for preschool booster
#define Tt0 (pt[parindex[20]]) // Time at which vaccination is started
#define Tt1 (pt[parindex[21]]) // Time at which the maximal vaccine coverage is reached
#define Tt2 (pt[parindex[22]]) // Date of first data point for vaccine coverage
#define Tdelta(J) (pt[parindex[23] + (J)]) // Vector of aging rates in each age group
#define TN(J) (pt[parindex[24] + (J)]) // Vector of population sizes in each age group
#define Tmodel_rho (pt[parindex[25]]) // Model code, determines which hypothesis is made for the age-specific reporting probabilities
#define Tmodel_q (pt[parindex[26]]) // Model code, determines which hypothesis is made for the age-specific q
#define Tmodel_vac (pt[parindex[27]]) // Model code, determines the hypothesis for vaccinal immunity (0: aP=wP; 1:wP=natural infection)
#define TIC(J) (pt[parindex[28] + (J)]) // Initial conditions

// Covariates
#define birth (covars[covindex[0]])     // annual birth rate
#define mu(J) (covars[covindex[1] + (J)]) // Age-specific annual migration rates
#define seas(J) (covars[covindex[2] + (J)]) // seasonality basis functions

// Transforms the parameters to the transformed scale
void par_untrans (double *pt, double *p, int *parindex){
	int i, j;
	int *(*get_pomp_userdata_int)(const char *);
	get_pomp_userdata_int = (int *(*)(const char *)) R_GetCCallable("pomp","get_pomp_userdata_int");
	int nseas = *(get_pomp_userdata_int("nseas")); // No of periodic spline bases
	int nages = 17, nages_cmat = 15, nstates = 8; // No of age classes and states

	Tiota = log(iota);
	for(i = 0; i < nages; i++) Tq1(i) = logit(q1(i));
	for(i = 0; i < nages; i++) Tq2(i) = log(q2(i));
	for(i = 0; i < nseas; i++) {
		TomegaC(i) = omegaC(i);
		TomegaT(i) = omegaT(i);
	}
	for(i = 0; i < nages_cmat; i++) for(j = 0; j < nages_cmat; j++) TCR(i, j) = CR(i, j);
	Ttheta = logit(theta);
	Tbeta_sd = log(beta_sd);
	Tsigma = log(sigma);
	Tgamma = log(gamma);
	TepsilonV = logit(epsilonV);
	TepsilonI = logit(epsilonI);
	TalphaV = log(alphaV);
	TalphaI = log(alphaI);
	for(i = 0; i < nages; i++) Trho1(i) = logit(rho1(i));
	Trho2 = logit(rho2);
	Ttau = log(tau);
	TepsA = logit(epsA);
	Tv1 = logit(v1);
	Tv2 = logit(v2);
	Tt0 = t0;
	Tt2 = t2;
	Tt1 = logitCons(t1, t0, t2);
	for(i = 0; i < nages; i++) {
		Tdelta(i) = delta(i);
		TN(i) = N(i);
	}
	Tmodel_rho = model_rho;
	Tmodel_q = model_q;
	Tmodel_vac = model_vac;

	// Apply log-barycentric transform for each age group
	double sum[nages]; // Sums of fractions in every age group
	for(i = 0; i < nages; i++) sum[i] = 0.0;
	for(i = 0; i < nages; i++){
		for(j = 0; j < nstates; j++) {
			sum[i] += IC(i + nages * j);
		}
	}
	for(i = 0; i < nages; i++){
		for(j = 0; j < nstates; j++) {
			TIC(i + nages * j) = log(IC(i + nages * j) / sum[i]);
		}
	}
}

// Transforms the parameters to the natural scale
void par_trans (double *pt, double *p, int *parindex){
	int i, j;
	int *(*get_pomp_userdata_int)(const char *);
	get_pomp_userdata_int = (int *(*)(const char *)) R_GetCCallable("pomp","get_pomp_userdata_int");
	int nseas = *(get_pomp_userdata_int("nseas")); // No of periodic spline bases
	int nages = 17, nages_cmat = 15, nstates = 8; // No of age classes and states

	Tiota = exp(iota);
	for(i = 0; i < nages; i++) Tq1(i) = expit(q1(i));
	for(i = 0; i < nages; i++) Tq2(i) = exp(q2(i));
	for(i = 0; i < nseas; i++) {
		TomegaC(i) = omegaC(i);
		TomegaT(i) = omegaT(i);
	}
	for(i = 0; i < nages_cmat; i++) for(j = 0; j < nages_cmat; j++) TCR(i, j) = CR(i, j);
	Ttheta = expit(theta);
	Tbeta_sd = exp(beta_sd);
	Tsigma = exp(sigma);
	Tgamma = exp(gamma);
	TepsilonV = expit(epsilonV);
	TepsilonI = expit(epsilonI);
	TalphaV = exp(alphaV);
	TalphaI = exp(alphaI);
	for(i = 0; i < nages; i++) Trho1(i) = expit(rho1(i));
	Trho2 = expit(rho2);
	Ttau = exp(tau);
	TepsA = expit(epsA);
	Tv1 = expit(v1);
	Tv2 = expit(v2);
	Tt0 = t0;
	Tt2 = t2;
	Tt1 = expitCons(t1, t0, t2);
	for(i = 0; i < nages; i++) {
		Tdelta(i) = delta(i);
		TN(i) = N(i);
	}
	Tmodel_rho = model_rho;
	Tmodel_q = model_q;
	Tmodel_vac = model_vac;

	// Apply log-barycentric transform for each age group
	double sum[nages]; // Sums of fractions in every age group
	for(i = 0; i < nages; i++) sum[i] = 0.0;
	for(i = 0; i < nages; i++){
		for(j = 0; j < nstates; j++) {
			TIC(i + nages * j) = exp(IC(i + nages * j));
		}
	}
	for(i = 0; i < nages; i++){
		for(j = 0; j < nstates; j++) {
			sum[i] += exp(IC(i + nages * j));
		}
	}
	for(i = 0; i < nages; i++){
		for(j = 0; j < nstates; j++) {
			TIC(i + nages * j) /= sum[i];
		}
	}
}
//  Observation model: negative binomial, with mean rho1 * (C1 + rho2 * C2) and overdispersion parameter tau
// The alternative parametrization is used (dnbinom_mu), cf http://stat.ethz.ch/R-manual/R-patched/library/stats/html/NegBinomial.html
void dObs (double *lik, double *y, double *x, double *p, int give_log,
		int *obsindex, int *stateindex, int *parindex, int *covindex,
		int ncovars, double *covars, double t) {
	int nages_data = 7; // No of age classes in the data 0-1, 1-5, 5-10, 10-15, 15-20, 20-40, and 40+
	double rho1_vec[nages_data]; // Vector of reporting probability in each age class
	double C1_vec[nages_data], C2_vec[nages_data]; // No of cases in the aggregated age groups
	int i, j;
	double f = 0.0;

	// Check state variables
	for(i = 0; i < 10; i++){
		for(j = 0; j < 17; j++) {
			if(x[stateindex[i] + j] < 0.0 || ISNAN(x[stateindex[i] + j]) || ISNA(x[stateindex[i] + j])) x[stateindex[i] + j] = 0.0;
		}
	}

	// Sum no of primary cases
	C1_vec[0] = C1(0) + C1(1); // 0-1 age class = 0-4 mo + 4-12 mo
	C1_vec[1] = C1(2); // 1-5 y
	C1_vec[2] = C1(3); // 5-10 y
	C1_vec[3] = C1(4); // 10-15 y
	C1_vec[4] = C1(5); // 15-20 y
	C1_vec[5] = C1(6) + C1(7) + C1(8) + C1(9); // 20-40y = 20-25 + 25-30 + 30-35 + 35-40
	C1_vec[6] = C1(10) + C1(11) + C1(12) + C1(13) + C1(14) + C1(15) + C1(16); // 40+ = 40-45+45-50+50-55+55-60+60-65+65-70+70-75

	// Sum no of secondary cases
	C2_vec[0] = C2(0) + C2(1); // 0-1 age class = 0-4 mo + 4-12 mo
	C2_vec[1] = C2(2); // 1-5 y
	C2_vec[2] = C2(3); // 5-10 y
	C2_vec[3] = C2(4); // 10-15 y
	C2_vec[4] = C2(5); // 15-20 y
	C2_vec[5] = C2(6) + C2(7) + C2(8) + C2(9); // 20-40y = 20-25 + 25-30 + 30-35 + 35-40
	C2_vec[6] = C2(10) + C2(11) + C2(12) + C2(13) + C2(14) + C2(15) + C2(16); // 40+ = 40-45+45-50+50-55+55-60+60-65+65-70+70-75

	switch((int) model_rho) {

	// Model 51: 5 different reporting probabilities for 0-1, 1-5, 5-10, 10-20, and over 20
	// Reporting probability in 20+ assumed to be lower than that in 10-20
	case 51:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1) * rho1(0); // 1--5 y
		rho1_vec[2] = rho1(2) * rho1(0); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1(3); // 10-20
		rho1_vec[5] = rho1_vec[6] = rho1(4) * rho1(3); // over 20
		break;

		// Model 41: 4 different reporting probabilities for 0--1, 1--5, 5--10, and over 10, parametrization 1
	case 41:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1); // 1--5 y
		rho1_vec[2] = rho1(2); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(3); // over 10y
		break;
		// Model 42: as model 41, but different parametrization (reporting probability in 1--5 and 5--10 relative to 0--1)
	case 42:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1) * rho1(0); // 1--5 y
		rho1_vec[2] = rho1(2) * rho1(0); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(3); // over 10y
		break;
		// Model 31: 3 different reporting probabilities for 0-1, 1-10, 10+, parametrization 1
	case 31:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1); // 1--5 y
		rho1_vec[2] = rho1(1); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(2); // over 10y
		break;
		// Model 32: as model 31, but different paramterization (reporting probability in 1?10 relative to that in 0?1)
	case 32:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1) * rho1(0); // 1--5 y
		rho1_vec[2] = rho1(1) * rho1(0); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(2); // over 10y
		break;
		// Model 21: 2 different reporting probabilities for 0-10 and 10+
	case 21:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(0); // 1--5 y
		rho1_vec[2] = rho1(0); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(1); // over 10y
		break;
	default:
		error("Unknown model, %d\n", (int) model_rho);
		break;
	}

	// Sum the log-likelihoods for each age group
	for(i = 0; i < nages_data; i++) {
		if(ISNA(y[obsindex[i]])) {
			f += 0.0;
		} else {
			f += dnbinom_mu(nearbyint(y[obsindex[i]]), 1.0 / tau, rho1_vec[i] * (C1_vec[i] + rho2 * C2_vec[i]), 1);
		}
	}
	*lik = (give_log) ? f : exp(f);
}

void rObs (double *y, double *x, double *p,
		int *obsindex, int *stateindex, int *parindex, int *covindex,
		int ncovars, double *covars, double t) {
	int nages_data = 7; // No of age classes in the data 0-1, 1-5, 5-10, 10-15, 15-20, 20-40, and 40+
	double rho1_vec[nages_data]; // Vector of reporting probability in each age class
	double C1_vec[nages_data], C2_vec[nages_data]; // No of cases in the aggregated age groups
	int i;

	// Sum no of primary cases
	C1_vec[0] = C1(0) + C1(1); // 0-1 age class = 0-4 mo + 4-12 mo
	C1_vec[1] = C1(2); // 1-5 y
	C1_vec[2] = C1(3); // 5-10 y
	C1_vec[3] = C1(4); // 10-15 y
	C1_vec[4] = C1(5); // 15-20 y
	C1_vec[5] = C1(6) + C1(7) + C1(8) + C1(9); // 20-40y = 20-25 + 25-30 + 30-35 + 35-40
	C1_vec[6] = C1(10) + C1(11) + C1(12) + C1(13) + C1(14) + C1(15) + C1(16); // 40+ = 40-45+45-50+50-55+55-60+60-65+65-70+70-75

	// Sum no of secondary cases
	C2_vec[0] = C2(0) + C2(1); // 0-1 age class = 0-4 mo + 4-12 mo
	C2_vec[1] = C2(2); // 1-5 y
	C2_vec[2] = C2(3); // 5-10 y
	C2_vec[3] = C2(4); // 10-15 y
	C2_vec[4] = C2(5); // 15-20 y
	C2_vec[5] = C2(6) + C2(7) + C2(8) + C2(9); // 20-40y = 20-25 + 25-30 + 30-35 + 35-40
	C2_vec[6] = C2(10) + C2(11) + C2(12) + C2(13) + C2(14) + C2(15) + C2(16); // 40+ = 40-45+45-50+50-55+55-60+60-65+65-70+70-75

	switch((int) model_rho) {
	// Model 51: 5 different reporting probabilities for 0-1, 1-5, 5-10, 10-20, and over 20
	// Reporting probability in 20+ assumed to be lower than that in 10-20
	case 51:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1) * rho1(0); // 1--5 y
		rho1_vec[2] = rho1(2) * rho1(0); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1(3); // 10-20
		rho1_vec[5] = rho1_vec[6] = rho1(4) * rho1(3); // over 20
		break;

		// Model 41: 4 different reporting probabilities for 0--1, 1--5, 5--10, and over 10, parametrization 1
	case 41:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1); // 1--5 y
		rho1_vec[2] = rho1(2); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(3); // over 10y
		break;
		// Model 42: as model 41, but different parametrization (reporting probability in 1?5 and 5?10 relative to 0?1)
	case 42:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1) * rho1(0); // 1--5 y
		rho1_vec[2] = rho1(2) * rho1(0); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(3); // over 10y
		break;
		// Model 31: 3 different reporting probabilities for 0-1, 1-10, 10+, parametrization 1
	case 31:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1); // 1--5 y
		rho1_vec[2] = rho1(1); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(2); // over 10y
		break;
		// Model 32: as model 31, but different paramterization (reporting probability in 1?10 relative to that in 0?1)
	case 32:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(1) * rho1(0); // 1--5 y
		rho1_vec[2] = rho1(1) * rho1(0); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(2); // over 10y
		break;
		// Model 21: 2 different reporting probabilities for 0-10 and 10+
	case 21:
		rho1_vec[0] = rho1(0); // 0--4 mo, and 4--12 mo
		rho1_vec[1] = rho1(0); // 1--5 y
		rho1_vec[2] = rho1(0); // 5--10 y
		rho1_vec[3] = rho1_vec[4] = rho1_vec[5] = rho1_vec[6] = rho1(1); // over 10y
		break;
	default:
		error("Unknown model, %d\n", (int) model_rho);
		break;
	}
	for(i = 0; i < nages_data; i++) {
		y[obsindex[i]] = rnbinom_mu(1.0 / tau, rho1_vec[i] * (C1_vec[i] + rho2 * C2_vec[i]));
	}
}

// Process model, stochastic variant
void rSim (double *x, const double *p,
		const int *stateindex, const int *parindex, const int *covindex,
		int covdim, const double *covars,
		double t, double dt) {
	int i, j;
	int nages = 17, nages_cmat = 15; // No of age classes in the model, and in the contact matrix (without extra infant classes)
	double lambda[nages]; // Force of infection in each age group
	int *(*get_pomp_userdata_int)(const char *);
	get_pomp_userdata_int = (int *(*)(const char *)) R_GetCCallable("pomp","get_pomp_userdata_int");
	int nseas = *(get_pomp_userdata_int("nseas")); // No of periodic spline bases
	double N_vec[nages]; // Total population in each age group
	double time_switch = 1996.75 - 1990.0; // Time of switch to aP vaccine (October 1996)
	double tstart_boost = 1967.0 - 1990.0; // Time of start of preschool booster dose (1967)

	// Calculate vaccine coverage for primary vaccine course
	// Assume a linear ramp-up between t0 and t1, then constant from t1
	double v1_t = (t < t0) ? 0.0 : ((t < t1) ? (v1 * (t - t0) / (t1 - t0)) : v1);

	// Calculate vaccine coverage for second vaccine course (preschool booster)
	// if t3 > t1 (end of linear ramp-up), then constant coverage from time t3; otherwise assume linear ramp-up as for the primary course
	double v2_t, v2_t_corr;
	if(tstart_boost >= t1) {
		v2_t  = (t < tstart_boost) ? 0.0 : v2;
	} else {
		v2_t = (t < tstart_boost) ? 0.0 : ((t < t1) ? (v2 * (t - t0) / (t1 - t0)) : v2);
	}

	// Calculate corrected vaccine coverage for second course
	// This is necessary because all S1 are vaccinated in the model
	// That is, unvaccinated chldren and children vaccinated but for whom the vaccine didn't take
	// No such correction is necessary for S2, because children in this class were vaccinated previously
	v2_t_corr = epsA * v1_t / (epsA * v1_t + 1.0 - v1_t) * v2_t;

	// Check state variables
	for(i = 0; i < 10; i++){
		for(j = 0; j < 17; j++) {
			if(x[stateindex[i] + j] < 0.0 || ISNAN(x[stateindex[i] + j]) || ISNA(x[stateindex[i] + j])) x[stateindex[i] + j] = 0.0;
		}
	}

	for(i = 0; i < nages; i++) {
		N_vec[i] = V(i) + S1(i) + E1(i) + I1(i) + S2(i) + E2(i) + I2(i) + R(i);
	}

	// Create augmented contact matrix with seasonality factor
	// Necessary because more age classes in the model than in the contact data
	double Cseas[nages_cmat][nages_cmat]; // Matrix of contact rates with seasonal forcing
	for(i = 0; i < nages_cmat; i++) for(j = 0; j < nages_cmat; j++) Cseas[i][j] = CR(i, j);

	// Compute seasonality factors in 5--10y and 10--20 y
	double seas_multC = 0.0, seas_multT = 0.0;
	double sum_C = 0.0, sum_T = 0.0;
	for(i = 0; i < (nseas - 1); i++) {
		sum_C += omegaC(i);
		sum_T += omegaT(i);
		seas_multC += omegaC(i) * seas(i);
		seas_multT += omegaT(i) * seas(i);
	}
	seas_multC -= sum_C * seas(nseas - 1);
	seas_multT -= sum_T * seas(nseas - 1);
	seas_multC = exp(seas_multC);
	seas_multT = exp(seas_multT);

	Cseas[1][1] *= seas_multC; // 5--10 y
	Cseas[2][2] *= seas_multT; // 10--15 y
	Cseas[3][3] *= seas_multT; // 15--20 y

	double Cseas_aug[nages][nages_cmat]; // "Augmented" contact matrix with 2 extra rows for infant classes
	double I1_tilde[nages_cmat], I2_tilde[nages_cmat], N_tilde[nages_cmat]; // No infected and population size without infant age classes
	I1_tilde[0] = I1(0) + I1(1) + I1(2); // Sum no infected in 0-4 mo, 4-12 mo, and 1-5 y
	I2_tilde[0] = I2(0) + I2(1) + I2(2); // Sum no infected in 0-4 mo, 4-12 mo, and 1-5 y
	N_tilde[0] = N_vec[0] + N_vec[1] + N_vec[2]; // Sum population in 0-4 mo, 4-12 mo, and 1-5 y
	for(i = 1; i < nages_cmat; i++) {
		I1_tilde[i] = I1(i + 2);
		I2_tilde[i] = I2(i + 2);
		N_tilde[i] = N_vec[i + 2];
	}

	// Add two rows in the augmented contact matrix
	for(i = 0; i < nages; i++) {
		for(j = 0; j < nages_cmat; j++){
			if(i <= 2) {
				Cseas_aug[i][j] = Cseas[0][j];
			} else {
				Cseas_aug[i][j] = Cseas[i - 2][j];
			}
		}
	}

	// Define age-specific susceptibility and infectiousness, depending on the model
	double q1_vec[nages], q2_vec[nages_cmat];
	switch((int) model_q) {

	// Model 1: 1 susceptibility factor and infectiousness (equal to 1) for all age groups
	case 1:
		for(i = 0; i < nages; i++) q1_vec[i] = q1(0);
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(0);
		break;
		// Model 34: age-independent infectiousness, 4 age-specific susceptibility factors (0-1, 1-5, 5-10, 10+)
	case 34:
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(0);
		q1_vec[0] = q1_vec[1] = q1(0); // 0-4 mo and 4-12?mo
		q1_vec[2] = q1(1); // 1-5 y
		q1_vec[3] = q1(2); // 5-10 y
		for(i = 4; i < nages; i++) q1_vec[i] = q1(3); // 10+
		break;
		// Model 33: age-independent infectiousness, 3 age-specific susceptibility factors (0-1, 1-10, 10+)
	case 33:
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(0);
		q1_vec[0] = q1_vec[1] = q1(0); // 0-4 mo and 4-12?mo
		q1_vec[2] = q1(1); // 1-5 y
		q1_vec[3] = q1(1); // 5-10 y
		for(i = 4; i < nages; i++) q1_vec[i] = q1(2); // 10+
		break;
		// Model 332: age-independent infectiousness, 3 age-specific susceptibility factors (0-10, 10-20, 20+)
		// Susceptibility is assumed to decrease with age
	case 332:
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(i);
		q1_vec[0] = q1_vec[1] = q1(0); // 0-4 mo and 4-12?mo
		q1_vec[2] = q1(0); // 1-5 y
		q1_vec[3] = q1(0); // 5-10 y
		q1_vec[4] = q1_vec[5] = q1(0) * q1(1); // 10-20
		for(i = 6; i < nages; i++) q1_vec[i] = q1(0) * q1(1) * q1(2); // 20+
		break;
		// Model 32: age-independent infectiousness, 2 age-specific susceptibilities (0-10, 10+)
	case 32:
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(0);
		q1_vec[0] = q1_vec[1] = q1(0); // 0-4 mo and 4-12?mo
		q1_vec[2] = q1(0); // 1-5 y
		q1_vec[3] = q1(0); // 5-10 y
		for(i = 4; i < nages; i++) q1_vec[i] = q1(1); // 10+
		break;
	default:
		error("Unknown model, %d\n", (int) model_q);
		break;
	}

	// Force of infection in each age group
	for(i = 0; i < nages; i++) {
		lambda[i] = 0.0; // Initialize values
		for(j = 0; j < nages_cmat; j++){
			lambda[i] += q2_vec[j] * Cseas_aug[i][j] * (I1_tilde[j] + theta * I2_tilde[j] + iota) / N_tilde[j];
		}
		lambda[i] *= q1_vec[i];
	}

	// Transitions between compartments
	double birth_no;
	double fromV[nages][3], fromS1[nages][3], fromE1[nages][3], fromI1[nages][3];
	double fromS2[nages][3], fromE2[nages][3], fromI2[nages][3], fromR[nages][3];
	double rate2[2], trans2[2];
	double rate3[3], trans3[3];
	double rRE2[nages], rVE2[nages];
	for(i = 0; i < nages; i++) {
		rRE2[i] = epsilonI * lambda[i]; // Rate of transition from R to E2
		rVE2[i] = epsilonV * lambda[i]; // Rate of transition from V to E2
	}

	// Initialize nof transitions
	for(i = 0; i < nages; i++) {
		for(j = 0; j < 3; j++){
			fromV[i][j] = 0.0;
			fromS1[i][j] = 0.0;
			fromE1[i][j] = 0.0;
			fromI1[i][j] = 0.0;
			fromS2[i][j] = 0.0;
			fromE2[i][j] = 0.0;
			fromI2[i][j] = 0.0;
			fromR[i][j] = 0.0;
		}
	}

	birth_no = rpois(birth * dt); // Births

	// Transitions in newborns, i=0
	// From S1
	rate3[0] = lambda[0]; // Infection
	rate3[1] = delta(0) * (1.0 - v1_t * (1.0 - epsA)); // Aging and not vaccinated
	rate3[2] = delta(0) * v1_t * (1.0 - epsA); // Aging and vaccinated
	reulermultinom(3, S1(0), &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) fromS1[0][j] = trans3[j];

	// From S2
	reulermultinom(3, S2(0), &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) fromS2[0][j] = trans3[j];

	// From E1
	rate2[0] = sigma; // Progression to I
	rate2[1] = delta(0); // Aging
	reulermultinom(2, E1(0), &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) fromE1[0][j] = trans2[j];

	// From E2
	reulermultinom(2, E2(0), &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) fromE2[0][j] = trans2[j];

	// From I1
	rate2[0] = gamma; // Recovery
	rate2[1] = delta(0); // Aging
	reulermultinom(2, I1(0), &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) fromI1[0][j] = trans2[j];

	// From I2
	reulermultinom(2, I2(0), &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) fromI2[0][j] = trans2[j];

	// From R
	rate3[0] = alphaI; // Waning of infection-derived immunity
	rate3[1] = rRE2[0]; // Leak of infection-derived immunity
	rate3[2] = delta(0); // Aging and not vaccinated
	reulermultinom(3, R(0), &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) fromR[0][j] = trans3[j];

	// Transitions in 1-5y, i=2
	// From V
	rate3[0] = alphaV; // Waning of vaccine-derived immunity
	rate3[1] = rVE2[2]; // Leak of vaccine-derived immunity
	rate3[2] = delta(2); // Aging
	reulermultinom(3, V(2), &rate3[0], dt, &trans3[0]);
	fromV[2][0] = trans3[0]; // From V to S2
	fromV[2][1] = trans3[1]; // From V to E2
	fromV[2][2] = trans3[2]; // From V_i to V_{i+1}

	// From S1
	rate3[0] = lambda[2]; // Infection
	rate3[1] = delta(2) * (1.0 - v2_t_corr * (1.0 - epsA)); // Aging and not vaccinated
	rate3[2] = delta(2) * v2_t_corr * (1.0 - epsA); // Aging and vaccinated
	reulermultinom(3, S1(2), &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) fromS1[2][j] = trans3[j];

	// From S2
	rate3[0] = lambda[2]; // Infection
	rate3[1] = delta(2) * (1.0 - v2_t * (1.0 - epsA)); // Aging and not vaccinated
	rate3[2] = delta(2) * v2_t * (1.0 - epsA); // Aging and vaccinated
	reulermultinom(3, S2(2), &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) fromS2[2][j] = trans3[j];

	// From E1
	rate2[0] = sigma; // Progression to I
	rate2[1] = delta(2); // Aging
	reulermultinom(2, E1(2), &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) fromE1[2][j] = trans2[j];

	// From E2
	reulermultinom(2, E2(2), &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) fromE2[2][j] = trans2[j];

	// From I1
	rate2[0] = gamma; // Recovery
	rate2[1] = delta(2); // Aging
	reulermultinom(2, I1(2), &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) fromI1[2][j] = trans2[j];

	// From I2
	reulermultinom(2, I2(2), &rate2[0], dt, &trans2[0]);
	for(j = 0; j < 2; j++) fromI2[2][j] = trans2[j];

	// From R
	rate3[0] = alphaI; // Waning of infection-derived immunity
	rate3[1] = rRE2[2]; // Leak of infection-derived immunity
	rate3[2] = delta(2); // Aging and not vaccinated
	reulermultinom(3, R(2), &rate3[0], dt, &trans3[0]);
	for(j = 0; j < 3; j++) fromR[2][j] = trans3[j];

	// Transitions in other age groups, i=1, i>=3
	for(i = 1; i < nages; i++){

		if(i != 2) {
			// From V
			rate3[0] = alphaV; // Waning of vaccine-derived immunity
			rate3[1] = rVE2[i]; // Leak of vaccine-derived immunity
			rate3[2] = delta(i); // Aging
			reulermultinom(3, V(i), &rate3[0], dt, &trans3[0]);
			fromV[i][0] = trans3[0]; // From V to S2
			fromV[i][1] = trans3[1]; // From V to E2
			fromV[i][2] = trans3[2]; // From V_i to V_{i+1}

			// From S1
			rate2[0] = lambda[i]; // Primary infection
			rate2[1] = delta(i); // Aging
			reulermultinom(2, S1(i), &rate2[0], dt, &trans2[0]);
			fromS1[i][0] = trans2[0]; // from S1 to E1
			fromS1[i][1] = trans2[1]; // From S1_i to S1_{i+1}

			// From S2
			reulermultinom(2, S2(i), &rate2[0], dt, &trans2[0]);
			fromS2[i][0] = trans2[0]; // from S2 to E2
			fromS2[i][1] = trans2[1]; // From S2_i to S2_{i+1}

			// From E1
			rate2[0] = sigma; // Progression to I
			rate2[1] = delta(i); // Aging
			reulermultinom(2, E1(i), &rate2[0], dt, &trans2[0]);
			fromE1[i][0] = trans2[0]; // from E1 to I1
			fromE1[i][1] = trans2[1]; // From E1_i to E1_{i+1}

			// From E2
			reulermultinom(2, E2(i), &rate2[0], dt, &trans2[0]);
			fromE2[i][0] = trans2[0]; // from E2 to I2
			fromE2[i][1] = trans2[1]; // From E2_i to E2_{i+1}

			// From I1
			rate2[0] = gamma; // Recovery
			rate2[1] = delta(i); // Aging
			reulermultinom(2, I1(i), &rate2[0], dt, &trans2[0]);
			fromI1[i][0] = trans2[0]; // from I1 to R
			fromI1[i][1] = trans2[1]; // From I1_i to I1_{i+1}

			// From I2
			reulermultinom(2, I2(i), &rate2[0], dt, &trans2[0]);
			fromI2[i][0] = trans2[0]; // from I2 to R
			fromI2[i][1] = trans2[1]; // From I2_i to I2_{i+1}

			// From R
			rate3[0] = alphaI; // Waning of infection-derived immunity
			rate3[1] = rRE2[i]; // Leak of infection-derived immunity
			rate3[2] = delta(i); // Aging
			reulermultinom(3, R(i), &rate3[0], dt, &trans3[0]);
			fromR[i][0] = trans3[0]; // From R to S2
			fromR[i][1] = trans3[1]; // From R to E2
			fromR[i][2] = trans3[2]; // From R_i to R_{i+1}
		}
	}

	// Balance the equations
	// Age class 0, newborns
	V(0) += 0.0; // Vaccinated
	S1(0) += birth_no - fromS1[0][0] - fromS1[0][1] - fromS1[0][2]; // Primary susceptible
	E1(0) += fromS1[0][0] - fromE1[0][0] - fromE1[0][1]; // Primary exposed
	I1(0) += fromE1[0][0] - fromI1[0][0] - fromI1[0][1]; // Primary infected
	S2(0) += fromR[0][0] - fromS2[0][0] - fromS2[0][1] - fromS2[0][2]; // Secondary susceptible
	E2(0) += fromS2[0][0] + fromR[0][1] - fromE2[0][0] - fromE2[0][1]; // Secondary exposed
	I2(0) += fromE2[0][0] - fromI2[0][0] - fromI2[0][1]; // Secondary infected
	R(0) += fromI1[0][0] + fromI2[0][0] - fromR[0][0] - fromR[0][1] - fromR[0][2]; // Recovered
	C1(0) += fromI1[0][0]; // Cumulative primary infections
	C2(0) += fromI2[0][0]; // Cumulative secondary infections

	// older age groups, i>0
	for(i = 1; i < nages; i++){
		V(i) += fromV[i-1][2] - fromV[i][0] - fromV[i][1] - fromV[i][2]; // Vaccinated
		S1(i) += fromS1[i-1][1] - fromS1[i][0] - fromS1[i][1] - fromS1[i][2]; // Primary susceptible
		E1(i) += fromE1[i-1][1] + fromS1[i][0] - fromE1[i][0] - fromE1[i][1]; // Primary exposed
		I1(i) += fromI1[i-1][1] + fromE1[i][0] - fromI1[i][0] - fromI1[i][1]; // Primary infected
		S2(i) += fromS2[i-1][1] + fromR[i][0] + fromV[i][0] - fromS2[i][0] - fromS2[i][1] - fromS2[i][2]; // Secondary susceptible
		E2(i) += fromE2[i-1][1] + fromS2[i][0] + fromR[i][1] + fromV[i][1] - fromE2[i][0] - fromE2[i][1]; // Secondary exposed
		I2(i) += fromI2[i-1][1] + fromE2[i][0] - fromI2[i][0] - fromI2[i][1]; // Secondary infected
		R(i) += fromR[i-1][2] +  fromI1[i][0] + fromI2[i][0] - fromR[i][0] - fromR[i][1] - fromR[i][2]; // Recovered
		C1(i) += fromI1[i][0]; // Cumulative primary infections
		C2(i) += fromI2[i][0]; // Cumulative secondary infections

		if((int) model_vac == 0) {
			V(i) += fromS1[i-1][2] + fromS2[i-1][2]; // wP=aP, all S move to V
		} else {
			if(t < time_switch) {
				R(i) += fromS1[i-1][2] + fromS2[i-1][2]; // wP era: vaccinated move from S to R
			} else {
				V(i) += fromS1[i-1][2] + fromS2[i-1][2]; // aP era: vaccinated move from S to V
			}
		}
	}

	double mig[nages]; // Variations due to migrationsin each age group
	for(i = 0; i < nages; i++) mig[i] = exp(mu(i) * dt);
	for(i = 0; i < nages; i++) {
		V(i) = nearbyint(V(i) * mig[i]);
		S1(i) = nearbyint(S1(i) * mig[i]);
		E1(i) = nearbyint(E1(i) * mig[i]);
		I1(i) = nearbyint(I1(i) * mig[i]);
		S2(i) = nearbyint(S2(i) * mig[i]);
		E2(i) = nearbyint(E2(i) * mig[i]);
		I2(i) = nearbyint(I2(i) * mig[i]);
		R(i) = nearbyint(R(i) * mig[i]);
		C1(i) = nearbyint(C1(i));
		C2(i) = nearbyint(C2(i));
	}
}

// Continuous-time deterministic skeleton
void skel_continuous (double *f, double *x, double *p,
		int *stateindex, int *parindex, int *covindex,
		int ncovars, double *covars, double t) {

	int i, j;
	int nages = 17, nages_cmat = 15; // No of age classes in the model, and in the contact matrix (without infant classes)
	double lambda[nages]; // Force of infection in each age group
	int *(*get_pomp_userdata_int)(const char *);
	get_pomp_userdata_int = (int *(*)(const char *)) R_GetCCallable("pomp","get_pomp_userdata_int");
	int nseas = *(get_pomp_userdata_int("nseas")); // No of periodic spline bases
	double time_switch = 1996.75 - 1990.0; // Time of switch to aP vaccine (October 1996)
	double tstart_boost = 1967.0 - 1990.0; // Time of start of preschool booster dose (1967)

	// Calculate vaccine coverage for primary vaccine course
	// Assume a linear ramp-up between t0 and t1, then constant from t1
	double v1_t = (t < t0) ? 0.0 : ((t < t1) ? (v1 * (t - t0) / (t1 - t0)) : v1);

	// Calculate vaccine coverage for second vaccine course (preschool booster)
	// if t3 > t1 (end of linear ramp-up), then constant coverage from time t3; otherwise assume linear ramp-up as for the primary course
	double v2_t, v2_t_corr;
	if(tstart_boost >= t1) {
		v2_t  = (t < tstart_boost) ? 0.0 : v2;
	} else {
		v2_t = (t < tstart_boost) ? 0.0 : ((t < t1) ? (v2 * (t - t0) / (t1 - t0)) : v2);
	}

	// Calculate corrected vaccine coverage for second course
	// This is necessary because all S1 are vaccinated in the model
	// That is, unvaccinated chldren and children vaccinated but for whom the vaccine didn't take
	// No such correction is necessary for S2, because children in this class were vaccinated previously
	v2_t_corr = epsA * v1_t / (epsA * v1_t + 1.0 - v1_t) * v2_t;

	// Check state variables
	for(i = 0; i < 10; i++){
		for(j = 0; j < 17; j++) {
			if(x[stateindex[i] + j] < 0.0 || ISNAN(x[stateindex[i] + j]) || ISNA(x[stateindex[i] + j])) x[stateindex[i] + j] = 0.0;
		}
	}

	// Total population in each group
	double N_vec[nages]; // Total population in each age group
	for(i = 0; i < nages; i++) {
		N_vec[i] = V(i) + S1(i) + E1(i) + I1(i) + S2(i) + E2(i) + I2(i) + R(i);
	}

	// Create augmented contact matrix with seasonality factor
	// Necessary because more age classes in the model than in the contact data
	double Cseas[nages_cmat][nages_cmat]; // Matrix of contact rates with seasonal forcing
	for(i = 0; i < nages_cmat; i++) for(j = 0; j < nages_cmat; j++) Cseas[i][j] = CR(i, j);

	// Compute seasonality factors in 5--10y and 10--20 y
	double seas_multC = 0.0, seas_multT = 0.0;
	double sum_C = 0.0, sum_T = 0.0;
	for(i = 0; i < (nseas - 1); i++) {
		sum_C += omegaC(i);
		sum_T += omegaT(i);
		seas_multC += omegaC(i) * seas(i);
		seas_multT += omegaT(i) * seas(i);
	}
	seas_multC -= sum_C * seas(nseas - 1);
	seas_multT -= sum_T * seas(nseas - 1);
	seas_multC = exp(seas_multC);
	seas_multT = exp(seas_multT);

	Cseas[1][1] *= seas_multC; // 5--10 y
	Cseas[2][2] *= seas_multT; // 10--15 y
	Cseas[3][3] *= seas_multT; // 15--20 y

	double Cseas_aug[nages][nages_cmat]; // "Augmented" contact matrix with 2 extra rows for infant classes
	double I1_tilde[nages_cmat], I2_tilde[nages_cmat], N_tilde[nages_cmat]; // No infected and population size without infant age classes
	I1_tilde[0] = I1(0) + I1(1) + I1(2); // Sum no infected in 0-4 mo, 4-12 mo, and 1-5 y
	I2_tilde[0] = I2(0) + I2(1) + I2(2); // Sum no infected in 0-4 mo, 4-12 mo, and 1-5 y
	N_tilde[0] = N_vec[0] + N_vec[1] + N_vec[2]; // Sum population in 0-4 mo, 4-12 mo, and 1-5 y
	for(i = 1; i < nages_cmat; i++) {
		I1_tilde[i] = I1(i + 2);
		I2_tilde[i] = I2(i + 2);
		N_tilde[i] = N_vec[i + 2];
	}

	// Add two rows in the augmented contact matrix
	for(i = 0; i < nages; i++) {
		for(j = 0; j < nages_cmat; j++){
			if(i <= 2) {
				Cseas_aug[i][j] = Cseas[0][j];
			} else {
				Cseas_aug[i][j] = Cseas[i - 2][j];
			}
		}
	}

	// Define age-specific susceptibility and infectiousness, depending on the model
	double q1_vec[nages], q2_vec[nages_cmat];
	switch((int) model_q) {

	// Model 1: 1 susceptibility factor and infectiousness (equal to 1) for all age groups
	case 1:
		for(i = 0; i < nages; i++) q1_vec[i] = q1(0);
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(0);
		break;
		// Model 23: Same susceptibility factor for all age groups; 3 age-specific infectiousness (0-5y, 5-10y, and 10+)
	case 23:
		for(i = 0; i < nages; i++) q1_vec[i] = q1(0);
		q2_vec[0] = q2(0); // 0--5y
		q2_vec[1] = q2(1); // 5--10y
		for(i = 2; i < nages_cmat; i++) q2_vec[i] = q2(2);
		break;
		// Model 34: age-independent infectiousness, 4 age-specific susceptibility factors (0-1, 1-5, 5-10, 10+)
	case 34:
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(0);
		q1_vec[0] = q1_vec[1] = q1(0); // 0-4 mo and 4-12?mo
		q1_vec[2] = q1(1); // 1-5 y
		q1_vec[3] = q1(2); // 5-10 y
		for(i = 4; i < nages; i++) q1_vec[i] = q1(3); // 10+
		break;
		// Model 33: age-independent infectiousness, 3 age-specific susceptibility factors (0-1, 1-10, 10+)
	case 33:
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(0);
		q1_vec[0] = q1_vec[1] = q1(0); // 0-4 mo and 4-12?mo
		q1_vec[2] = q1(1); // 1-5 y
		q1_vec[3] = q1(1); // 5-10 y
		for(i = 4; i < nages; i++) q1_vec[i] = q1(2); // 10+
		break;
		// Model 332: age-independent infectiousness, 3 age-specific susceptibility factors (0-10, 10-20, 20+)
		// Susceptibility is assumed to decrease with age
	case 332:
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(i);
		q1_vec[0] = q1_vec[1] = q1(0); // 0-4 mo and 4-12?mo
		q1_vec[2] = q1(0); // 1-5 y
		q1_vec[3] = q1(0); // 5-10 y
		q1_vec[4] = q1_vec[5] = q1(0) * q1(1); // 10-20
		for(i = 6; i < nages; i++) q1_vec[i] = q1(0) * q1(1) * q1(2); // 20+
		break;
		// Model 32: age-independent infectiousness, 2 age-specific susceptibilities (0-10, 10+)
	case 32:
		for(i = 0; i < nages_cmat; i++) q2_vec[i] = q2(0);
		q1_vec[0] = q1_vec[1] = q1(0); // 0-4 mo and 4-12?mo
		q1_vec[2] = q1(0); // 1-5 y
		q1_vec[3] = q1(0); // 5-10 y
		for(i = 4; i < nages; i++) q1_vec[i] = q1(1); // 10+
		break;
	default:
		error("Unknown model, %d\n", (int) model_q);
		break;
	}

	// Force of infection in each age group
	for(i = 0; i < nages; i++) {
		lambda[i] = 0.0; // Initialize values
		for(j = 0; j < nages_cmat; j++){
			lambda[i] += q2_vec[j] * Cseas_aug[i][j] * (I1_tilde[j] + theta * I2_tilde[j] + iota) / N_tilde[j];
		}
		lambda[i] *= q1_vec[i];
	}

	// Derivatives

	// Model equations in newborns (0-4 mo, index 0)
	DV(0) = 0.0; // Vaccinated
	DS1(0) = birth - (lambda[0] + delta(0) - mu(0)) * S1(0); // Primary susceptible
	DE1(0) = lambda[0] * S1(0) - (sigma + delta(0) - mu(0)) * E1(0); // Primary exposed
	DI1(0) = sigma * E1(0) - (gamma + delta(0) - mu(0)) * I1(0); // Primary infected
	DS2(0) = alphaV * V(0) + alphaI * R(0) - (lambda[0] + delta(0) - mu(0)) * S2(0); // Secondary susceptible
	DE2(0) = lambda[0] * (S2(0) + epsilonV * V(0) + epsilonI * R(0)) - (sigma + delta(0) - mu(0)) * E2(0); // Secondary exposed
	DI2(0) = sigma * E2(0) - (gamma + delta(0) - mu(0)) * I2(0); // Secondary infected
	DR(0) = gamma * (I1(0) + I2(0)) - (alphaI + epsilonI * lambda[0] + delta(0) - mu(0)) * R(0); // Recovered
	DC1(0) = gamma * I1(0); // Cumulative primary infections
	DC2(0) = gamma * I2(0); // Cumulative repeat infections

	// Model equations in older age groups (i>0)
	for(i = 1; i < nages; i++) {
		DV(i) = delta(i - 1) * V(i - 1) - (alphaV + epsilonV * lambda[i] + delta(i) - mu(i)) * V(i); // Vaccinated
		DS1(i) = delta(i - 1) * S1(i - 1) - (lambda[i] + delta(i) - mu(i)) * S1(i); // Primary susceptible
		DE1(i) = delta(i - 1) * E1(i - 1) + lambda[i] * S1(i) - (sigma + delta(i) - mu(i)) * E1(i); // Primary exposed
		DI1(i) = delta(i - 1) * I1(i - 1) + sigma * E1(i) - (gamma + delta(i) - mu(i)) * I1(i); // Primary infected
		DS2(i) = delta(i - 1) * S2(i - 1) + alphaV * V(i) + alphaI * R(i) - (lambda[i] + delta(i) - mu(i)) * S2(i); // Secondary susceptible
		DE2(i) = delta(i - 1) * E2(i - 1) + lambda[i] * (S2(i) + epsilonV * V(i) + epsilonI * R(i)) - (sigma + delta(i) - mu(i)) * E2(i); // Secondary exposed
		DI2(i) = delta(i - 1) * I2(i - 1) + sigma * E2(i) - (gamma + delta(i) - mu(i)) * I2(i); // Secondary infected
		DR(i) = delta(i - 1) * R(i - 1) + gamma * (I1(i) + I2(i)) - (alphaI + epsilonI * lambda[i] + delta(i) - mu(i)) * R(i); // Recovered
		DC1(i) = gamma * I1(i); // Cumulative primary infections
		DC2(i) = gamma * I2(i); // Cumulative repeat infections
	}

	// Model vaccination: transfer of susceptibles from S to V or R class upon aging
	if((int) model_vac == 0){
		// aP=wP: susceptibles move to V
		DV(1) += v1_t * (1.0 - epsA) * delta(0) * (S1(0) + S2(0)); // 0-4 to 4-12 mo, primary course
		DS1(1) -= v1_t * (1.0 - epsA) * delta(0) * S1(0);
		DS2(1) -= v1_t * (1.0 - epsA) * delta(0) * S2(0);
		DV(3) += (1.0 - epsA) * delta(2) * (v2_t_corr * S1(2) + v2_t * S2(2)); // 1-5 y to 5-10, preschool booster
		DS1(3) -= v2_t_corr * (1.0 - epsA) * delta(2) * S1(2);
		DS2(3) -= v2_t * (1.0 - epsA) * delta(2) * S2(2);
	} else {
		// wP=natural infection: suceptibles move to R before the switch to aP, to V after
		if(t < time_switch) {
			DR(1) += v1_t * (1.0 - epsA) * delta(0) * (S1(0) + S2(0)); // 0-4 to 4-12 mo, primary course
			DS1(1) -= v1_t * (1.0 - epsA) * delta(0) * S1(0);
			DS2(1) -= v1_t * (1.0 - epsA) * delta(0) * S2(0);
			DR(3) += (1.0 - epsA) * delta(2) * (v2_t_corr * S1(2) + v2_t * S2(2)); // 1-5 y to 5-10, preschool booster
			DS1(3) -= v2_t_corr * (1.0 - epsA) * delta(2) * S1(2);
			DS2(3) -= v2_t * (1.0 - epsA) * delta(2) * S2(2);
		} else {
			DV(1) += v1_t * (1.0 - epsA) * delta(0) * (S1(0) + S2(0)); // 0-4 to 4-12 mo, primary course
			DS1(1) -= v1_t * (1.0 - epsA) * delta(0) * S1(0);
			DS2(1) -= v1_t * (1.0 - epsA) * delta(0) * S2(0);
			DV(3) += (1.0 - epsA) * delta(2) * (v2_t_corr * S1(2) + v2_t * S2(2)); // 1-5 y to 5-10, preschool booster
			DS1(3) -= v2_t_corr * (1.0 - epsA) * delta(2) * S1(2);
			DS2(3) -= v2_t * (1.0 - epsA) * delta(2) * S2(2);
		}
	}
}

#undef V
#undef S1
#undef E1
#undef I1
#undef S2
#undef E2
#undef I2
#undef R
#undef C1
#undef C2

#undef DV
#undef DS1
#undef DE1
#undef DI1
#undef DS2
#undef DE2
#undef DI2
#undef DR
#undef DC1
#undef DC2

#undef iota
#undef q1
#undef q2
#undef omegaC
#undef omegaT
#undef CR
#undef theta
#undef beta_sd
#undef sigma
#undef gamma
#undef epsilonV
#undef epsilonI
#undef alphaV
#undef alphaI
#undef rho1
#undef rho2
#undef tau
#undef epsA
#undef v1
#undef v2
#undef t0
#undef t1
#undef t2
#undef delta
#undef N
#undef model_rho
#undef model_q
#undef model_vac
#undef IC

#undef Tiota
#undef Tq1
#undef Tq2
#undef TomegaC
#undef TomegaT
#undef TCR
#undef Ttheta
#undef Tbeta_sd
#undef Tsigma
#undef Tgamma
#undef TepsilonV
#undef TepsilonI
#undef TalphaV
#undef TalphaI
#undef Trho1
#undef Trho2
#undef Ttau
#undef TepsA
#undef Tv1
#undef Tv2
#undef Tt0
#undef Tt1
#undef Tt2
#undef Tdelta
#undef TN
#undef Tmodel_rho
#undef Tmodel_q
#undef tmodel_vac
#undef TIC

#undef birth
#undef mu
#undef seas
