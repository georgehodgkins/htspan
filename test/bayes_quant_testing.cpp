#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <fstream>

#include "htspan/bayes_orient_bias_quant.hpp"
#include "htspan/nucleotide.hpp"
#include "frontend/simul_writer.hpp"

#include "stograd/src/stograd/stograd.hpp"

#define NEPOCHS 3000
#define NLOCS 75000
#define EPS 1e-6

#define RD_MIN 50
#define RD_MAX 150

using namespace std;

template <typename bayes_stepper_t>
void core_testing_method (hts::frontend::simul_writer &output, const char* test_name, const char* stepper_name, const size_t bsize,
		const double lrate, const double beta_th, const double beta_ph, const double alpha0, const double beta0, int sim_seed) {

	hts::bayes_orient_bias_quant_f bobquant(nuc_G, nuc_T);
	bobquant.simulate(NLOCS, 1, beta_th, 1, beta_ph, RD_MIN, RD_MAX, sim_seed);

	output << test_name << "\t" << stepper_name << "\t" << bsize << "\t" << lrate << 
		"\t" << alpha0 << "\t" << beta0 << "\t...";

	time_t start, finish;
	try {

		hts::hparams result;
		start = time(NULL);

		if (alpha0 < 0) { // use freq estimate for initial alpha and beta
			result = bobquant.operator()<bayes_stepper_t>(bsize, NEPOCHS, lrate, EPS);
		} else { // use fixed estimate
			result = bobquant.operator()<bayes_stepper_t>(bsize, NEPOCHS, lrate, EPS,
				alpha0, beta0, alpha0, beta0);
		}

		finish = time(NULL);

		int elapsed = finish - start;
		output << "\b\b\b" << bobquant.n_theta_epochs << '\t' << bobquant.n_phi_epochs << '\t' << result.alpha_theta <<
			'\t' << result.beta_theta << '\t' << result.alpha_phi << '\t' << result.beta_phi << '\t' <<
			elapsed / 60 << ':' << elapsed % 60 << '\n';

	} catch (exception &e) {
		output << "\b\b\bException thrown: " << e.what() << '\n';
	}
}


int main(int argc, char** argv) {

	if (argc < 2) {
		throw runtime_error("Seed argument reqd!");
	}
	int sim_seed = atoi(argv[1]);
	if (sim_seed == 0) {
		throw runtime_error("Seed cannot be zero or invalid!");
	}

	// optionally write to both stdout and a file
	hts::frontend::simul_writer output;
	output.use_cout(true);
	output.setprecision(5);//set to five digits of output, no scientific notation
	// any value here that needs to be reported in scientific notation is oob anyway
	cout.setf(ios_base::fixed);
	if (argc > 2) {
		output.add_file(argv[2]);
	}

	output << "Testing various configurations of stograd for Bayesian quant.\n" <<
		"Fixed params:\nEpoch count: " << NEPOCHS << "\nLocus count: " << NLOCS <<
		"\nEps for convergence: " << EPS << "\nRead count (min, max): " << RD_MIN << ", " << RD_MAX << "\n" <<
		"Random seed: " << sim_seed << "\n" <<
		"Parameter values for tests, fmt (alpha, beta): high=(1,100), low=(1,1000), very low=(1,10000)\n\n" <<

		"Test results:\n" <<
		"test\tstepper\tbsize\tlrate\talpha0\tbeta0\tth_ep\tph_ep\ta_th\tb_th\ta_ph\tb_ph\ttime\n";

	// set of batch sizes to test, indexed by n_b
	const size_t BSIZES[] = {100, 1000, 10000};
	// set of learning rates to test, indexed by n_lr
	const double LRATES[] = {1e-2, 1e-3, 1e-4}; 

	// set of initial alpha values to test (-1 = use freq to estimate), indexed by n_ips
	const double IALPHAS[] = {-1, 1, 1, .01};
	// set of initial beta values to test, indexed by n_ips
	const double IBETAS[] = {-1, 1000, 10000, 100000, 100000}; 

	// set of test case names, indexed by n_tst
	const char* TESTS[] = {"ls_hd", "ls_ld", "ls_vld", "hs_ld"};
	// set of beta_theta values for test cases, params for simulation, indexed by n_tst
	const double BETA_THETAS[] = {1000, 1000, 1000, 100}; 
	// set of beta_phi values for test cases, params for simulation, indexed by n_tst
	const double BETA_PHIS[] = {100, 1000, 10000, 1000}; 

	// select batch size
	for (size_t n_b = 0; n_b < sizeof(BSIZES)/sizeof(int); ++n_b) {
		// select test
		for (size_t n_tst = 0; n_tst < sizeof(TESTS)/sizeof(char*); ++n_tst) {
			// select initial parameter
			for (size_t n_ips = 0; n_ips < sizeof(IALPHAS)/sizeof(double); ++n_ips) {
				// yamadam does not use a learning rate
				//core_testing_method<stograd::stepper::yamadam<double> > (output, TESTS[n_tst], "yamadam", BSIZES[n_b], 0,
				//	BETA_THETAS[n_tst], BETA_PHIS[n_tst], IALPHAS[n_ips], IBETAS[n_ips], sim_seed);

				// select learning rate
				for (size_t n_lr = 0; n_lr < sizeof(LRATES)/sizeof(double); ++n_lr) {
					core_testing_method<stograd::stepper::constant<double> > (output, TESTS[n_tst], "const", BSIZES[n_b], LRATES[n_lr],
						BETA_THETAS[n_tst], BETA_PHIS[n_tst], IALPHAS[n_ips], IBETAS[n_ips], sim_seed);
					core_testing_method<stograd::stepper::adam<double> > (output, TESTS[n_tst], "adam", BSIZES[n_b], LRATES[n_lr],
						BETA_THETAS[n_tst], BETA_PHIS[n_tst], IALPHAS[n_ips], IBETAS[n_ips], sim_seed);
					core_testing_method<stograd::stepper::rmsprop<double> > (output, TESTS[n_tst], "rmsprop", BSIZES[n_b], LRATES[n_lr],
						BETA_THETAS[n_tst], BETA_PHIS[n_tst], IALPHAS[n_ips], IBETAS[n_ips], sim_seed);
					core_testing_method<stograd::stepper::amsgrad<double> > (output, TESTS[n_tst], "amsgrad", BSIZES[n_b], LRATES[n_lr],
						BETA_THETAS[n_tst], BETA_PHIS[n_tst], IALPHAS[n_ips], IBETAS[n_ips], sim_seed);
				}
			}
		}
	}
	return 0;
}
