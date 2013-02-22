/**
 * SplitOperator.cpp
 *
 * For summary, see SplitOperator.h file
 *
 * Copyright © 2010 Anne van Rossum <anne@almende.com>
 *
 * @author 	Anne C. van Rossum
 * @author	Ivan Dornic
 * @date	Sep 16, 2011
 * @project	Replicator FP7
 * @company	Almende B.V.
 * @license	LGPL v3 or newer
 * @case	Self-organised criticality
 */

#include <SplitOperator.h>

#include <math.h>
#include <iomanip>
#include <sstream>
#include <helix.h>

#include <boost/random.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <boost/thread.hpp>
#include <boost/bind.hpp>

using namespace std;

static int seed1=437;
static boost::mt19937 rnd_gen1(seed1);   //Mersenne Twister generator

static int seed2=4357;
static boost::mt19937 rnd_gen2(seed2);

#define SIMPLE_TEST
#undef SIMPLE_TEST

#define MULTI_THREADED
#undef MULTI_THREADED

TYPE rgamma( TYPE mean, TYPE variance, boost::mt19937& rng = rnd_gen1) {
	// return 0 if mean is 0
	if (mean <= 0) return 0.0;

	// check for NaN
	if (!(mean > 0)) cout << "Do not try the gamma distribution with mean: " << mean << endl;

	const TYPE shape = ( mean*mean )/variance;
	TYPE scale = variance/mean;

	boost::gamma_distribution<TYPE> gd( shape );
	boost::variate_generator<boost::mt19937&, boost::gamma_distribution<TYPE> > var_gamma( rng, gd );

	TYPE result = scale*var_gamma();
	return result;
}

TYPE rgamma( TYPE mean, boost::mt19937& rng = rnd_gen1) {
	if (mean <= 0) return 0.0;
	if (!(mean > 0)) cout << "Do not try the gamma distribution with mean: " << mean << endl;
	boost::gamma_distribution<TYPE> gd( mean );
	boost::variate_generator<boost::mt19937&,boost::gamma_distribution<TYPE> > var_gamma( rng, gd );
	TYPE result = var_gamma();
	return result;
}

TYPE xpoisson( TYPE mean, boost::mt19937& rng = rnd_gen2) {
	if (mean <= 0) {
		cout << "Do not try poisson distribution with mean <= 0: " << mean << endl;
		assert (mean > 0);
	}
	boost::poisson_distribution<int,TYPE> pd( mean );
	boost::variate_generator<boost::mt19937&,boost::poisson_distribution<int,TYPE> > var_poisson( rng, pd );

	TYPE result = var_poisson();
	return result;
}

/**
 * Normal destructor. Update the parameters here at your wish.
 */
SplitOperator::SplitOperator() {
	// See: caption Fig.1 in dornic2005
	a = 1.75623;		// Fig.1, however see email conversation, 1.85...
	b = 1.0;            // or -cc
	D = 0.25;
	sigma = sqrt(2.0);
	dx = 1.0;
	dt = 0.2;
	dD = D;
	m = 17;

	// overwrite Fig.1 caption values
	// seems not to be the right (non-universal) critical value
	a = 1.84701;
	// a = bb
	//	a = 1.75623;
	// becomes slow when dt < 0.1
	dt = 0.1;

	N = pow(2,m);
	p = new helix<TYPE>(N);

//	timespan = 10e4;
	timespan = 10e3;

	// create new file per "a" value
	ostringstream filename; filename.clear();
	filename << "integration_" << a << ".log";
	file.open (filename.str().c_str());
}

/**
 * Default constructor just removes "helix" data structure that contains the
 * probabilities over time.
 */
SplitOperator::~SplitOperator() {
	file.close();
	delete p;
}

#ifndef SIMPLE_TEST
/**
 * Run the algorithm over a limited range
 */
void SplitOperator::RunRange(int start, int end) {
	int niter = timespan / dt;

	// start with "1" so that we do not have a "0" coordinate (irritating with log plotting)
	for (int iter = 1; iter < niter; ++iter) {
		TYPE p_star, poisson_arg, alpha, mu_plus_1, gamma_arg;
		for (int i = start; i < end; ++i) {
			// See: Eq.6
			poisson_arg 	= poisson_arg_const * p->get(i);
			// See: Eq.4 (second part is the sum over p_left + p_right)
			alpha       	= alpha_const * ( p->left(i) + p->right(i) );
			// See: 2nd line under Eq.3 (which contains -1, but in Eq.6 we will need mu+1)
			mu_plus_1		= 2.0*alpha / (sigma*sigma); //for default sigma: mu+1=alpha

			if (poisson_arg == 0) {
				// xpoisson(0) is set to 0, choice...
				gamma_arg 		= mu_plus_1; //exp(-poisson_arg);
			} else {
				// Eq.6 (argument)
				gamma_arg		= mu_plus_1 + xpoisson(poisson_arg);
			}
			// Eq.6 finally
			p_star=rgamma(gamma_arg)/lambda;

			// 5th line below Eq.4, a bit hidden...
			p->set(p_star/(1.0+b*dt*p_star), i);
		}

#ifdef MULTI_THREADED
		// barrier, only one thread will return "true" here
		if (rendezvous->wait()) {
			p->update();
			// break out...
			if (Sum(iter)) stop_run = true;
		}

		// wait all again, so that we know stop_run has been set
		rendezvous->wait();
#else
		p->update();
		if (Sum(iter)) stop_run = true;
#endif
		if (stop_run) break;
	}
}
#endif

/**
 * The only function to be called. This will do the entire integration.
 */
bool SplitOperator::Run() {

	// See: line above Eq.4 in dornic2005
	beta = a-2.0*dD/(dx*dx);

	if (beta < 1e-5) {
		// use a factor "dt": limit beta->0 for beta / (exp(beta*dt)-1) = 1/dt
		lambda = 2.0 / (sigma*sigma * dt);
	} else {
		// See: 2nd line below Eq.3 in dornic2005
		lambda = 2.0*beta / (sigma*sigma * (exp(beta*dt)-1.0));
	}
	//	TYPE lambda      		= 2.0*beta / (sigma*sigma * (1.0-exp(-beta*dt))); //why this in Fortran code?

	// See: Eq.6 (excluding p_0, only constant part)
	poisson_arg_const	= lambda * exp(beta*dt);

	// See: Eq.4 (first part before the summation)
	alpha_const		= D/(dx*dx);

	cout << "beta: " << beta << endl; //", " << "alpha const: " << alpha_const << endl;
	cout << "coeff: alpha: " << alpha_const << endl;
	cout << "lambda: " << lambda << ", " << "poisson arg const: " << poisson_arg_const << endl;

	// initial conditions: perfectly homogeneous (does not change anything)
	for (int i = 0; i < N; ++i) p->set(1.0, i);
	p->update();

	stop_run = false;

#ifndef MULTI_THREADED
	time.Start();
	RunRange(0, N);
	return true;
#endif
	// we use a group of 8 threads
	boost::thread_group group;
	int nof_threads = 8;

	// we use a barrier for synchronisation purposes
	rendezvous = new boost::barrier(nof_threads);

	// for now only allow N that can be perfectly divided amongst X threads
	int range_size = N / nof_threads;
	assert (range_size * nof_threads == N);

	time.Start();

	// start all threads with each their own section of the x-axis
	for (int i = 0; i < nof_threads; ++i) {
		int start = i * range_size;
		int end = start + range_size;
		group.create_thread(boost::bind(&SplitOperator::RunRange, this, start, end));
	}
	// blocks till all threads are finished
	group.join_all();

	delete rendezvous;
	return true;
}

/**
 * Sum over all the values and check if the result is 0. If that is the case we do
 * not need to continue anymore. The result has apparently converged to 0.
 */
bool SplitOperator::Sum(int iter) {
	int iunsdt=ceil(1.0/dt);
	TYPE unsn=1.0/(TYPE)N;

	TYPE sum1;
	// print and check for sum1==0, granularity decreases over time
	if (((iter <= 100*iunsdt) && ((iter % (5*iunsdt)) == 0)) ||
			((iter <= iunsdt*10000) && ((iter % (50*iunsdt)) == 0)) ||
			((iter > iunsdt*10000) && ((iter % (500*iunsdt)) == 0))) {

		time.Stop();
		time.Print();

		sum1=0.0;
		for (int i = 0; i < N; ++i) sum1+=p->get(i);

		if (sum1 < 1.0e-7) {
			cout << "End... (sum1 == " << sum1 << ")" << endl;
			return true;
		} else {
			// print to stdout
			cout << setprecision(0) << "[dt=" << iter*dt << "] ";
			cout << setiosflags(ios::fixed) << setw(12) << setprecision(11) << sum1*unsn << endl;
			//			cout << resetiosflags(ios::fixed);

			// print to file
			file << iter*dt << ", " << sum1*unsn << endl;
			time.Start();
		}
	}
	return false;
}


#ifdef SIMPLE_TEST
/**
 * Just sum of left and right value. To test the "helix" data structure
 * and its get and sets functions.
 */
void SplitOperator::RunRange(int start, int end) {
	int niter = timespan / dt;

	for (int iter = 1; iter < niter; ++iter) {
		for (int i = start; i < end; ++i) {
			p->set( p->left(i) + p->right(i), i);
		}

#ifdef MULTI_THREADED
		// have a boost barrier, we wait till all threads arrive here
		if (rendezvous->wait()) {
			p->update();

			cout << niter << ":";
			for (int i = 0; i < N; ++i) cout << p->get(i) << " ";
			cout << endl;

			if (iter == 3) {
				stop_run = true;
			}
		}

		// wait all again, so that we know stop_run has been set
		rendezvous->wait();
#else
		if (iter == 3) {
			stop_run = true;
		}
#endif
		if (stop_run) break;
	}
}
#endif

