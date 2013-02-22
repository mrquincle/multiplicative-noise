/*
 * SplitOperator.h
 *
 * Integrate stochastic differential equation according to article in Physical Review
 * Letters: "Integration of Langevin equations with multiplicative noise and the
 * viability of field theories for absorbing phase transitions" by Dornic, Chaté, and
 * Muñoz (2005)
 *
 * Steps (we used shorthand p for \rho, u for \mu):
 *
 *  1. Given: dt=0.25, dx=1, a=1.75623, b=1, D=0.25, sigma=sqrt(2), m=22
 *  2. Calculate p*(t) using current state p(t) as p0
 *    - beta        = a - 2dD / (dx)^2
 *    - lambda      = 2*beta / sigma^2 * (exp(beta*dt)-1)
 *    - poisson_arg = lambda * p(x,t) * exp(beta*dt)
 *    - alpha       = D/(dx)^2 * ( p(x-dx,t) + p(x+dx,t) )
 *    - u           = -1 + 2*alpha / sigma^2
 *    - u + 1       = 2*alpha / sigma^2
 *    - gamma_arg   = u + 1 + poisson[poisson_arg]
 *    - p*          = gamma[gamma_arg]
 *  3. Calculate p at next time step
 *    - p(t+dt)     = p*(t) / (1 + p*(t) b dt)
 *
 * This code is written after the kindly provided Fortran code by Ivan Dornic concerning
 * the original article.
 *
 * This software is published under the GNU Lesser General Public license (LGPL).
 *
 * It is not possible to add usage restrictions to an open-source license. Nevertheless,
 * we personally strongly object against this software used by the military, in the
 * bio-industry, for animal experimentation, or anything that violates the Universal
 * Declaration of Human Rights.
 *
 * Copyright © 2011 Anne van Rossum <anne@almende.com>
 *
 * @author 	Anne C. van Rossum
 * @author  Ivan Dornic
 * @date	Sep 16, 2011
 * @project	Replicator FP7
 * @company	Almende B.V.
 * @license LGPL v3 or newer
 * @case	Self-organised criticality
 */

#ifndef SPLITOPERATOR_H_
#define SPLITOPERATOR_H_

#include <helix.h>

#include <fstream>

#include <Time.h>

#include <boost/thread/barrier.hpp>

//typedef long double TYPE;
// In the original Fortran code only double precision seems to be used
typedef double TYPE;

class SplitOperator {
public:
	// Constructor
	SplitOperator();

	~SplitOperator();

	//! Run this function
	bool Run();
private:
	void RunRange(int start, int end);

	bool Sum(int iter);
public:
	// just play stupid and make everything public for now
	// we will use the syntax of the mentioned article, not of the fortran code kindly
	// provided by Ivan Dornic.

	// SDE = \rho^dot = D \nabla^2 \rho + a \rho - b \rho^2 + sigma sqrt(\rho) \mu

	// a*\rho in SDE
	TYPE a;

	TYPE b;

	TYPE sigma;

	TYPE D;

	TYPE mu;

	TYPE rho;

	TYPE dD;

	TYPE dx;

	TYPE dt;

	TYPE rho_0;

	//! Total time span
	int timespan;

	//! There are 2^m sites
	int m;

	// There are N sites (2^m)
	int N;

	//! The array with all state values
	helix<TYPE> * p;

	std::ofstream file;
private:

	TYPE poisson_arg_const;

	TYPE alpha_const;

	TYPE lambda;

	TYPE beta;

	Time time;

	bool stop_run;

	boost::barrier *rendezvous;
};

#endif /* SPLITOPERATOR_H_ */
