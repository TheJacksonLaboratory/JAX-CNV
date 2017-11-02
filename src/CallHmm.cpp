#include <iostream>
#include <cstring>
#include <list>

#include "umdhmm-v1.02/nrutil.h"
#include "umdhmm-v1.02/hmm.h"

#include "CallHmm.h"

namespace CallHmm { 
bool HmmAndViterbi (const std::list <SReadDepth>& read_depth) {
	if (read_depth.empty()) return false;

	// Init HMM
	int T = read_depth.size();
	int* O = new int [T]; // observation sequence O[1..T]
	for (std::list <SReadDepth>::const_iterator ite = read_depth.begin(); ite != read_depth.end(); ++ite) {
		O[std::distance(read_depth.begin(), ite)] = ite->count;
	}

	HMM hmm;
	hmm.N  = hmm_N;
	hmm.M  = hmm_M;
	hmm.A  = new double* [hmm.N];
	hmm.B  = new double* [hmm.N];
	hmm.pi = new double [hmm.N];

	for (int i = 0; i < hmm.N; ++i) {
		hmm.A[i] = new double [hmm.N];
		hmm.B[i] = new double [hmm.M];
		std::memcpy(hmm.A[i], hmm_A[i], sizeof hmm.N);
		std::memcpy(hmm.B[i], hmm_B[i], sizeof hmm.M);
	}
	std::memcpy(hmm.pi, hmm_pi, sizeof hmm.N);

	int* q = new int [T]; // resultant states
	int** psi = new int* [T];
	double **delta = new double* [T];
	for (int i = 0; i < T; ++i) {
		psi[i] = new int [hmm.N];
		delta[i] = new double [hmm.N];
	}
	double logproba = 0;
	// End of Init HMM
	
	
	ViterbiLog(&hmm, T, O, delta, psi, q, &logproba);

	// Clean up
	for (int i = 0; i < hmm.N; ++i) {
		delete hmm.A[i];
		delete hmm.B[i];
	}
	
	for (int i = 0; i < T; ++i) {
		delete psi[i];
		delete delta[i];
	}
	
	delete O;
	delete hmm.A;
	delete hmm.B;
	delete hmm.pi;
	delete q;
	delete psi;
	delete delta;
	return true;
}
} // namespace CallHmm
