#include <iostream>
#include <cstring>
#include <vector>
#include <string>

#include "umdhmm-v1.02/nrutil.h"
#include "umdhmm-v1.02/hmm.h"

#include "DataStruct.h"
#include "CallHmm.h"

namespace {
void PrintHmm (const HMM& hmm, const int& T, const int* O) {
	std::cerr << "DEBUG: CallHmm::HmmAndViterbi" << std::endl;
	std::cerr << "=====Start HMM table printing=====" << std::endl;
	std::cerr << "N, M: " << hmm.N << ", " << hmm.M << std::endl;
	std::cerr << "A" << std::endl;
	for (int i = 0; i < hmm.N; ++i) {
		std::cerr << "A" << i << std::endl;
		for (int j = 0; j < hmm.N; ++j)
			std::cerr << hmm.A[i][j] << " ";
		std::cerr << std::endl;
	}

	std::cerr << "B" << std::endl;
	for (int i = 0; i < hmm.N; ++i) {
		std::cerr << "B" << i << std::endl;
		for (int j = 0; j < hmm.M; ++j)
			std::cerr << hmm.B[i][j] << " ";
		std::cerr << std::endl;
	}

	std::cerr << "pi" << std::endl;
	for (int i = 0; i < hmm.N; ++i) {
		std::cerr << hmm.pi[i] << " ";
	}

	std::cerr << std::endl;
	std::cerr << "T: " << T << std::endl;
	std::cerr << "O" << std::endl;
	for (int i = 0; i < T; ++i) {
		std::cerr << O[i] << " ";
	}
	std::cerr << std::endl;
	std::cerr << "=====End HMM table printing=====" << std::endl;
}

void SmoothStats(std::vector<SHmmStats> & cnvs, const std::string & ref_name, 
			const std::vector <SReadDepth>& read_depth, const int bin_size, const int* q, const int T) {
	if (read_depth.size() != T) {
		std::cerr << "ERROR: HMM read_depth's size does not match with the number of stats." << std::endl;
		return;
	}

	std::vector <SHmmStats> result;
	std::vector <SReadDepth>::const_iterator rd_ite = read_depth.begin();
	// Ccollapse stats.
	for (int i = 1; i <= T; ++i, ++rd_ite) {
		// If there are >50% N's in the region, the region won't be taken in account so we set the stats to NORMAL.
#ifdef DEBUG
		std::cerr << rd_ite->pos << "\t" << rd_ite->n_count << "\t" << (((rd_ite->n_count * 2) > bin_size) ? 3 : q[i]) << std::endl;
#endif
		const int cur_stat = (rd_ite->n_count * 2) > bin_size ? 3 : q[i];
		
		//std::cerr << rd_ite->pos << "\t" << rd_ite->n_count << "\t" << cur_stat << "\t" << q[i] << std::endl;
		if (result.empty() || cur_stat != result.back().stats) { // Create the init hmm_stats.
			SHmmStats tmp(rd_ite->pos, cur_stat, 0);
			result.push_back(tmp);
		}
		result.back().length += bin_size;
	}

#ifdef DEBUG
	std::cerr << "HMM before smoothing" << std::endl;
	for (std::vector <SHmmStats>::const_iterator ite = result.begin(); ite != result.end(); ++ite) {
		std::cerr << ite->pos << "\t" << ite->stats << "\t" << ite->length << std::endl;
	}
#endif

	std::vector <SHmmStats> smooth_result;
	smooth_result.push_back(result.front());
	for (std::vector <SHmmStats>::const_iterator ite = std::next(result.begin()); ite != result.end(); ++ite) {
		if (ite->length < 5000 || ite->stats == smooth_result.back().stats)
			smooth_result.back().length += ite->length;
		else
			smooth_result.push_back(*ite);
	}

#ifdef DEBUG
	std::cerr << "HMM after smoothing" << std::endl;
	for (std::vector <SHmmStats>::const_iterator ite = smooth_result.begin(); ite != smooth_result.end(); ++ite) {
		std::cerr << ite->pos << "\t" << ite->stats << "\t" << ite->length << std::endl;
	}
#endif

	for (std::vector <SHmmStats>::const_iterator ite = smooth_result.begin(); ite != smooth_result.end(); ++ite) {
		if (ite->stats != 3 && ite->length > 250000) {
			cnvs.push_back(*ite);
			cnvs.back().chr = ref_name;
		}
	}
}
} // namespace

namespace CallHmm { 
bool HmmAndViterbi (std::vector<SHmmStats> & cnvs, const std::string & ref_name, const std::vector <SReadDepth>& read_depth, const int & bin_size, const int & coverage) {
	if (read_depth.empty()) return false;

	// Init HMM
	int T = read_depth.size();
	int* O = new int [T + 1]; // observation sequence O[1..T]
	for (std::vector <SReadDepth>::const_iterator ite = read_depth.begin(); ite != read_depth.end(); ++ite) {
		//const int rd_diff = (ite->count - coverage) / static_cast<double>(ite->kmer_score);
		//int tmp_o = round((coverage + rd_diff) / coverage * 50);
		//if (tmp_o < 0) tmp_o = 0;
		const int tmp_o = round(ite->count / static_cast<float>(coverage) * 50);
		O[std::distance(read_depth.begin(), ite) + 1] = std::min(tmp_o, 180);
	}

	HMM hmm;
	hmm.N  = hmm_N;
	hmm.M  = hmm_M;
	hmm.A  = new double* [hmm.N + 1];
	hmm.B  = new double* [hmm.N + 1];
	hmm.pi = new double [hmm.N + 1];

	for (int i = 1; i <= hmm.N; ++i) {
		hmm.A[i] = new double [hmm.N + 1];
		hmm.B[i] = new double [hmm.M + 1];
		std::memcpy(hmm.A[i] + 1, hmm_A[i - 1], sizeof(double)* hmm.N);
		std::memcpy(hmm.B[i] + 1, hmm_B[i - 1], sizeof(double) * hmm.M);
	}
	std::memcpy(hmm.pi + 1, hmm_pi, sizeof(double) * hmm.N);

	int* q = new int [T + 1]; // resultant states
	int** psi = new int* [T + 1];
	double **delta = new double* [T + 1];
	for (int i = 1; i <= T; ++i) {
		psi[i] = new int [hmm.N + 1];
		delta[i] = new double [hmm.N + 1];
	}
	double logproba = 0;
	// End of Init HMM
	
#ifdef DEBUG
PrintHmm(hmm, T, O);
#endif
	
	ViterbiLog(&hmm, T, O, delta, psi, q, &logproba);
	SmoothStats(cnvs, ref_name, read_depth, bin_size, q, T);

	// Clean up
	for (int i = 1; i <= hmm.N; ++i) {
		delete hmm.A[i];
		delete hmm.B[i];
	}
	
	for (int i = 1; i <= T; ++i) {
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
