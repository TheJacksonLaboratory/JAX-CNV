#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
// htslib include
#include "htslib/sam.h"
#include "fastaq/fasta.h"
#include "EstimateCoverage.h"

namespace {
bool LoadKmer(std::string & chr_name_prefix, std::string & kmer_seq, const char * kmer_table, const std::string & chr_name) {
	bool load_chr = false;
	const bool quiet = true;
	if (chr_name_prefix.empty()) {
		if (Fastaq::FastaLoad(kmer_seq, kmer_table, false, chr_name.c_str(), quiet)) {
			load_chr = true;
		} else if (Fastaq::FastaLoad(kmer_seq, kmer_table, false, ("chr" + chr_name).c_str()), quiet) {
			load_chr = true;
			chr_name_prefix = "chr";
		} else {
			std::cerr << "WARNING: Cannot load kmer seqeunces " << chr_name << " from " << kmer_table << std::endl;
			load_chr = false;
		}
	} else {
		if (Fastaq::FastaLoad(kmer_seq, kmer_table, false, (chr_name_prefix + chr_name).c_str()), quiet) {
			load_chr = true;
		} else {
			std::cerr << "WARNING: Cannot load kmer seqeunces " << chr_name << " from " << kmer_table << std::endl;
			load_chr = false;
		}
	}

	return load_chr;
}
void CalculateChrCoverage(std::vector<float> & coverages, std::string & chr_name_prefix, 
				samFile * bam_reader, bam_hdr_t *header, hts_idx_t * idx, 
				const char * kmer_table, const char * char_chr_name, const int & min_region_size) {
	const std::string chr_name = char_chr_name;
	std::string kmer_seq;
	const bool load_kmer = LoadKmer(chr_name_prefix, kmer_seq, kmer_table, chr_name);

	size_t pos = 0;

	coverages.clear();

	while (pos != std::string::npos) {
		pos = kmer_seq.find('"', pos);
		if (pos != std::string::npos) {
			int length = 0;
			while(kmer_seq[pos] == '"') {
				++pos;
				++length;
			}

			if (length > min_region_size) {
				const std::string cat_region = chr_name_prefix + chr_name + ":" + std::to_string(pos - length) + '-' +  std::to_string(pos);
				hts_itr_t * ite = sam_itr_querys(idx, header, cat_region.c_str());
				bam1_t * aln = bam_init1();
				unsigned int base_count = 0;
				while (ite && sam_itr_next(bam_reader, ite, aln) >= 0) {
					if (!(aln->core.flag & BAM_FUNMAP)) {
						const uint32_t* pCigar = bam_get_cigar(aln);
						for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
							const uint32_t op = bam_cigar_op(*(pCigar + i));
							if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
								base_count += bam_cigar_oplen(*(pCigar + i));
						}
					}
				}
				//std::cerr << chr_name_prefix + cat_region << "\t" << base_count / static_cast<float>(length) << std::endl;
				hts_itr_destroy(ite);
				bam_destroy1(aln);
				coverages.push_back(base_count / static_cast<float>(length));
			}
			
		}
	}
}

bool CoverageSort(const std::pair<int, float> & t1, const std::pair<int, float> & t2) {
	return t1.second < t2.second;
}

void CalculateCoverage(int & all_chr_cov, std::vector<unsigned int> & aneuploidies, 
	const std::vector<float> & coverages, unsigned int begin_chr_id, unsigned int end_chr_id) {

	if (end_chr_id > coverages.size() - 1) end_chr_id = coverages.size() - 1;
	if (begin_chr_id > end_chr_id) return; // Invalid ids.
	if (end_chr_id - begin_chr_id < 3) return; // Not enough element for IQR calculation.

	std::vector<std::pair<int, float> > covs; // id and coverages
	for (unsigned int i = begin_chr_id; i <= end_chr_id; ++i) {
		covs.push_back(std::pair<int, float>(begin_chr_id + i, coverages[i]));
	}
	std::sort(covs.begin(), covs.end(), CoverageSort);

	const unsigned int Q2_id = (end_chr_id - begin_chr_id + 1) / 2;
	const unsigned int Q1_id = (end_chr_id - begin_chr_id + 1) / 4;
	const unsigned int Q3_id = Q2_id + (Q1_id == 0 ? 1 : Q1_id);
	const float IQR = covs[Q3_id].second - covs[Q1_id].second;

	// Based on IQR, we caluculate coverage.
	float total_coverage = 0.0;
	int total_coverage_count = 0;
	for (std::vector<std::pair<int, float> >::const_iterator ite = covs.begin(); ite != covs.end(); ++ite) {
		//std::cout << "chr:" << ite->first << "\t" << ite->second << std::endl;
		if (ite->second > covs[Q1_id].second && ite->second < covs[Q3_id].second) {
			total_coverage += ite->second;
			++total_coverage_count;
		}
	}
	all_chr_cov = std::round(total_coverage / static_cast<float>(total_coverage_count));

	// Based on the coverage, we detect aneuploidies
	aneuploidies.clear();
	for (std::vector<std::pair<int, float> >::const_iterator ite = covs.begin(); ite != covs.end(); ++ite) {
		if (ite->second < (all_chr_cov * 0.6) || ite->second > (all_chr_cov * 1.4)) 
			aneuploidies.push_back(ite->first);
	}
	std::sort(aneuploidies.begin(), aneuploidies.end());
}
}

namespace EstimateCoverage {
int EstimateCoverage(std::vector<float> & coverages, const char * bam_filename, const char * kmer_table) {
	std::string chr_name_prefix;

	samFile * bam_reader = sam_open(bam_filename, "r");
	bam_hdr_t *header;
	header = sam_hdr_read(bam_reader);
	hts_idx_t * idx = sam_index_load(bam_reader,  bam_filename);

	coverages.resize(Human::HumanAutosomeSize + Human::HumanAllosomeSize, 0);

	for (int i  = 0; i < Human::HumanAutosomeSize; ++i) {
		std::vector<float> chr_cov;
		int min_region_size = 20000;
		while (chr_cov.size() < 10 && min_region_size > 2000) {
			min_region_size = min_region_size >> 1;
			chr_cov.clear();
			CalculateChrCoverage(chr_cov, chr_name_prefix, bam_reader, header, idx, kmer_table, Human::HumanAutosome[i], min_region_size);
		}
		float cov_chr_total = 0.0;
		for (std::vector<float>::const_iterator cov_ite = chr_cov.begin(); cov_ite != chr_cov.end(); ++cov_ite) 
			cov_chr_total += *cov_ite;
		coverages[i] = cov_chr_total / static_cast<float>(chr_cov.size());
#ifdef DEBUG
		std::cerr << Human::HumanAutosome[i] << "\t" << coverages[i] << std::endl;
#endif
	}

	for (int i  = 0; i < Human::HumanAllosomeSize; ++i) {
		std::vector<float> chr_cov;
		int min_region_size = 20000;
		while (chr_cov.size() < 10 && min_region_size > 2000) {
			min_region_size = min_region_size >> 1;
			chr_cov.clear();
			CalculateChrCoverage(chr_cov, chr_name_prefix, bam_reader, header, idx, kmer_table, Human::HumanAllosome[i], min_region_size);
		}
		float cov_chr_total = 0.0;
		for (std::vector<float>::const_iterator cov_ite = chr_cov.begin(); cov_ite != chr_cov.end(); ++cov_ite) 
			cov_chr_total += *cov_ite;
		coverages[i + Human::HumanAutosomeSize] = cov_chr_total / static_cast<float>(chr_cov.size());
#ifdef DEBUG
		std::cerr << Human::HumanAllosome[i] << "\t" << coverages[i + Human::HumanAutosomeSize] << std::endl;
#endif
	}

	int all_chr_cov = 0;
	std::vector<unsigned int> aneuploidies;
	// Calculate IQR and then coverage as well as detect aneuploidies.
	CalculateCoverage(all_chr_cov, aneuploidies, coverages, 0, Human::HumanAutosomeSize - 1);
	if (aneuploidies.size() > 0) std::cerr << "Aneuploidies:";
	for (std::vector<unsigned int>::const_iterator ite = aneuploidies.begin(); ite != aneuploidies.end(); ++ite)
		std::cerr << "\t" << Human::HumanAutosome[*ite]; 
	if (aneuploidies.size() > 0) std::cerr << std::endl;

	// Clean up
	bam_hdr_destroy(header);
	sam_close(bam_reader);

	//return std::round(all_chr_cov / static_cast<float>(Human::HumanAutosomeSize));
	return all_chr_cov;
}
} //namespace EstimateCoverage
