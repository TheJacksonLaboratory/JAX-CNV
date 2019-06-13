#ifndef _ESTIMATECOVERAGE_H_
#define _ESTIMATECOVERAGE_H_

namespace EstimateCoverage {

int EstimateCoverage(std::vector<float> & coverages, bool & female, bool & male, 
			const char * bam_filename, const char * kmer_table); // kmer_table is in FASTA format. 

namespace Human {
static const int HumanAutosomeSize = 22;
static const char* HumanAutosome[HumanAutosomeSize] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"};
static const int HumanAllosomeSize = 2;
static const char* HumanAllosome[HumanAllosomeSize] = {"X","Y"};
} //namespace Human
} // namespace EstimateCoverage

#endif
