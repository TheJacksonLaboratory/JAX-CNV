#include <string>
#include <iostream>
#include <fstream>
#include <vector>

// Self include
#include "CountKmer.h"
#include "GetCnvSignal.h"


// FASTAQ include
#include "fastaq/fasta.h"
#include "fastaq/reference.h"
#include "fastaq/region.h"

// Jellyfish include
#include "jellyfish/file_header.hpp"
#include "jellyfish/jellyfish.hpp"
#include "jellyfish/mapped_file.hpp"

// htslib include
#include "htslib/sam.h"

namespace {
struct SBamData {
	unsigned int total_read = 0;
	unsigned int proper_pairs = 0;
	unsigned int inproper_pairs = 0;
	unsigned int mate_unmapped = 0;
	std::vector<unsigned int> poss;
	std::vector<unsigned int> rds;
	std::vector<unsigned int> isizes;

	void Clean() {
		total_read = 0;
		proper_pairs = 0;
		inproper_pairs = 0;
		mate_unmapped = 0;
		poss.clear();
		rds.clear();
		isizes.clear();
	}
};

void PrintBamData (const SBamData & bam_data) {
	if (bam_data.total_read == 0) {
		std::cout << "0\t0\t0\t0" << std::endl;
	} else {
		std::cout << bam_data.proper_pairs / static_cast<double>(bam_data.total_read) << "\t"
				<< bam_data.inproper_pairs / static_cast<double>(bam_data.total_read) << "\t"
				<< bam_data.mate_unmapped / static_cast<double>(bam_data.total_read) << "\t";
	
		uint64_t sum = 0;
		for (unsigned int i = 0; i < bam_data.isizes.size(); ++i)
			sum += bam_data.isizes[i];
		std::cout << sum / static_cast<double>(bam_data.total_read) << std::endl;
	}
};

void ProcessAlignment (SBamData & bam_data, const bam1_t * aln) {
	if (aln->core.flag & BAM_FUNMAP) return;
	++bam_data.total_read;
	if (aln->core.flag & BAM_FPROPER_PAIR) ++bam_data.proper_pairs;
	else ++bam_data.inproper_pairs;
	if (aln->core.flag & BAM_FMUNMAP) ++bam_data.mate_unmapped;
	else bam_data.isizes.push_back(aln->core.isize < 0 ? -aln->core.isize : aln->core.isize);

	//int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
	//char *chr = header->target_name[aln->core.tid] ; //contig name (chromosome)
	//uint32_t len = aln->core.l_qseq; //length of the read.
	//uint8_t *q = bam_get_seq(aln); //quality string
	//uint32_t q2 = aln->core.qual ; //mapping quality
}

void ProcessBam (const char * bam_filename, const Fastaq::SRegion & region, const int & bin) {
	samFile * bam_reader = sam_open(bam_filename, "r");

	bam_hdr_t *header;
	header = sam_hdr_read(bam_reader);
	bam1_t * aln = bam_init1();

	SBamData bam_data;

	if (region.chr.empty()) { // the region is not set
		int pre_bin = 0;
		int pre_tid = 0;
		while (sam_read1(bam_reader, header, aln) >= 0) {
			// Enter a new chr.
			if (aln->core.tid != pre_tid) {
				pre_bin = 0;
				pre_tid = aln->core.tid;
			}
			const int cur_bin = aln->core.pos / bin;
			if (cur_bin != pre_bin) {
				for (int i = pre_bin; i < cur_bin - 1; ++i) {
					SBamData dummy_data;
					PrintBamData(dummy_data);
				}
				PrintBamData(bam_data);
				bam_data.Clean();
				pre_bin = cur_bin;
			}
			ProcessAlignment(bam_data, aln);
		}
	} else { // the region is given.
		hts_idx_t * idx = sam_index_load(bam_reader,  bam_filename);
		bool load_index = true;
		if (idx == NULL) {
			if (sam_index_build(bam_filename, 0) < 0) { // Try to build bam index
				std::cerr << "ERROR: The region givin but bam index cannot be built and loaded." << std::endl;
				load_index = false;
			} else {
				idx = sam_index_load(bam_reader,  bam_filename);
			}
			
		}

		if (load_index) {
			const std::string cat_region = region.chr + ":" + std::to_string(region.begin) + '-' +  std::to_string(region.end);
			hts_itr_t * iter = sam_itr_querys(idx, header, cat_region.c_str());
			int pre_bin = region.begin / bin;
			while (iter && sam_itr_next(bam_reader, iter, aln) >= 0) {
				const int cur_bin = aln->core.pos / bin;
				if ((cur_bin > pre_bin) && (cur_bin != pre_bin)) {
					for (int i = pre_bin; i < cur_bin - 1; ++i) {
						SBamData dummy_data;
						PrintBamData(dummy_data);
					}
					PrintBamData(bam_data);
					bam_data.Clean();
					pre_bin = cur_bin;
				}
				ProcessAlignment(bam_data, aln);
			}
	
			// Clean up
			hts_itr_destroy(iter);
		}
		
	}

	PrintBamData(bam_data);

	// Clean up
	bam_destroy1(aln);
	bam_hdr_destroy(header);
	sam_close(bam_reader);
}

void PrintResults(std::stringstream & bam_signal_out, std::stringstream & count_kmer_out) {
	while (!bam_signal_out.eof() || !count_kmer_out.eof()) {
		std::string tmp;
		std::getline(bam_signal_out, tmp);
		std::cout << tmp << "\t";
		std::getline(count_kmer_out, tmp);
		std::cout << tmp;
		std::cout << std::endl;
	}
}
} // namespace

GetCnvSignal::GetCnvSignal(int argc, char** argv)
	: cmdline(argc, argv)
{
}

int GetCnvSignal::Run () const {
	if (!cmdline.CheckArg()) {
		std::cerr << cmdline.Help("GetCnvSignal");
		return 1;
	}

	std::streambuf * coutbuf = std::cout.rdbuf(); //save old buf
	std::stringstream count_kmer_out;
	// Perform CountKmer
	if (!cmdline.input_jfdb.empty() && !cmdline.fasta.empty()) {
		std::cout.rdbuf(count_kmer_out.rdbuf()); //redirect std::cout to count_kmer_out
		// rle = false
		// output is controlled by count_kmer_out
		CountKmer count_kmer(cmdline.input_jfdb.c_str(), cmdline.fasta.c_str(), NULL, 
					cmdline.region.c_str(), cmdline.bin, cmdline.ascii, false);
		count_kmer.Run();
		std::cout.rdbuf(coutbuf); //reset to standard output again
	}

	// Parse region.
	Fastaq::SRegion region;
	if (!cmdline.region.empty()) {
		if (!region.Parse(cmdline.region)) {
			std::cerr << "ERROR: The given region is not valid." << std::endl;
			return 1;
		}
	}

	coutbuf = std::cout.rdbuf(); //save old buf
	std::stringstream bam_signal_out;
	std::cout.rdbuf(bam_signal_out.rdbuf()); //redirect std::cout to bam_signal_out
	ProcessBam(cmdline.bam.c_str(), region, cmdline.bin);
	std::cout.rdbuf(coutbuf); //reset to standard output again

	// Open output if given.
	if (!cmdline.output.empty()) {
		coutbuf = std::cout.rdbuf(); //save old buf
		std::ofstream ofs;
		ofs.open(cmdline.output, std::ofstream::out | std::ofstream::app);
		std::cout.rdbuf(ofs.rdbuf()); //redirect std::cout to file;
		PrintResults(bam_signal_out, count_kmer_out);
		std::cout.rdbuf(coutbuf); //reset to standard output again
		ofs.close();
	} else {
		PrintResults(bam_signal_out, count_kmer_out);
	}

	return 0;
}
