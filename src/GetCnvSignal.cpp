#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <limits>
#include <cstdint>

// Self include
#include "CountKmer.h"
#include "GetCnvSignal.h"
#include "ReadDepth.h"
#include "CallHmm.h"


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
	unsigned int paired_reads = 0;
	unsigned int proper_pairs = 0;
	unsigned int inproper_pairs = 0;
	unsigned int mate_unmapped = 0;
	std::list<unsigned int> isizes;
	std::list<unsigned int> softclips;
	std::list<unsigned int> mismatches;
	std::list <SReadDepth> read_depth;

	void Clean() {
		total_read = 0;
		paired_reads = 0;
		proper_pairs = 0;
		inproper_pairs = 0;
		mate_unmapped = 0;
		isizes.clear();
		softclips.clear();
		mismatches.clear();
		read_depth.clear();
	}
};

// Func: Once alignments in a bin have been processed, then dump the info we collect for this bin.
//       Also, clean the info for the next bin.
void PrintCleanBamData (SBamData & bam_data, std::list <SReadDepth> & hmm_rd, const int & max_pos) {
	// If the max_pos means the last element.
	const int cur_pos = !bam_data.read_depth.empty() && max_pos == std::numeric_limits<std::int32_t>::max()
				? bam_data.read_depth.back().pos : max_pos;
	std::cout << cur_pos << "\t";

	if (bam_data.total_read == 0) {
		std::cout << "0\t0\t0\t0\t0\t0\t0\t";
	} else {
		std::cout << bam_data.total_read << "\t" << bam_data.paired_reads << "\t" // total reads and total paired-end reads
				<< bam_data.proper_pairs << "\t" // proper pairs
				<< bam_data.inproper_pairs << "\t" // inproper pairs
				<< bam_data.mate_unmapped << "\t"; // mate unmapped
	
		// Isize
		uint64_t sum = 0;
		for (std::list<unsigned int>::const_iterator ite = bam_data.isizes.begin(); ite != bam_data.isizes.end(); ++ite)
			sum += *ite;
		std::cout << sum / static_cast<double>(bam_data.total_read) << "\t";

		// Softclip
		sum = 0;
		for (std::list<unsigned int>::const_iterator ite = bam_data.softclips.begin(); ite != bam_data.softclips.end(); ++ite)
			sum += *ite;
		std::cout << sum / static_cast<double>(bam_data.total_read) << "\t";
	}

	// Read depth
	uint64_t sum = 0;
	unsigned int pos_count = 0;
	while (!bam_data.read_depth.empty() && bam_data.read_depth.front().pos <= max_pos) {
		++pos_count;
		sum += bam_data.read_depth.front().count;
		bam_data.read_depth.pop_front();
	}
	std::cout << (pos_count == 0 ? 0 : sum / static_cast<double>(pos_count)) << std::endl;

	SReadDepth rd_tmp(cur_pos, round(sum / static_cast<double>(pos_count)));
	hmm_rd.push_back(rd_tmp);

	// Clean
	//  bam_data.read_depth has been cleaned in the while loop.
	bam_data.total_read = 0;
	bam_data.paired_reads = 0;
	bam_data.proper_pairs = 0;
	bam_data.inproper_pairs = 0;
	bam_data.mate_unmapped = 0;
	bam_data.isizes.clear();
	bam_data.softclips.clear();
	bam_data.mismatches.clear();
};

// Func: Locate the SBamData in the list by using pos for read depth calculation.
//       If pos is not in SBamData::read_depth, then push new continues elements.
//       ProcessAlignment calls this function.
// TODO: The performance of this function needs to be improved.
inline std::list<SReadDepth>::iterator GetRdListIte (std::list<SReadDepth> & read_depth, const int & pos) {
	if (!read_depth.empty() && pos <= read_depth.back().pos) {
		for (std::list<SReadDepth>::iterator ite = read_depth.begin();
			ite != read_depth.end(); ++ite) {
			if (pos == ite->pos)
				return ite; // Get the match pos and return the iterator.
		}
	} else {
	
		// The pos is larger than read_depth.end()
		SReadDepth tmp_data((read_depth.empty() ? pos : read_depth.back().pos + 1), 0); 
		while (tmp_data.pos <= pos) {
			read_depth.push_back(tmp_data);
			++tmp_data.pos;
		}
	}

	// Return the last ite.
	std::list<SReadDepth>::iterator ite = read_depth.begin();
	std::advance(ite, read_depth.size() - 1); //ite is set to last element
	return ite;
}

// Func: Extract flag and read depth information and collect them in bam_data.
// @bam_data: The info of the alignment will be kept in this data.
// @aln: hts_lib bam alignment.
void ProcessAlignment (SBamData & bam_data, const bam1_t * aln) {
	// The alignment is not mapped.
	if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSECONDARY || aln->core.flag & BAM_FQCFAIL 
		|| aln->core.flag & BAM_FDUP || aln->core.flag & BAM_FSUPPLEMENTARY) 
		return;

	++bam_data.total_read;
	// Paired-end read
	if (aln->core.flag & BAM_FPAIRED) ++bam_data.paired_reads;
	// Proper pairs
	if (aln->core.flag & BAM_FPROPER_PAIR) ++bam_data.proper_pairs;
	else ++bam_data.inproper_pairs;
	// Is mate unmapped?
	// If it is, isize will be collected.
	if (aln->core.flag & BAM_FMUNMAP) ++bam_data.mate_unmapped;
	else bam_data.isizes.push_back(aln->core.isize < 0 ? -aln->core.isize : aln->core.isize);
	// Softclip
	const uint32_t* pCigar = bam_get_cigar(aln);
	uint32_t sc = 0;
	if (bam_cigar_op(*pCigar) == BAM_CSOFT_CLIP)
		sc += bam_cigar_oplen(*pCigar);
	// Check the last CIGAR
	if (bam_cigar_op(*(pCigar + aln->core.n_cigar - 1)) == BAM_CSOFT_CLIP)
		sc += bam_cigar_oplen(*(pCigar + aln->core.n_cigar - 1));
	bam_data.softclips.push_back(sc);

	// Read depth
	if (!bam_data.read_depth.empty() && aln->core.pos < bam_data.read_depth.front().pos) {
		std::cerr << "ERROR: The pos of " << bam_get_qname(aln) << " " << aln->core.pos << " is inconsistent." << std::endl
				<< "\t\tThe smallest pos of the current bin is " << bam_data.read_depth.front().pos << std::endl;
	} else {
		int32_t pos = aln->core.pos;
		// Locate the SBamData for RD in the list by using pos.
		std::list<SReadDepth>::iterator ite = GetRdListIte(bam_data.read_depth, pos);
		for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
			const uint32_t op = bam_cigar_op(*(pCigar + i));
			if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL || op == BAM_CREF_SKIP) {
				for (uint32_t j = 0; j < bam_cigar_oplen(*(pCigar + i)); ++j) {
					if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) ++(ite->count);
					++ite;
					++pos;
					// GetRdListIte cannot locate SBamData in the list for the given pos so append a new SBamData.
					if (ite == bam_data.read_depth.end()) { // Need to add new element in the list.
						SReadDepth tmp_data(pos, 0);
						bam_data.read_depth.push_back(tmp_data);
						ite =  bam_data.read_depth.begin();
						std::advance(ite, bam_data.read_depth.size() - 1); //iter is set to last element
					}
				} // end of for loop
			} // end of if
		} // end of for loop
	}

}

// Func: Process bam by giving bam filename, bam region and bin size.
//       For each alignment, Function ProcessAlignment will extract info the alignment.
//       HMM using read depths is also embedded.
// @ref: We can access bases of the entire chromosome from ref.
void ProcessBam (std::list <SReadDepth> & hmm_rd, const char * bam_filename, const Fastaq::SRegion & region, const int & bin, const std::string & ref) {
	samFile * bam_reader = sam_open(bam_filename, "r");

	bam_hdr_t *header;
	header = sam_hdr_read(bam_reader);
	bam1_t * aln = bam_init1();

	SBamData bam_data;
	// idx must be okay. We have checked in Run().
	hts_idx_t * idx = sam_index_load(bam_reader,  bam_filename);
	const bool load_index = idx == NULL ? false : true;


	if (load_index) {
		const std::string cat_region = region.chr + ":" + std::to_string(region.begin) + '-' +  std::to_string(region.end);
		hts_itr_t * ite = sam_itr_querys(idx, header, cat_region.c_str());
		int pre_bin = region.begin / bin;
		while (ite && sam_itr_next(bam_reader, ite, aln) >= 0) {
			const int cur_bin = aln->core.pos / bin;
			// If the cur_bin is not the same as pre_bin, we clean up the pre_bin.
			if ((cur_bin > pre_bin) && (cur_bin != pre_bin)) {
				// Every bin in the region between pre_bin and cur_bin will be padded.
				for (int i = pre_bin; i < cur_bin; ++i){
					PrintCleanBamData(bam_data, hmm_rd, (i + 1) * bin - 1); // (i + 1) * bin - 1 for giving the max pos of the bin.
					// Calculate the number of N's in this region.
					hmm_rd.back().n_count = 0;
					//TODO: The for loop seems slow.
					for (std::string::const_iterator s_ite = std::next(ref.begin(), i * bin); 
						s_ite != ref.end() && s_ite != std::next(ref.begin(), (i + 1) * bin - 1); ++s_ite) {
						if (*s_ite == 'N') 
							++(hmm_rd.back().n_count);
					}
				}
				pre_bin = cur_bin;
			}
			ProcessAlignment(bam_data, aln);
		}

		// Clean up
		hts_itr_destroy(ite);
	}
	
	//PrintCleanBamData(bam_data, std::numeric_limits<std::int32_t>::max());

	// Clean up
	bam_destroy1(aln);
	bam_hdr_destroy(header);
	sam_close(bam_reader);
	bam_data.Clean();
}

void PrintResults(std::ofstream & log, std::stringstream & bam_signal_out, std::stringstream & count_kmer_out, const bool have_count_kmer_out) {
	log << "#POS\tREADS\tPAIRED\tPROPER_PAIRS\tINPROPER_PAIRS\tMATE_UNMAPPED\tISIZE\tSOFTCLIPS\tREAD_DEPTH"
			<< (!have_count_kmer_out ? "\n" : "\tKMER_COUNT\n");
	while (!bam_signal_out.eof()) { // The bam_signal_out is the major player here.
		std::string tmp;
		if (std::getline(bam_signal_out, tmp).eof()) break;
		log << tmp;
		if (have_count_kmer_out && !std::getline(count_kmer_out, tmp).eof()) {
			log << "\t" << tmp;
		}
		log << std::endl;
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

	// coutbuf for storing the result of count_kmer.
	// Since count_kmer prints out result on std::cout, we need to redirect std::cout buffer.
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
	// Parse region from the command line or parse regions from the bam header.
	std::list<Fastaq::SRegion> regions;
	if (!cmdline.region.empty()) { // Parse region from the command line.
		Fastaq::SRegion tmp_region;
		if (!tmp_region.Parse(cmdline.region)) {
			std::cerr << "ERROR: The given region is not valid." << std::endl;
			return 1;
		}
		// Only chromosome name is given.
		// Need to find the length of the chromosome.
		if (tmp_region.begin == 0 && tmp_region.end == 0) {
			// Load bam header
			samFile * bam_reader = sam_open(cmdline.bam.c_str(), "r");
			bam_hdr_t *header;
			header = sam_hdr_read(bam_reader);
			for (int32_t i = 0; i < header->n_targets; ++i) {
				if (tmp_region.chr.compare(header->target_name[i]) == 0) {
					tmp_region.end = (header->target_len[i]) - 1;
					break;
				}
			}
		}
		regions.push_back(tmp_region);
	} else { // Parse regions from the bam header.
		// Load bam header
		samFile * bam_reader = sam_open(cmdline.bam.c_str(), "r");
		bam_hdr_t *header;
		header = sam_hdr_read(bam_reader);
		for (int32_t i = 0; i < header->n_targets; ++i) {
			Fastaq::SRegion tmp_region;
			tmp_region.chr = (header->target_name[i]);
			tmp_region.begin = 0;
			tmp_region.end = (header->target_len[i]) - 1;
			regions.push_back(tmp_region);
		}
		bam_hdr_destroy(header);
		sam_close(bam_reader);
	}
	// Divide regions into 5M block if it is larger than 5M.
	// HMM seems to get much faster performance for smaller regions.
	for (std::list<Fastaq::SRegion>::iterator ite = regions.begin(); ite != regions.end(); ++ite) {
		if (ite->end - ite->begin + 1 > 5000000) {
			Fastaq::SRegion tmp;
			tmp = *ite;
			ite->end = ite->begin + 5000000 - 1;
			tmp.begin = ite->end + 1;
			regions.insert(std::next(ite), tmp);
		}
	}

	// Check BAI
	samFile * bam_reader = sam_open(cmdline.bam.c_str(), "r");
	hts_idx_t * idx = sam_index_load(bam_reader,  cmdline.bam.c_str());
	if (idx == NULL && sam_index_build(cmdline.bam.c_str(), 0) < 0) { // Try to build bam index
		std::cerr << "ERROR: The region givin but bam index cannot be built and loaded." << std::endl;
		sam_close(bam_reader);
		return 1;
	}
	sam_close(bam_reader);

	// Re-direct cout
	coutbuf = std::cout.rdbuf(); //save old buf
	std::stringstream bam_signal_out;
	std::cout.rdbuf(bam_signal_out.rdbuf()); //redirect std::cout to bam_signal_out

	// Process BAM by regions
	std::string ref_seq;
	std::string ref_name;
	for (std::list<Fastaq::SRegion>::const_iterator ite = regions.begin(); ite != regions.end(); ++ite) {
		// The chromosome is not in ref. Load it from fasta.
		if (ref_name != ite->chr) {
			std::cerr << "chr: " << ite->chr << std::endl;
			ref_name = ite->chr; // Keep the new chr name.
			// Load a complete seq of the chromosome.
			if (!Fastaq::FastaLoad(ref_seq, cmdline.fasta.c_str(), true, ite->chr.c_str())) {
				std::cerr << "ERROR: Cannot load chromosome " << ite->chr << " from " << cmdline.fasta << std::endl;
				return 1;
			}
		}
		std::list <SReadDepth> hmm_rd; // The list to collect read depth info for HMM.
		ProcessBam(hmm_rd,cmdline.bam.c_str(), *ite, cmdline.bin, ref_seq);
		// Perform HMM	
		CallHmm::HmmAndViterbi(hmm_rd, cmdline.bin, cmdline.coverage);
		
	}
	std::cout.rdbuf(coutbuf); //reset to standard output again

	// Open a file for outputing log
	if (!cmdline.log.empty()) {
		std::ofstream log;
		log.open(cmdline.log, std::ofstream::out);
		PrintResults(log, bam_signal_out, count_kmer_out, (!cmdline.input_jfdb.empty() && !cmdline.fasta.empty()));
		log.close();
	}

	return 0;
}
