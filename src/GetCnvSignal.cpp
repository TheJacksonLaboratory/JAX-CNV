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
	struct SReadDepth {
		SReadDepth(const int & i_pos, const unsigned int & i_count) : pos(i_pos), count(i_count){}
		int pos = 0;
		unsigned int count = 0;
	};
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

void PrintCleanBamData (SBamData & bam_data, const int & max_pos) {
	// If the max_pos means the last element.
	std::cout << (!bam_data.read_depth.empty() && max_pos == std::numeric_limits<std::int32_t>::max()
		? bam_data.read_depth.back().pos
		: max_pos)
		<< "\t";


	if (bam_data.total_read == 0) {
		std::cout << "0\t0\t0\t0\t0\t0\t0\t";
	} else {
		std::cout << bam_data.total_read << "\t" << bam_data.paired_reads << "\t"
				<< bam_data.proper_pairs / static_cast<double>(bam_data.total_read) << "\t" // proper pairs
				<< bam_data.inproper_pairs / static_cast<double>(bam_data.total_read) << "\t" // inproper pairs
				<< bam_data.mate_unmapped / static_cast<double>(bam_data.total_read) << "\t"; // mate unmapped
	
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

inline std::list<SBamData::SReadDepth>::iterator GetRdListIte (std::list<SBamData::SReadDepth> & read_depth, const int & pos) {
	for (std::list<SBamData::SReadDepth>::iterator ite = read_depth.begin();
		ite != read_depth.end(); ++ite) {
		if (pos == ite->pos)
			return ite;
	}

	// The pos is larger than read_depth.end()
	SBamData::SReadDepth tmp_data((read_depth.empty() ? pos : read_depth.back().pos + 1), 0); 
	while (tmp_data.pos <= pos) {
		read_depth.push_back(tmp_data);
		++tmp_data.pos;
	}

	// Return the last ite.
	std::list<SBamData::SReadDepth>::iterator ite = read_depth.begin();
	std::advance(ite, read_depth.size() - 1); //ite is set to last element
	return ite;
}

void ProcessAlignment (SBamData & bam_data, const bam1_t * aln) {
	// The alignment is not mapped.
	if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSECONDARY || aln->core.flag & BAM_FQCFAIL 
		|| aln->core.flag & BAM_FDUP || aln->core.flag & BAM_FSUPPLEMENTARY) 
		return;

/*
std::cerr << "aln pos: " << aln->core.pos << "\t";
for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
	std::cerr << bam_cigar_oplen(*(bam_get_cigar(aln) + i)) << "\t" << bam_cigar_op(*(bam_get_cigar(aln) + i)) << "\t";
}
std::cerr << std::endl;

//int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
//char *chr = header->target_name[aln->core.tid] ; //contig name (chromosome)
//uint32_t len = aln->core.l_qseq; //length of the read.
//uint8_t *q = bam_get_seq(aln); //quality string
//uint32_t q2 = aln->core.qual ; //mapping quality
*/
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
		std::list<SBamData::SReadDepth>::iterator ite = GetRdListIte(bam_data.read_depth, pos);
		for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
			const uint32_t op = bam_cigar_op(*(pCigar + i));
			if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL || op == BAM_CREF_SKIP) {
				for (uint32_t j = 0; j < bam_cigar_oplen(*(pCigar + i)); ++j) {
					if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) ++(ite->count);
					++ite;
					++pos;
					if (ite == bam_data.read_depth.end()) { // Need to add new element in the list.
						SBamData::SReadDepth tmp_data(pos, 0);
						bam_data.read_depth.push_back(tmp_data);
						ite =  bam_data.read_depth.begin();
						std::advance(ite, bam_data.read_depth.size() - 1); //iter is set to last element
					}
				} // end of for loop
			} // end of if
		} // end of for loop
	}

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
				for (int i = pre_bin; i < cur_bin; ++i)
					PrintCleanBamData(bam_data, (i + 1) * bin - 1); // (i + 1) * bin - 1 for giving the max pos of the bin.
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
			hts_itr_t * ite = sam_itr_querys(idx, header, cat_region.c_str());
			int pre_bin = region.begin / bin;
			while (ite && sam_itr_next(bam_reader, ite, aln) >= 0) {
				const int cur_bin = aln->core.pos / bin;
				// If the cur_bin is not the same as pre_bin, we clean up the pre_bin.
				if ((cur_bin > pre_bin) && (cur_bin != pre_bin)) {
					for (int i = pre_bin; i < cur_bin; ++i)
						PrintCleanBamData(bam_data, (i + 1) * bin - 1); // (i + 1) * bin - 1 for giving the max pos of the bin.
					pre_bin = cur_bin;
				}
				ProcessAlignment(bam_data, aln);
			}
	
			// Clean up
			hts_itr_destroy(ite);
		}
		
	}

	//PrintCleanBamData(bam_data, std::numeric_limits<std::int32_t>::max());

	// Clean up
	bam_destroy1(aln);
	bam_hdr_destroy(header);
	sam_close(bam_reader);
	bam_data.Clean();
}

void PrintResults(std::stringstream & bam_signal_out, std::stringstream & count_kmer_out) {
	while (!bam_signal_out.eof() || !count_kmer_out.eof()) {
		std::string tmp;
		std::getline(bam_signal_out, tmp);
		std::cout << tmp << std::endl;
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

	//coutbuf = std::cout.rdbuf(); //save old buf
	//std::stringstream bam_signal_out;
	//std::cout.rdbuf(bam_signal_out.rdbuf()); //redirect std::cout to bam_signal_out
	ProcessBam(cmdline.bam.c_str(), region, cmdline.bin);
	//std::cout.rdbuf(coutbuf); //reset to standard output again

	// Open output if given.
/*
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
*/
	return 0;
}
