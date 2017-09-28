#include <string>
#include <iostream>
#include <fstream>
#include <vector>

// Self include
#include "CountKmer.h"

// FASTAQ include
#include "fastaq/fasta.h"
#include "fastaq/reference.h"
#include "fastaq/region.h"

// Jellyfish include
#include "jellyfish/file_header.hpp"
#include "jellyfish/jellyfish.hpp"
#include "jellyfish/mapped_file.hpp"

namespace {

// Calculate log2 score.
char CeilLog2 (const uint64_t in)
{
	uint64_t x = in;
	static const unsigned long long t[6] = {
		0xFFFFFFFF00000000ull,
		0x00000000FFFF0000ull,
		0x000000000000FF00ull,
		0x00000000000000F0ull,
		0x000000000000000Cull,
		0x0000000000000002ull
	};

	//int y = (((x & (x - 1)) == 0) ? 0 : 1);
	int y = (x > 0) ? 1 : 0;

	for (int i = 0, j = 32; i < 6; i++) {
		int k = (((x & t[i]) == 0) ? 0 : j);
		y += k;
		x >>= k;
		j >>= 1;
	}

	// Use 125 as the max.
	return (y > 92) ? 125 : y + 33;
}
void GetKmerCount (const Fastaq::CReference & ref, const Fastaq::SRegion & region, const unsigned int & kmer_size, 
			const jellyfish::file_header & header, const binary_query & bq, const bool running_length_encoding,
			const bool ascii, const int & bin) {


	std::vector<std::string> ref_names;
	ref.GetReferenceNames(&ref_names);

	// If a region is given, there will be only one reference in the vector.
	for (unsigned int i = 0; i < ref_names.size(); i++) {
		const unsigned int target_end = (!region.chr.empty() && region.end > 0) // Use region.end when we set region
						? std::min(ref.GetReferenceLength(ref_names[i].c_str()), region.end + kmer_size) 
						: ref.GetReferenceLength(ref_names[i].c_str());
		const unsigned int target_begin = !region.chr.empty() ? std::max(static_cast<unsigned int>(0), region.begin) : 0;

#ifdef DEBUG
	std::cerr << "GetKmerCount:" << std::endl;
	std::cerr << "\tProcessing region " << ref_names[i] << ":" << target_begin << "-" << target_end << std::endl;
#endif
		// Cannot proceed when target_len < kmer_size
		if ((target_end - target_begin) < kmer_size) break;

		unsigned int score_count = 0;
		unsigned int score_sum = 0;
		char score = '\0';
		for (unsigned int j = target_begin; j < target_end - kmer_size; ++j) {
			jellyfish::mer_dna m;
			m = ref.GetSubString(ref_names[i], j, kmer_size).c_str();
			if (header.canonical()) m.canonicalize();
			if (running_length_encoding) { // running_length_encoding must use ascii for reporting.
				if (CeilLog2(bq.check(m)) == score) {
					++score_count;
				} else {
					if (score != '\0')
						std::cout << score << "\t" << score_count;
					score_count = 1;
					score = CeilLog2(bq.check(m));
				}
			} else { // not running_length_encoding
				score_sum += bq.check(m);
				++score_count;
				if (((j-target_begin + 1) % bin) == 0) {
					if (ascii)
						std::cout << static_cast<char>(CeilLog2((std::lround(score_sum / static_cast<float>(score_count)))));
					else
						std::cout << score_sum / static_cast<float>(score_count) << std::endl;
					score_sum = 0;
					score_count = 0;
				}
			}
		}
		// Output the last score. Only running_length_encoding mode will use this.
		if (running_length_encoding) {
			if (score_count > 0)
				std::cout << score << "\t" << score_count;
			std::cout << std::endl;
		} else { // not running_length_encoding
			if (score_count > 0) {
				if (ascii)
					std::cout << static_cast<char>(CeilLog2((std::lround(score_sum / static_cast<float>(score_count)))));
				else
					std::cout << score_sum / static_cast<float>(score_count) << std::endl;
			}
			if (ascii)
				std::cout << std::endl;
		}
	}
}
} // namespace

int CountKmer::Run () const {
	if (!cmdline.CheckArg()) {
		std::cerr << cmdline.Help("CountKmer");
		return 1;
	}

	// Read jellyfish database
	std::ifstream db(cmdline.input_jfdb, std::ios::in|std::ios::binary);
	jellyfish::file_header header(db);
	if(!db.good()) { // The jellyfish database is broken.
		std::cerr << "ERROR: Cannot open " << cmdline.input_jfdb << std::endl;
		return 1;
	}
	if (header.format() != binary_dumper::format) {
		std::cerr << "ERROR: Cannot process jellyfish database built by bloom filter." << std::endl;
		return 1;
	}
	
	// Read kmer size from the header of jellyfish database
	const unsigned int kmer_size = header.key_len() / 2; // The kmer size is key_len() / 2.
	jellyfish::mer_dna::k(kmer_size);
	if (kmer_size < 1) { // Cannot proceed if kmer size is not larger than zero.
		std::cerr << "ERROR: The kmer size (" << kmer_size << ") should be larger than 0." << std::endl;
		return 1;
	}

	// Load jellyfish database as query db.
	jellyfish::mapped_file binary_map(cmdline.input_jfdb.c_str());
	binary_map.load(); // Load in memory for speedup.
	binary_query bq(binary_map.base() + header.offset(), header.key_len(), header.counter_len(), header.matrix(),
				header.size() - 1, binary_map.length() - header.offset());

	// Parse region.
	Fastaq::SRegion region;
	if (!cmdline.region.empty()) {
		if (!region.Parse(cmdline.region)) {
			std::cerr << "ERROR: The given region is not valid." << std::endl;
			return 1;
		}
	}

	// Load reference from FASTA.
	Fastaq::CReference ref; // fastaq lib.
	if (!cmdline.region.empty()) {
		if (!Fastaq::FastaLoad(ref, cmdline.fasta.c_str(), true, region.chr.c_str())) {
			std::cerr << "ERROR: Cannot load chromosome " << region.chr << " from FASTA." << std::endl;
			return 1;
		}
	} else {
		if (!Fastaq::FastaLoad(ref, cmdline.fasta.c_str())) {
			std::cerr << "ERROR: Cannot load FASTA." << std::endl;
			return 1;
		}
	}

	// Open output if given.
	std::ofstream ofs;
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	if (!cmdline.output.empty()) {
		ofs.open(cmdline.output, std::ofstream::out | std::ofstream::app);
		std::cout.rdbuf(ofs.rdbuf()); //redirect std::cout to file;
	}

	GetKmerCount(ref, region, kmer_size, header, bq, cmdline.rle, cmdline.ascii, cmdline.bin);

	// Clean up
	if (!cmdline.output.empty()) {
		std::cout.rdbuf(coutbuf); //reset to standard output again
		ofs.close();
	}

	db.close();

	return 0;

}

CountKmer::CountKmer(int argc, char** argv)
	: cmdline(argc, argv)
{
}

CountKmer::CountKmer(
	const char * pInput_jfdb, const char * pInput_fasta, const char * pOutput,
	const char * pRegion, const int input_bin, const bool input_ascii, const bool input_rle)
{
	SetParameters(pInput_jfdb, pInput_fasta, pOutput, pRegion, input_bin, input_ascii, input_rle);
}

void CountKmer::SetParameters(const SCountKmerCml & cml) {
	cmdline = cml;
}

void CountKmer::SetParameters(
	const char * pInput_jfdb, const char * pInput_fasta, const char * pOutput,
	const char * pRegion, const int input_bin, const bool input_ascii, const bool input_rle) {

	cmdline.input_jfdb = pInput_jfdb;
	cmdline.fasta = pInput_fasta;
	if (pOutput) cmdline.output = pOutput;
	if (pRegion) cmdline.region = pRegion;
	cmdline.bin = input_bin;
	cmdline.ascii = input_ascii;
	cmdline.rle = input_rle;
	
}
