#include <string>
#include <iostream>
#include <fstream>
#include <vector>

// Self include
#include "include/CommandLine.h"

// FASTAQ include
#include "fasta.h"
#include "reference.h"
#include "region.h"

// Jellyfish include
#include "jellyfish/file_header.hpp"
#include "jellyfish/jellyfish.hpp"
#include "jellyfish/mapped_file.hpp"

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
			const jellyfish::file_header & header, const binary_query & bq, const bool & running_length_encoding) {

	std::vector<std::string> ref_names;
	ref.GetReferenceNames(&ref_names);
	// If a region is given, there will be only one reference in the vector.
	for (unsigned int i = 0; i < ref_names.size(); i++) {
		const unsigned int target_len = std::min(ref.GetReferenceLength(ref_names[i].c_str()), region.end);
		const unsigned int target_begin = std::max(static_cast<unsigned int>(0), region.begin);

		// Cannot proceed when target_len < kmer_size
		if (target_len < kmer_size) break;

		//if (cmdline.region.empty()) // Ouput refernce name only if a region is not given.
		//	std::cout << ">" << ref_names[i] << std::endl;

		int score_count = 0;
		char score = '\0';
		for (unsigned int j = target_begin; j < target_len - kmer_size; ++j) {
			jellyfish::mer_dna m;
			m = ref.GetSubString(ref_names[i], j, kmer_size).c_str();
			if (header.canonical()) m.canonicalize();
			if (running_length_encoding) {
				if (CeilLog2(bq.check(m)) == score) {
					++score_count;
				} else {
					if (score != '\0')
						std::cout << score << "\t" << score_count << std::endl;
					score_count = 1;
					score = CeilLog2(bq.check(m));
				}
			} else { // not running_length_encoding
				std::cout << CeilLog2(bq.check(m));
			}
			//std::cout << ref.GetSubString(ref_names[i], j, kmer_size) << std::endl;
		}
		// Output the last score. Only running_length_encoding mode will use this.
		if (score_count > 0)
			std::cout << score << "\t" << score_count << std::endl;
		if (!running_length_encoding)
			std::cout << std::endl;
	}
}

int main (int argc, char** argv) {
	const SCmdLine cmdline(argc, argv);
	if (!cmdline.CheckArg()) {
		std::cerr << cmdline.Help(argv[0]);
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

	GetKmerCount(ref, region, kmer_size, header, bq, cmdline.running_length_encoding);

	// Clean up
	if (!cmdline.output.empty()) {
		std::cout.rdbuf(coutbuf); //reset to standard output again
		ofs.close();
	}

}
