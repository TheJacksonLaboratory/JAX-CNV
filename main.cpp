#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "fasta.h"
#include "reference.h"
#include "jellyfish/file_header.hpp"
#include "jellyfish/jellyfish.hpp"
#include "jellyfish/mapped_file.hpp"

int main (int argc, char** argv) {
	if (argc != 3) {
		std::cerr << "USAGE: " << argv[0] << " db.jf FASTA" << std::endl;
		return 1;
	}

	// Read jellyfish database
	std::ifstream db(argv[1], std::ios::in|std::ios::binary);
	jellyfish::file_header header(db);
	if(!db.good()) { // The jellyfish database is broken.
		std::cerr << "ERROR: Cannot open " << argv[1] << std::endl;
		return 1;
	}
	if (header.format() != binary_dumper::format) {
		std::cerr << "ERROR: Cannot process jellyfish database built by bloom filter." << std::endl;
		return 1;
	}
	
	// Read kmer size from the header of jellyfish database
	const int kmer_size = header.key_len() / 2; // The kmer size is key_len() / 2.
	jellyfish::mer_dna::k(kmer_size);
	if (kmer_size < 1) { // Cannot proceed if kmer size is not larger than zero.
		std::cerr << "ERROR: The kmer size (" << kmer_size << ") should be larger than 0." << std::endl;
		return 1;
	}

	// Load jellyfish database as query db.
	jellyfish::mapped_file binary_map(argv[1]);
	binary_map.load(); // Load in memory for speedup.
	binary_query bq(binary_map.base() + header.offset(), header.key_len(), header.counter_len(), header.matrix(),
				header.size() - 1, binary_map.length() - header.offset());

	CReference ref; // fastaq lib.
	Fasta::Load(ref, argv[2]); // fastaq lib.
	std::vector<std::string> ref_names;
	ref.GetReferenceNames(&ref_names);
	for (unsigned int i = 0; i < ref_names.size(); i++) {
		const int contig_len = ref.GetReferenceLength(ref_names[i].c_str());
		for (int j = 0; j < contig_len - kmer_size; ++j) {
			jellyfish::mer_dna m;
			m = ref.GetSubString(ref_names[i], j, kmer_size).c_str();
			if (header.canonical()) m.canonicalize();
			std::cout << bq.check(m) << std::endl;
			//std::cout << ref.GetSubString(ref_names[i], j, kmer_size) << std::endl;
		}
	}
}
