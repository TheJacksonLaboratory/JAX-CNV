#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "fasta.h"
#include "reference.h"
#include "jellyfish/file_header.hpp"

int main (int argc, char** argv) {
	if (argc != 3) {
		std::cerr << "USAGE: " << argv[0] << " db.jf FASTA" << std::endl;
		return 1;
	}

	std::ifstream db(argv[1], std::ios::in|std::ios::binary);
	jellyfish::file_header header(db);
	if(!db.good()) {
		std::cerr << "Cannot open " << argv[1] << std::endl;
		return 1;
	}

	const int kmer_size = header.key_len() / 2;
	std::cerr << kmer_size << std::endl;

	CReference ref;
	Fasta::Load(ref, argv[2]);
	std::vector<std::string> ref_names;
	ref.GetReferenceNames(&ref_names);
	for (unsigned int i = 0; i < ref_names.size(); i++) {
		const int contig_len = ref.GetReferenceLength(ref_names[i].c_str());
		for (int j = 0; j < contig_len - kmer_size; ++j) {
			std::cout << ref.GetSubString(ref_names[i], j, kmer_size) << std::endl;
		}
	}
}
