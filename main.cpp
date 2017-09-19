#include <string>
#include <iostream>
#include <vector>

#include "lib/fastaq/include/fasta.h"
#include "lib/fastaq/include/reference.h"

int main (int argc, char** argv) {
	if (argc != 4) {
		std::cerr << "USAGE: " << argv[0] << " data.jf FASTA kmer_size" << std::endl;
		return 1;
	}

	const int kmer_size = atoi(argv[3]);

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
