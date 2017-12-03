#include <vector>
#include <list>
#include <iostream>
#include <fstream>

// Self include
#include "GenerateKmer.h"

// FASTAQ include
#include "fastaq/fasta.h"
#include "fastaq/reference.h"

namespace {
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
}

int GenerateKmer::Run() const
{
	if (!cmdline.CheckArg()) {
		std::cerr << cmdline.Help("GenerateKmer");
		return 1;
	}

	std::vector<std::string> ref_names;
	std::list<std::list<uint64_t> > kmer;

	Fastaq::CountKmer(cmdline.fasta.c_str(), cmdline.kmer_size, ref_names, kmer);

	if (ref_names.size() != kmer.size()) {
		std::cerr << "ERROR: The # of ref_names does not match the # of kmer_table." << std::endl;
		return 1;
	}
	
	std::list<std::list<uint64_t> >::const_iterator table_ite = kmer.begin();


	// Re-direct cout
	std::ofstream ofs;
	std::streambuf * coutbuf = std::cout.rdbuf(); //save old buf
	if (!cmdline.output.empty()) {
		ofs.open(cmdline.output, std::ofstream::out);
		std::cout.rdbuf(ofs.rdbuf()); //redirect std::cout to file;
	}
	
	for (std::vector<std::string>::const_iterator name_ite = ref_names.begin(); name_ite != ref_names.end(); ++name_ite) {
		std::cout << *name_ite << std::endl;
		for (std::list<uint64_t>::const_iterator kmer_ite = table_ite->begin(); kmer_ite != table_ite->end(); ++kmer_ite) {
			std::cout << CeilLog2(*kmer_ite);
		}
		std::cout << std::endl;
		++table_ite;
	}

	// Clean up
	if (!cmdline.output.empty()) {
		std::cout.rdbuf(coutbuf); //reset to standard output again
		ofs.close();
	}

	return 0;
}

GenerateKmer::GenerateKmer(int argc, char** argv)
	: cmdline(argc, argv)
{
}

GenerateKmer::GenerateKmer(const char * pInput_fasta, const char * pOutput, const int size)
{
	SetParameters(pInput_fasta, pOutput, size);
}

void GenerateKmer::SetParameters(const SGenerateKmerCml & cml) {
	cmdline = cml;
}

void GenerateKmer::SetParameters(const char * pInput_fasta, const char * pOutput, const int size)
{
	cmdline.fasta = pInput_fasta;
	if (pOutput) cmdline.output = pOutput;
	cmdline.kmer_size = size;
}
