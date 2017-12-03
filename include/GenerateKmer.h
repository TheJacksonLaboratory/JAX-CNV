#ifndef _GENERATEKMER_H_
#define _GENERATEKMER_H_

#include <getopt.h>
#include <string>


struct SGenerateKmerCml {
	SGenerateKmerCml(){}
	SGenerateKmerCml(const int argc, char** const argv){Parse(argc, argv);}

	bool help = false;

	// i/o
	std::string fasta;      // -f --fasta
	std::string output;     // -o --output

	// operation parameters
	int kmer_size = 15;	// -s --size
	
	// command line
	std::string cmd;

	const char* short_option = "hf:o:s:";

	// Help list
	const std::string Help (const char* program) const { return
		std::string("\n") +
		std::string("USAGE: ") + program + std::string(" -i <jellyfish_db> -f <FASTA>\n\n") +
		std::string("	-h --help			Print this help list.\n") +
		std::string("\n") +
		std::string("Input & Output:\n") +
		std::string("	-f --fasta <FASTA>	FASTA for kmer lookup.\n") +
		std::string("	-o --output <FILE>	Output file.\n")+
		std::string("operation:\n") +
		std::string("   -s --size <FILE>	Kmer size [15].\n");
	}

	// Check the required arguments.
	bool CheckArg () const {
		bool ok = true;

		return ok && !fasta.empty();
	}

	bool Parse (const int argc, char** const argv) {
		// Record the input command line.
		for (int i = 0; i < argc; i++) cmd += std::string(argv[i]) + " ";
		const struct option long_option[] = {
			{"help", required_argument, NULL, 'h'},
			// i/o
			{"fasta", required_argument, NULL, 'f'},
			{"output", required_argument, NULL, 'o'},

			// operation parameters
			{"size", required_argument, NULL, 's'},
			{0,0,0,0}
		};
		int option_index = 0;
		int c = -1;
		while ((c = getopt_long(argc, argv, short_option, long_option, &option_index)) != -1) {
			switch (c) {
				case 'h': help = true; break;
				case 'f': fasta = optarg; break;
				case 'o': output = optarg; break;
				case 's': kmer_size = atoi(optarg); break;
				default: std::cerr << "WARNING: Unkonw parameter: " << long_option[option_index].name << std::endl; break;
			}
		}

		if (help) {
			Help(argv[0]);
			return false;
		}

		return CheckArg();
	}
}; // SCountKmerCml

class GenerateKmer {
 public:
	// Constructors
	GenerateKmer();
	GenerateKmer(int argc, char** argv);
	GenerateKmer(const char * pInput_fasta, const char * pOutput = NULL, const int size = 15);

	// The function will report kmer count according to the parameter setting.
	// Return: 0 is successful.
	int Run() const;

	// If files are not assinged when declaring the class, you may use the function to assign them.
	void SetParameters(const SGenerateKmerCml & cml);
	void SetParameters(const char * pInput_fasta, const char * pOutput = NULL, const int size = 15);
 private:
	SGenerateKmerCml cmdline;
	// Not allow to use copy and assign constructors.
	GenerateKmer(const GenerateKmer&);
	GenerateKmer& operator= (const GenerateKmer&);
};
#endif
