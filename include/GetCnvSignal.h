#ifndef _GETCNVSIGNAL_H_
#define _GETCNVSIGNAL_H_

#include <getopt.h>
#include <string>

struct SGetCnvSignalCml {
	SGetCnvSignalCml(){}
	SGetCnvSignalCml(const int argc, char** const argv){Parse(argc, argv);}

	bool help = false;

	// i/o
	std::string input_jfdb; // -i --input
	std::string fasta;      // -f --fasta
	std::string bam;	// -b --bam
	std::string output;     // -o --output

	// operation parameters
	double coverage = 0.0;	// -c --coverage
	std::string region;     // -r --region
	bool ascii = false;	// --ascii
	int bin = 1;		// --bin
	std::string log;	// --log
	
	// command line
	std::string cmd;

	const char* short_option = "hb:i:f:o:c:r:";

	// Help list
	const std::string Help (const char* program) const { return
		std::string("USAGE: ") + program + std::string(" -i <jellyfish_db> -f <FASTA>\n\n") +
		std::string("	-h --help			Print this help list.\n") +
		std::string("\n") +
		std::string("Input & Output:\n") +
		std::string("	-b --bam <BAM>			Input BAM; required.\n") +
		std::string("	-i --input <jellyfish_db>	Jellyfish created count database.\n") +
		std::string("	-f --fasta <FASTA>		FASTA for kmer lookup.\n") +
		std::string("	-o --output <FILE>		Output file.\n") +
		std::string("\n") +
		std::string("Operations:\n") +
		std::string("   -c --coverage <FLOAT>		Specify the coverage.\n") +
		std::string("	-r --region chr:begin-end	Specify a target region.\n") +
		std::string("   --ascii				Report kmer count in ASCII: (log2(#) + 1) + 33.\n") +
		std::string("	--bin <INT>			Report a result for each # bp. [1]\n") +
		std::string("   --log <FILE>			Log output.\n");
	}

	// Check the required arguments.
	bool CheckArg () const {
		bool ok = true;
		if (bin < 1) {
			std::cerr << "ERROR: --bin <INT> should not smaller than 1." << std::endl;
			ok = false;
		}
		if (bam.empty()) {
			std::cerr << "ERROR: -b <BAM> is required." << std::endl;
			ok = false;
		}
		if (fasta.empty()) {
			std::cerr << "ERROR: -f <FASTA> is required." << std::endl;
			ok = false;
		}
		if (coverage == 0.0) {
			std::cerr << "ERROR: -c <FLOAT> is required. Please specify the coverage." << std::endl;
			ok = false;
		}

		return ok;
	}

	bool Parse (const int argc, char** const argv) {
		// Record the input command line.
		for (int i = 0; i < argc; i++) cmd += std::string(argv[i]) + " ";
		const struct option long_option[] = {
			{"help", required_argument, NULL, 'h'},
			// i/o
			{"bam", required_argument, NULL, 'b'},
			{"input", required_argument, NULL, 'i'},
			{"fasta", required_argument, NULL, 'f'},
			{"output", required_argument, NULL, 'o'},

			// operation parameters
			{"coverage", required_argument, NULL, 'c'},
			{"region", required_argument, NULL, 'r'},
			{"ascii", no_argument, NULL, 1},
			{"bin", required_argument, NULL, 2},
			{"log", required_argument, NULL, 3},
			{0,0,0,0}
		};
		int option_index = 0;
		int c = -1;
		while ((c = getopt_long(argc, argv, short_option, long_option, &option_index)) != -1) {
			switch (c) {
				case 'h': help = true; break;
				case 'i': input_jfdb = optarg; break;
				case 'f': fasta = optarg; break;
				case 'b': bam = optarg; break;
				case 'o': output = optarg; break;
				case 'c': coverage = atof(optarg); break;
				case 'r': region = optarg; break;
				case 1: ascii = true; break;
				case 2: bin = atoi(optarg); break;
				case 3: log = optarg; break;
				default: std::cerr << "WARNING: Unkonw parameter: " << long_option[option_index].name << std::endl; break;
			}
		}

		if (help) {
			Help(argv[0]);
			return false;
		}

		return CheckArg();
	}
}; // SGetCnvSignalCml

class GetCnvSignal {
 public:
	// Constructors
	GetCnvSignal();
	GetCnvSignal(int argc, char** argv);

	// The function will report kmer count according to the parameter setting.
	// Return: 0 is successful.
	int Run() const;

	// If files are not assinged when declaring the class, you may use the function to assign them.
	//void SetParameters(const SGetCnvSignalCml & cml );

 private:
	SGetCnvSignalCml cmdline;
	GetCnvSignal(const GetCnvSignal&);
	GetCnvSignal& operator= (const GetCnvSignal&);
};
#endif // _GETCNVSIGNAL_H_
