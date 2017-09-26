#ifndef _COMMANDLINE_H_
#define _COMMANDLINE_H_

#include <getopt.h>
#include <string>

struct SCmdLine {
	SCmdLine(const int argc, char** const argv){Parse(argc, argv);}

	bool help = false;

	// i/o
	std::string input_jfdb; // -i --input
	std::string fasta;      // -f --fasta
	std::string bam;	// -b --bam
	std::string output;     // -o --output

	// operation parameters
	std::string region;     // -r --region
	int bin = 1;
	bool running_length_encoding = false; // --rle
	
	// command line
	std::string cmd;

	const char* short_option = "hi:f:b:o:r:";

	// Help list
	const std::string Help (const char* program) const { return
		std::string("USAGE: ") + program + std::string(" -i <jellyfish_db> -f <FASTA>\n\n") +
		std::string("	-h --help			Print this help list.\n") +
		std::string("\n") +
		std::string("Input & Output:\n") +
		std::string("	-i --input <jellyfish_db>	Jellyfish created count database.\n") +
		std::string("	-f --fasta <FASTA>		FASTA for kmer lookup.\n") +
		std::string("	-b --bam <BAM>			Input BAM.\n") +
		std::string("	-o --output <FILE>		Output file.\n") +
		std::string("\n") +
		std::string("Operations:\n") +
		std::string("	-r --region chr:begin-end	Specify a target region.\n") +
		std::string("	--bin <INT>			Report a result for each # bp. [1]\n") +
		std::string("	--rle				Ouput by running length encoding.\n");
	}

	// Check the required arguments.
	bool CheckArg () const {
		bool ok = true;
		if (bin < 1) {
			std::cerr << "ERROR: --bin <INT> should not smaller than 1." << std::endl;
			ok = false;
		}
		if (bin > 1 && running_length_encoding) {
			std::cerr << "ERROR: --rle only work for --bin 1." << std::endl;
			ok = false;
		}

		return ok && ((!input_jfdb.empty() && !fasta.empty()) || !bam.empty());
	}

	bool Parse (const int argc, char** const argv) {
		// Record the input command line.
		for (int i = 0; i < argc; i++) cmd += std::string(argv[i]) + " ";
		const struct option long_option[] = {
			{"help", required_argument, NULL, 'h'},
			// i/o
			{"input", required_argument, NULL, 'i'},
			{"fasta", required_argument, NULL, 'f'},
			{"bam", required_argument, NULL, 'b'},
			{"output", required_argument, NULL, 'o'},

			// operation parameters
			{"region", required_argument, NULL, 'r'},
			{"bin", required_argument, NULL, 1},
			{"rle", no_argument, NULL, 2},
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
				case 'r': region = optarg; break;
				case 1: bin = atoi(optarg); break;
				case 2: running_length_encoding = true; break;
				default: std::cerr << "WARNING: Unkonw parameter: " << long_option[option_index].name << std::endl; break;
			}
		}

		if (help) {
			Help(argv[0]);
			return false;
		}

		return CheckArg();
	}
};

#endif
