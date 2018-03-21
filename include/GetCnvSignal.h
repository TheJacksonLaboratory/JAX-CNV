#ifndef _GETCNVSIGNAL_H_
#define _GETCNVSIGNAL_H_

#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <string>
#include <algorithm>

struct SGetCnvSignalCml {
	SGetCnvSignalCml(){}
	SGetCnvSignalCml(const int argc, char** const argv){Parse(argc, argv);}

	bool help = false;

	// i/o
	std::string kmer_table; // -k --kmer
	std::string fasta;      // -f --fasta
	std::string bam;	// -b --bam
	std::string output;     // -o --output

	// operation parameters
	int coverage = 0;	// -c --coverage
	std::string region;     // -r --region
	uint8_t aln_qual = 40;	// -q --aln_qual
	unsigned int minimum_report_size = 45000;	// -m --minimum_report_size
	int bin = 50;		// --bin
	std::string log;	// --log
	float unique_kmer = 0.6;	// --unique_kmer
	float kmer_score = 0.1;	// --kmer_score
	
	// command line
	std::string cmd;

	const char* short_option = "hb:k:f:o:c:r:q:m:";

	// Help list
	const std::string Help (const char* program) const { return
		std::string("USAGE: ") + program + std::string(" -f <FASTA> -k <kmer_table> -b <BAM>\n\n") +
		std::string("	-h --help			Print this help list.\n") +
		std::string("\n") +
		std::string("Input & Output:\n") +
		std::string("	-b --bam <BAM>			Input BAM; required.\n") +
		std::string("	-k --kmer <kmer_table>		Kmer table.\n") +
		std::string("	-f --fasta <FASTA>		FASTA for kmer lookup.\n") +
		std::string("	-o --output <FILE>		Output file.\n") +
		std::string("\n") +
		std::string("Operations:\n") +
		std::string("	-c --coverage <INT>		The expected coverage.\n") +
		std::string("	-r --region chr:begin-end	A target region.\n") +
		std::string("	-q --aln_qual			A mapping quality filter for alignments. [40]\n") +
		std::string("	-m --minimum_report_size	The minimum report SV size. [45K]\n") +
		std::string("	--bin <INT>			Report a result for each # bp. [50]\n") +
		std::string("	--log <FILE>			Log output.\n" +
		std::string("	--unique_kmer <FLOAT>		Require percentage of unique kmer to report a CNV. [0.6]\n") +
		std::string("	--kmer_score <FLOAT>		Score for log2(kmer count) = 2 positions. [0.1]\n"));
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
		} else if (access(bam.c_str(), F_OK) == -1) {
			std::cerr << "ERROR: Cannot open " << bam << std::endl;
			ok = false;
		}

		if (kmer_table.empty()) {
			std::cerr << "ERROR: -k <kmer_table> is required." << std::endl;
			ok = false;
		} else if (access(kmer_table.c_str(), F_OK) == -1) {
			std::cerr << "ERROR: Cannot open " << kmer_table << std::endl;
			ok = false;
		}

		if (fasta.empty()) {
			std::cerr << "ERROR: -f <FASTA> is required." << std::endl;
			ok = false;
		} else if (access(fasta.c_str(), F_OK) == -1) {
			std::cerr << "ERROR: Cannot open " << fasta << std::endl;
			ok = false;
		}

		if (unique_kmer > 1) {
			std::cerr << "ERROR: --unique_kmer <FLOAT> should not larger than 1." << std::endl;
			ok = false;
		}
		if (kmer_score > 1) {
			std::cerr << "ERROR: --kmer_score <FLOAT> should not larger than 1." << std::endl;
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
			{"kmer", required_argument, NULL, 'k'},
			{"fasta", required_argument, NULL, 'f'},
			{"output", required_argument, NULL, 'o'},

			// operation parameters
			{"coverage", required_argument, NULL, 'c'},
			{"region", required_argument, NULL, 'r'},
			{"aln_qual", required_argument, NULL, 'q'},
			{"minimum_report_size", required_argument, NULL, 'm'},
			{"bin", required_argument, NULL, 2},
			{"log", required_argument, NULL, 3},
			{"unique_kmer", required_argument, NULL, 4},
			{"kmer_score", required_argument, NULL, 5},
			{0,0,0,0}
		};
		int option_index = 0;
		int c = -1;
		bool error = false;
		while ((c = getopt_long(argc, argv, short_option, long_option, &option_index)) != -1) {
			std::string tmp;
			switch (c) {
				case 'h': help = true; break;
				case 'k': kmer_table = optarg; break;
				case 'f': fasta = optarg; break;
				case 'b': bam = optarg; break;
				case 'o': output = optarg; break;
				case 'c': coverage = atoi(optarg); break;
				case 'r': region = optarg; break;
				case 'q': aln_qual = atoi(optarg); break;
				case 'm':
					tmp = optarg;
					if (std::all_of(tmp.begin(), tmp.end(), ::isdigit)) { // all digits
						minimum_report_size = atoi(optarg);
					} else {
						if ((tmp.back() == 'K' || tmp.back() == 'k') && (std::all_of(tmp.begin(), tmp.end() - 1, ::isdigit))) {
							minimum_report_size = atoi(tmp.substr(0, tmp.size()-1).c_str()) * 1000;
						} else if ((tmp.back() == 'M' || tmp.back() == 'm') && (std::all_of(tmp.begin(), tmp.end() - 1, ::isdigit))) {
							minimum_report_size = atoi(tmp.substr(0, tmp.size()-1).c_str()) * 1000000;
						} else {
							std::cerr << "Error: Cannot parse " << tmp << ". Please use all digits, or K or M as suffix." << std::endl;
							error = true;
						}
					}
					break;
				case 2: bin = atoi(optarg); break;
				case 3: log = optarg; break;
				case 4: unique_kmer = atof(optarg); break;
				case 5: kmer_score = atof(optarg); break;
				default: std::cerr << "WARNING: Unkonw parameter: " << long_option[option_index].name << std::endl; break;
			}
		}

		if (help) {
			Help(argv[0]);
			return false;
		}

		return error && CheckArg();
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
