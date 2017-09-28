#include <string.h>

#include <iostream>
#include <string>

// Self include
#include "CountKmer.h"
#include "GetCnvSignal.h"

struct SubCommand {
	const unsigned int no_sub_commands = 2;
	const std::string sub_commands[2] = {"CountKmer", "GetCnvSignal"};

	const std::string Help (const char* program) const { return
		std::string("\n") +
		std::string("USAGE: ") + program + std::string(" <command> [options]\n\n") +
		std::string("Commands:\n") +
		std::string("\tCountKmer	Report the count of kmer giving Jellyfish database and a FASTA.\n") +
		std::string("\tGetCnvSignal	Report CNV signals such as read depth and kmer count.\n");
	}
};

int main (int argc, char** argv) {
	SubCommand cml_option;
	std::string command;
	if(argc < 2) {
		std::cerr << cml_option.Help(argv[0]) << std::endl;
		return 1;
	} else {
		for (unsigned int i = 0; i < cml_option.no_sub_commands; ++i) {
			if (strcmp(argv[1], cml_option.sub_commands[i].c_str()) == 0) {
				// Get the valid subcommand
				command = cml_option.sub_commands[i];
				break;
			}
		}
	}

	if (command.empty()) { // The given command is not in the list
		std::cerr << "ERROR: The given command (" << argv[1] << ") is not valid." << std::endl;
		std::cerr << cml_option.Help(argv[0]) << std::endl;
		return 1;
	} else { // Get the valid subcommand
		if (command == "CountKmer") {
			CountKmer count_kmer(argc - 1, argv + 1);
			count_kmer.Run();
		} else if (command == "GetCnvSignal") {
			GetCnvSignal get_cnv_signal(argc - 1, argv + 1);
			get_cnv_signal.Run();
		}
	}

	return 0;
}
