#ifndef _COMMANDLINE_H_
#define _COMMANDLINE_H_

#include <string>

struct SubCommand {
	const unsigned int no_sub_commands = 2;
	const char sub_commands[] = {"CountKmer", "GetCnvSignal"};

	const std::string Help (const char* program) const { return
		std::string("USAGE: ") + program + std::string(" <command> [options]\n\n") +
		std::string("Commands:\n") +
		std::string("\tCountKmer	Report the count of kmer giving Jellyfish database and a FASTA.\n") +
		std::string("\tGetCnvSignal	Report CNV signals such as read depth and kmer count.\n");
	}
}

#endif // _COMMANDLINE_H_
