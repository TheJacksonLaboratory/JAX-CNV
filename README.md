# genome_similarity

## Overview
genome_similarity is a tool for fast quantifying whole genome-wide DNA sequential similarity. 

## Dependencies
To run the tool, the following dependencies are required
  * **Linux/Unix**    64-bit
  * **gcc**    version 4.9 or higher
  * **autoconf**    version 2.69 or higher

genome_similarity depends on the following tools, which are already included in the genome_similarity/include/
  * **fastaq**    [fastaq Github](https://github.com/wanpinglee/fastaq)
  * **jellyfish**    [jellyfish Github release](https://github.com/gmarcais/Jellyfish/releases)
  * **htslib**    [htslib Github](https://github.com/samtools/htslib)

## Download and Installation
```Shell
git clone git@github.com:wanpinglee/genome_similarity.git
cd genome_similarity
make
```

## Command line options
### jellyfish
Please refer to [jellyfish](https://github.com/gmarcais/Jellyfish) for illustartion. 
### GetKmerCount
It collects the k-mer counts for each DNA position. The command can be used as
```
genome_similarity/bin/GetKmerCount -i <jellyfish_db> -f <FASTA>
```
And the command line options are:
```
USAGE: ./GetKmerCount -i <jellyfish_db> -f <FASTA>

	-h --help			Print this help list.

Input & Output:
	-i --input <jellyfish_db>	Jellyfish created count database.
	-f --fasta <FASTA>		FASTA for kmer lookup.
	-b --bam <BAM>			Input BAM.
	-o --output <FILE>		Output file.

Operations:
	-r --region chr:begin-end	Specify a target region.
	--bin <INT>			Report a result for each # bp. [1]
	--rle				Ouput by running length encoding.
```


