# genome_similarity

## Overview
Clinical Laboratory Improvement Amendments- (CLIA-)graded copy number variation detector.

## Dependencies
To run the tool, the following dependencies are required
  * **Linux/Unix**    64-bit
  * **gcc**    version 4.9 or higher
  * **zlib**
  * **autoconf**    version 2.69 or higher

genome_similarity depends on the following tools, which are already included in the genome_similarity/lib/
  * **fastaq**    [Github](https://github.com/wanpinglee/fastaq/tree/990d69bffe24a2ea2adf823052bddcf25ea71017)
  * **jellyfish-2.2.6**    [Github](https://github.com/gmarcais/Jellyfish/releases/tag/v2.2.6)
  * **htslib**    [Github](https://github.com/samtools/htslib/tree/b8941e42e1962a026ff1f742df1a66c7edddf89c)

## Download and Installation
```Shell
git clone https://github.com/TheJacksonLaboratory/genome_similarity.git
cd genome_similarity
make
```

## Usages
### Kmer FASTA file preparation
We employ jellyfish to check 25-mer counts and GrabJellyfishKmer to dump a kmer FASTA file.
Please check [jellyfish](https://github.com/gmarcais/Jellyfish/releases/tag/v2.2.6) for more jellyfish options, such as --threads/-t (Number of threads) and --Files/-F (Number files open simultaneously). -s is for Initial hash size and please adjust it for your machine.
```
bin/jellyfish count -m 25 -s <INT> -o <FASTA>.jf [-t <INT> -F <INT>] <FASTA>
bin/clia_cnv GrabJellyfishKmer --ascii -i <FASTA>.jf -f <FASTA> -o <FASTA>.kmer
```

### Detect CNVs
A sorted BAM, FASTA and Kmer are required. The results will be printed on stdout or use -o to specify an output file.
```
bin/clia_cnv GetCnvSignal -f <FASTA> -k <FASTA>.kmer -b <BAM> [-o <OUTPUT>]
```
The complete command line options are:
```
USAGE: GetCnvSignal -f <FASTA> -k <kmer_table> -b <BAM>

        -h --help                       Print this help list.

Input & Output:
        -b --bam <BAM>                  Input BAM; required.
        -k --kmer <kmer_table>          Kmer table.
        -f --fasta <FASTA>              FASTA for kmer lookup.
        -o --output <FILE>              Output file.

Operations:
        -c --coverage <INT>             The expected coverage.
        -r --region chr:begin-end       A target region.
        -q --aln_qual                   A mapping quality filter for alignments. [40]
        --bin <INT>                     Report a result for each # bp. [50]
        --log <FILE>                    Log output.
        --unique_kmer <FLOAT>           Require percentage of unique kmer to report a CNV. [0.6]
        --kmer_score <FLOAT>            Score for log2(kmer count) = 2 positions. [0.1]
```


