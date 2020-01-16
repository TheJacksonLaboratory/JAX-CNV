# JAX-CNV: clinical-graded copy number variation detector

## Overview
Here we present JAX-CNV, a newly developed NGS-based CNV algorithm and its performance on WGS data. We focused on deletions and duplications that are >50Kb.

## Dependencies
To run the tool, the following dependencies are required
  * **Linux/Unix**    64-bit
  * **gcc**    version 4.9 or higher
  * **zlib**
  * **autoconf**    version 2.69 or higher

JAX-CNV depends on the following tools, which are already included in JAX-CNV/lib/
  * **fastaq**    [Github](https://github.com/wanpinglee/fastaq/tree/990d69bffe24a2ea2adf823052bddcf25ea71017)
  * **jellyfish-2.2.6**    [Github](https://github.com/gmarcais/Jellyfish/releases/tag/v2.2.6)
  * **htslib**    [Github](https://github.com/samtools/htslib/tree/b8941e42e1962a026ff1f742df1a66c7edddf89c)

## Download and Installation
```Shell
git clone --recursive https://github.com/TheJacksonLaboratory/genome_similarity.git
cd genome_similarity
make
```

## Usages
### Kmer FASTA file preparation
We employ jellyfish to check 25-mer counts and GrabJellyfishKmer to dump a kmer FASTA file.
Please check [jellyfish](https://github.com/gmarcais/Jellyfish/releases/tag/v2.2.6) for more jellyfish options, such as --threads/-t (Number of threads) and --Files/-F (Number files open simultaneously). -s is for Initial hash size and please adjust it for your machine.
```
bin/jellyfish count -m 25 -s <INT> -o <FASTA>.jf [-t <INT> -F <INT>] <FASTA>
bin/JAX-CNV GrabJellyfishKmer --ascii -i <FASTA>.jf -f <FASTA> -o <FASTA>.kmer
```

### Detect CNVs
A sorted BAM, FASTA and Kmer are required. The results will be printed on stdout or use -o to specify an output file.
```
bin/JAX-CNV GetCnvSignal -f <FASTA> -k <FASTA>.kmer -b <BAM> [-o <OUTPUT>]
Rscript --vanilla JaxCNVMerge.R -i <OUTPUT>.bed
```
<OUTPUT>.bed.merge.bed is the final result. JaxCNVMerge.R could be also applied for other tools'' bed files.

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

## For Dcoker users
[Dockerfile](Dockerfile) is provided. Please notice that sudo may be required for docker usages depending on your machine setting.
```
cd JAX-CNV
docker build .
```
JAX-CNV wnd jellyfish will be built on /tools in docker.
Or, Pull docker image from [wanpinglee/jax-cnv](https://cloud.docker.com/repository/docker/wanpinglee/jax-cnv).
```
docker pull wanpinglee/jax-cnv:latest
```
