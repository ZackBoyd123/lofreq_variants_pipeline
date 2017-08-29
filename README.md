# ViralVariantAnalysisPipeline


A Collection of bash and python scripts used for my MSc Projects.

Below is a description each script and an example of it being run on the command line.

**Diversitools2VCF.py** 

A python script which takes the raw output “_entropy” file from diversitools and
converts it to a standard VCF 4.0 output. A file with only SNPs is produced
“.raw.txt”, a file with only indels “.indels.txt” and a stat file are also produced.
The stat file contains percentage of reads within a coverage threshold; number
of genome positions where mutation are within a threshold and number of
mutations between a certain threshold frequency.

Script Should be run as:

*python3 diveritools2VCF.py input_entropy output.txt*

**DiversitoolsStats.py**

A python script which outputs basic stats from a raw diversitools file (_entropy)
to screen. These stats are: total number of bases aligned; total number of
mutations; total number of insertions and total number of deletions. 

Run the script as:

*python3 DiversitoolsStats.py input_entropy*

**Freebayes2VCF.py**

A python script which converts the output from freebayes to a standard VCF 4.0
format. This script will also represent allele frequency as a real value rather than
0.5 or 1.0 as is standard in freebayes. An indel, raw and SNP file is produced.

Run the script as:

*python3 Freebayes2VCF.py inputfile output.txt*

**VarScan2VCF.py**

A python script which converts the output from Varscan to a standard VCF 4.0
format. 

Run the script as:

*python3 VarScan2VCF.py inputfile output.txt*

**VPhaser2VCF.py**

A python script which converts the output from VPhaser2 to a standard VCF 4.0
format.

Run as:

*python3 VPhaser2VCF.py inputfile output.txt*

**VCFDistanceMeasure.py**

A python script which takes two of the converted VCF files, generated from the
“2VCF” scripts- or a Lofreq vcf- and finds similar SNPs between the two. The
output is all mutations, their locations and allele frequencies and a sum of
squares score for each position and a total sum of squares score. Mutations
which are shared are at the top of the file, and those which are unique are
below shared variants.

Run as:

*python3 VCFDistanceMeasure.py -1 firstinput -2 secondinput
-O output.txt*

**VCFStat.py**

A python script which generates stats from a file converted from any of the
“2VCF” scripts or a Lofreq vcf file. Output to screen are the following:

• Total mutations


• Allele frequency of mean mutation

• Allele frequency of median mutations

• Lowest allele frequency reported

• Highest allele frequency reported

Furthermore, the total number of mutations which fall into the following
categories are reported:


• Allele frequency >= 10%


• Allele frequency >= 1% and < 10%


• Allele frequency >= 0.1% and < 1%


• Allele frequency >= 0.01% and < 0.1%


• Allele frequency >= 0.001% and < 0.01%


• Allele frequency < 0.001%


Run the script as:

*python3 VCFStat.py inputfile*

**VCFilter.py**

A python script written to filter variants from a VCF file. This script will work on
any of the files produced by the above”2VCF” scripts and a Lofreq* vcf file.
For example VCFilter.py could be used on the output
of “Diversitools2VCF.py” to only see high frequency mutations occurring above
0.5%. The script is run on the command line, and requires some flags to be set:

• -C: Minimum coverage of mutations.


• -F: Minimum allele frequency of mutations.


• -Q: Minimum quality of mutations.


• --strandbias: Should mutation be observed on both forward and reverse
strand? (yes/no)


• -I: Input file.


• -O: Output file.


Run as:

*python3 VCFilter.py -C 50000 -F 0.1 -Q 0 --strandbias yes -I inputfile -O output file*

**ValVs.sh**

Viral Alignment and Variants (ValVs) is a bash script which essentially
automates a portion of a bioinformatician’s workflow. The script takes an input,
raw fastq file(s) and a reference sequence, trims low quality bases at a
specified value, using TrimGalore
( https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) then aligns
the reads using one of three user specified aligners (BWA,Bowtie2,Tanoti). The
SAM files generated from aligning reads are automatically sorted, indexed and
converted to BAM files. Along with the BAM files generated, a consensus fastq
file and an mpileup file are produced. Users also have the ability to also run one
or multiple variant callers on the BAM files generated previously, choosing
between: Diversitools, VPhaser2, Lofreq and VarScan. The script requires
multiple flags to be set and are as follows:

• [-r]: A reference fasta file.


• [-U]: A raw fastq file with unpaired reads.


• [-1]: First raw fastq file of a paired read sample.


• [-2]: Second raw fastq file of a paired sample.


Note: Either use -U or -1 and -2, not both.


• [-q]: Int Quality score for trim galore to use to trim bases (INT).


• [-l]: Reads below this length will be trimmed from data (INT).


• [-c]: Remove this many bases pairs from 5’ end of reads (INT)


Note: -c only required for paired reads.


• [-a]: Which aligner to use [bowtie2] / [bwa] / [tanoti]


• [-m]: If bowtie2 is specified run in [local] or [end] alignment mode.


• [-v]: Which variant caller to use [diversitools] / [lofreq] / [varscan] /
[vphaser] / [all]


• [-t]: Number of threads to run capable programs on (INT)(Not required).


• [-h]: Get help for script.

Help can be invoked on the script by using only the -h flag.

An example of the script being run on unpaired data would be:

*./ValVs.sh -r ref.fasta -U unpaired.fastq -q 30 -l 50 -a
bowtie2 -m local -v diversitools -t 12*

An example of the script being run on paired data with a different aligner and
variant caller would be:

*./ValVs.sh -r ref.fasta -1 firstpair.fastq -2
secondpair.fastq -q 30 -l 50 -c 14 -a bwa -v vphaser -t 12*




