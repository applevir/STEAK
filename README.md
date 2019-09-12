# STEAK
Specific Transposable Element Aligner (HERV-K)
====================
STEAK is a tool designed to validate and discover transposable element (TE) and retroviral integrations in a variety of high-throughput sequencing data. By using trimmed reads both reference and non-reference integrations can be marked. Furthermore, if the input HTS data is paired-end, then respective paired-end reads will also be extracted for downstream validations.


Dependencies
-------------
STEAK relies on boost-libraries [http://www.boost.org], MPI [https://www.open-mpi.org/software/ompi/v2.1/], gcc (4.9.2 or higher), bedtools2 [https://github.com/arq5x/bedtools2/], python (2.7.x), and various common unix tools such as sort and awk. 

SAMtools [http://www.htslib.org] or biobambam2 [https://github.com/gt1/biobambam2/releases] can be used to decompress while collating (paired reads together) BAM files.

Trimmed reads can be processed with an aligner of choice - we recommend a sensitive mapper such as Novoalign [http://www.novocraft.com/support/download/], bwa (mem) [https://sourceforge.net/projects/bio-bwa/files/] or Stampy [http://www.well.ox.ac.uk/stampy].

Compilation
-------------
`cd STEAK/`

`make`

If it does not compile, make sure you indeed have boost-libraries and MPI. Moreover, make sure to be consistent in the compiler used for MPI and when you compile STEAK.

Using STEAK
================
I. Detecting Chimeric Reads
---------------------------
Here the input HTS file is parsed for chimeric reads. 

TE and retrovirus references should be in FASTA format. At this time, STEAK accepts only **one TE sequence per FASTA**. 

Please make sure to check that your TE and retrovirus sequences are in **upper case**. STEAK does not accept lower case FASTA sequences which are meant for genome masking.

STEAK accepts FASTQ, SAM, and piped BAM/CRAM as HTS inputs. 

Reads that partially match the TE or retrovirus reference are trimmed of the match portion. The host flank and TE are output into respective FASTQs. 

Parameters can be adjusted to fine-tune chimeric read detection. 

If you have HTS files that are paired-end, the paired reads must be collated if you wish to retain respective mate information (this is useful for guided detection). For a SAM/BAM file that means its must be **sorted by name or collated in pairs** so that a read is followed by its mate. Paired end FASTQs (e.g. `reads_1.fq & reads_2.fq`) must be in a single interleaved file. 

Otherwise to treat data in an unpaired fashion, or as singletons, your reads can can be in any order.

Usage: `steak options [default values]:`
         
            --input               Input NGS file  (mandatory without the --pipe option)
            --output              Output file (mandatory with the  --pipe option)
            --pipe                Takes the input from a pipe. Works only with BAM/CRAM and as a single process.
            --paired              Input NGS file(s) are paired-end reads.
            --TE-reference        FASTA with sequence of TE/virus of interest.
            --alignment-quality   Maximum proportion of 'M's in the CIGAR value of a read. 
                                  Only reads of this proportion or lower will be considered.   [.99]
            --match-quality       Proportion of the match needed between the TE reference and the read [.95]
            --transposon-length   Length of the TE reference that will be searched for within each read. [15]
            --host-length         The minimum length that a trimmed read (host flank) can be. 
                                  Trimmed reads smaller than this will be ignored. [25]
            --aligned             Uses an aligned genome. Requires a BAM or SAM file. 
                                  Without this option, STEAK will consider that the input consists 
                                  of a fasta or fastq file and the alignment-quality parameter will be ignored.
            


II. Integration Detection Module
-----------------------------------------------------

An aligner of choice can be used to map trimmed reads back to the original host. However, we do recommend a sensitive mapper (see dependencies for options). For our purposes, we have used Novoalign. This is an example of a basic mapping; you can adjust Novoalign parameters for your respective data:

Create Novoindex for host reference:

`novoindex hg19reference hg19.fasta`

Novoalign parameters for Trimmed Read Detection:

`novoalign -d hg19reference -o SAM -f trimmed_reads.1.fastq > trimmed_read_detection.sam`

Novoalign parameters for Guided Detection:

`novoalign -d hg19reference -o SAM -o FullNW -f reads.1.fastq reads.2.fastq > guided_detection.sam`

After the trimmed reads have been mapped to the host, for detection of reference and non-reference integrations the following commands can be used to produce the locations of candidate TE integrations.`

Convert BAM to BED:

`bedtools bamtobed -i  trimmed_read_detection.bam > mapped_trims.bed`

Mark number of mapped reads around known integration sites (works for reference and/or non-reference). This is using a 100 bp window around an integration site:

`bedtools window -w 100 -c -a list_of_known_integration_sites.bed -b mapped_trims.bed > readCount_for_integrations.bed`

Search for de novo integrations:

`bedtools window -w 2000 -v -a mapped_trims.bed -b list_of_known_integration_sites.bed | bedtools merge -c 4 -o count -d 1000 -i - > novel_integrations.bed`

Filter de novo candidates with minimum threshold (e.g. "$4 >= readCount"). Here we have used a minimum of 5 reads at least 2000 bp away from known integration sites:

`bedtools window -w 2000 -v -a mapped_trims.bed -b list_of_known_integration_sites.bed| bedtools merge -c 4 -o count -d 1000 -i - | awk '{if ($4 >=5) print($0);}'- > filtered_novel_ints.bed`



III. Outputs
-----------------------------------------------------
STEAK produces FASTQs and a tab-delimited file with additional read information.

If you chose to detect for chimeric reads in paired HTS data you will have outputs of trimmed reads and their corresponding mates in paired FASTQ files which can be used in the Guided Detection module.

All STEAK processing produces a FASTQ with the TE/viral match that has been trimmed off and an read information file.

The read information tab-delimited file is organised in the following manner:

| Column Number | Read Information | Example |
| :------------ | :------------    | :-------|
| 1  | Read name| `AWESOMEREAD`|
| 2  | TE/virus reference| `MOBILELEMENT`|
| 3  | Host read length |  `61` |
| 4  | TE/Viral match read length| `39` |
| 5  | Match proportion with TE reference | `1.00` |
| 6  | Read sequence: Upper case for host flank and lower case for TE portion| `ccctaagacccttttagtcagtgtggaaaatctctagcaCTCTCTCCAGGGGCTT...`|
| 7  | Mate sequence: Upper case for host flank and lower case for TE portion| `TCAGACTGCCTCCAGTCTGCCAACCTCACCACCTCCTGCCCCACCTCTGGCCTGCAGACA...`|
| 8  | Read mapping info : chr, pos, MAPQ, CIGAR | `1\|223662700\|60\|39S61M `|
| 9  | Mate mapping info: chr, pos, MAPQ, CIGAR  | `1\|223662700\|60\|39S61M` |

Example usage and expected results can be found in the test_dataset.

STEAK is fully parallelised for SAM and FASTQ processing with MPI. Piped BAMs and CRAMs can be handled in a multithreaded fashion using OpenMP. To take advantage of this, see usage under Tips & FAQ.


Tips & FAQ
================
**1. How do I make an interleaved FASTQ from paired FASTQs?**

Paired FASTQ files must be merged into a single interleaved file. If your FASTQ pairs are in two separate files, you can merge them into an interleaved file with this command:

`cat awesome_genome.1.fq awesome_genome.2.fq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > awesome_genome.interleaved.fastq`

Alternatively, if you have a single FASTQ file with pairs but it is not collated, you can use this version of the command:

`cat awesome_genome.unsorted.fq  | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > awesome_genome.interleaved.fastq`


**2. I have BAM/SAM files that are a mix of paired-end and single-end reads. Which mode should I use?**

If you want to use the guided detection you should separate the paired-reads from the singletons and run paired-read CRAM/BAM/SAM with the paired mode and the singletons in the unpaired mode.

Otherwise, you can also run your CRAM/BAM/SAM in the unpaired mode and output only trimmed reads without corresponding mate information.

**3. I would like to use the guided detection but I have a huge BAM file and I cannot be bothered to name sort or name sort + decompress. What should I do?**

STEAK accepts piped BAMs but for paired mode they must be collated by name. We have found biobambam2's bamcollate2 very useful in that you can pipe collated reads to STEAK. An example of this would be:

`bamcollate2 filename=awesome_genome.messy.bam outputformat=sam | steak --pipe --TE-reference MOBILELEMENT.fasta --paired --aligned --output awesome_genome_mobi`

However, we do encourage you make sure that your reads are truly collated in pairs to avoid errors from STEAK.

**4. How do I make STEAK process my SAM/FASTQ in a parallel fashion?**

Make sure you have MPI first. Then, for example, if you were to run with 200 cores:

`mpirun -n 200 steak --input ALL_DA_GENOMES.sam --TE-reference MOBILELEMENT.fasta --paired --aligned`

Note: mpirun is **only** for SAMs and FASTQs. Compressed files (BAM/CRAM) can only be read through pipe with one process. 

If you wish to handle the BAM/CRAM in a multithread fashion you can do the following:

`export $OMP_NUM_THREADS=4`

`bamcollate2 filename=awesome_genome.bam outputformat=sam | steak --pipe --TE-reference MOBILELEMENT.fasta --paired --aligned --output awesome_genome_mobi`

**5. Where can I get TE annotations for my genome of interest?**

If your HTS data corresponds to the genome of a well annotated species, then you can find reference TE annotations on RepeatMasker(http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html). The *.fa.out(.gz) files for each build are tabular files with the respective TE information (chromosome, position, strand, name, etc).

If you want to search for known non-reference integrations, it's best to trust the respective literature for your TE of interest.

For an updated database of human retrotransposons see dbRIP: http://dbrip.brocku.ca

**6. How can I find which read causes a problem?**

The name of the read being processed can be displayed. This allows the easy identification of the read that causes a crash but can create a large amount of data. This can be achieved by recompling STEAK to enable its debug mode and subsequently running it normally. For clarity, it is best to use this feature with a single process and a single thread.

`make clean`

`make debug`

Questions or help
-----------------
For questions or support contact Cindy Santander (cindy.santander@zoo.ox.ac.uk) or Philippe Gambron (philippe.gambron@stfc.ac.uk)


 
Copyright (C) 2015 Cindy Santander, Philippe Gambron
Software distributed under the GNU General Public License version 3 
(see file GPLv3.txt or https://www.gnu.org/licenses/gpl-3.0.html)

