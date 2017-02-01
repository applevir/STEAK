# STEAK
Specific Transposable Element Aligner (HERV-K)
====================
STEAK is a tool designed to validate and discover transposable element (TE) and retroviral integrations in a variety of high-throughput sequencing data. By using trimmed reads both reference and non-reference integrations can be marked. Furthermore, if the input HTS data is paired-end, then respective paired-end reads will also be extracted for downstream validations.

STEAK consists of a multi-threaded software written in C/C++ and a set of scripts for downstream analysis. It accepts input SAM files or the stdout produced by SAMtools suite if you are inputing BAM (e.g. samtools view). It searches for TEs using a transposable element reference provided by the user in FASTA format. STEAK then takes the 5’ and 3’ of the TE reference, creates reverse complements, to produce four TE baits to search out matches between the transposable element of interest and HTS reads.

Although STEAK was programmed to accept SAM for detection on already reference mapped data, it can also detect on unmapped projects as long as the sequences are in SAM format (e.g. picard-tools FastqToSam).

It can produce the paired-end (PE) or single-end (SE) outputs in FASTQ format along with two tab-delimited text files; one file with full-length reads and the other with trimmed reads and corresponding mate sequence. These files also contain previous mapping information, a useful attribute when processing pre-mapped reads like those found in public data. If the paired-end option has been specified in the parameters of STEAK, the files give information on the mates of the trimmed reads as well.

Dependencies
-------------
STEAK relies on boost-libraries [http://www.boost.org], OpenMP [http://www.openmp.org], gcc (4.9.2 or higher), bedtools2 [https://github.com/arq5x/bedtools2/], python (2.7.x), and various common unix tools such as sort and awk. 

SAMtools [http://www.htslib.org] or Bamtools [https://github.com/pezmaster31/bamtools] can be used to decompress BAM files. Picard Tools [https://broadinstitute.github.io/picard/] can be used to convert FASTQ raw reads to SAM format.

Trimmed reads can be processed with an aligner of choice – we recommend a sensitive mapper such as Novoalign [http://www.novocraft.com/support/download/], bwa (mem) [https://sourceforge.net/projects/bio-bwa/files/] or Stampy [http://www.well.ox.ac.uk/stampy].

Using STEAK
================
I. Trimming Chimeric Reads
---------------------------

Here the input HTS file (BAM or SAM) is parsed for chimeric reads. Reads that partially match the TE reference are trimmed of the TE and only host is left. Parameters can be adjusted to fine-tune chimeric read detection. 

If you wish to parse HTS PE reads and keep relevant mate-pair information the SAM or BAM must be sorted in pairs (e.g. this can usually be achieved sorting by name).  Otherwise to treat data as SE, SAM or BAM can be in any order. For example usage and execution, please see or use steak.sh

Usage: `steak <paired | unpaired> <sam> <cigar_ceiling> <smith_waterman_threshold> <minimal_trimmed_length> <reference_file> <bait_length>`

	mode						Pass input HTS file in PE (paired) or SE (unpaired) mode
    sam         				Name of SAM file
    cigar_ceiling      			Maximum percentage of ‘M’s in the CIGAR value of a read. Only reads of this percentage or lower will be considered (e.g. CIGAR of 92M2I6M is 0.98)
    smith_waterman_threshold    Percentage of matching needed between TE reference and TE matching part of the read (e.g. 0.90 in a TE match of 20 bp allows for 2 mismatches).
	minimal_trimmed_length   	The minimum length that a trimmed read (host flank) can be. Trimmed reads smaller than this will be ignored.
    reference_file         		FASTA with sequence of TE/virus of interest.
    bait_length        			Length of the TE reference that will be searched for each read (e.g. 40 means that the first 40 and last 40 bp of the TE reference will be searched for in each read).

    
If in paired mode, STEAK will produce four outputs in total. Two FASTQ files: the first with the trimmed reads and the second with the mates of the trimmed reads. It will also output two tab delimited text files with paired read information: a file with the original chimeric sequences with mate sequences and a file with the trimmed sequences and corresponding mate sequences. 

In unpaired mode, STEAK will produce the same outputs minus the second FASTQ file.


II. Mapping of trimmed reads & integration detection
-----------------------------------------------------

An aligner of choice can be used to map trimmed reads back to the original host. However, we do recommend a sensitive mapper (see dependencies for options). For our purposes have used Novoalign with parameters that can be found in the example_aligners.txt file.

After the trimmed reads have been mapped to the host, for detection of reference and non-reference integrations the commands found in the integration_detection.txt file can be used to produce the locations of candidate TE integrations.


Questions or help
-----------------
For questions or support contact Cindy Santander (cindy.santander@zoo.ox.ac.uk) or Philippe Gambron (philippe.gambron@stfc.ac.uk)
