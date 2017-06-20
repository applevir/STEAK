# STEAK
Specific Transposable Element Aligner (HERV-K)
====================
STEAK is a tool designed to validate and discover transposable element (TE) and retroviral integrations in a variety of high-throughput sequencing data. By using trimmed reads both reference and non-reference integrations can be marked. Furthermore, if the input HTS data is paired-end, then respective paired-end reads will also be extracted for downstream validations.


Dependencies
-------------
STEAK relies on boost-libraries [http://www.boost.org], OpenMP [http://www.openmp.org], gcc (4.9.2 or higher), bedtools2 [https://github.com/arq5x/bedtools2/], python (2.7.x), and various common unix tools such as sort and awk. 

SAMtools [http://www.htslib.org] or biobambam2 [https://github.com/gt1/biobambam2/releases] can be used to decompress while collating (pairs reads together) BAM files.

Trimmed reads can be processed with an aligner of choice – we recommend a sensitive mapper such as Novoalign [http://www.novocraft.com/support/download/], bwa (mem) [https://sourceforge.net/projects/bio-bwa/files/] or Stampy [http://www.well.ox.ac.uk/stampy].

Using STEAK
================
I. Trimming Chimeric Reads
---------------------------



II. Mapping of trimmed reads & integration detection
-----------------------------------------------------

An aligner of choice can be used to map trimmed reads back to the original host. However, we do recommend a sensitive mapper (see dependencies for options). For our purposes have used Novoalign with parameters that can be found in the example_aligners.txt file.

After the trimmed reads have been mapped to the host, for detection of reference and non-reference integrations the commands found in the integration_detection.txt file can be used to produce the locations of candidate TE integrations.


Questions or help
-----------------
For questions or support contact Cindy Santander (cindy.santander@zoo.ox.ac.uk) or Philippe Gambron (philippe.gambron@stfc.ac.uk)
