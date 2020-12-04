# ampliconFilter v1.0

## What does this app do?
This app runs [AmpliconFilter v1.0](https://github.com/moka-guys/ampliconfilter/tree/v1.0) to remove any alignments which are not the product of an amplicon.


## What are typical use cases for this app?
Amplicons can be formed from pairs of primers that were designed to create an amplicon and from 'super' amplicons, created from primers from neighbouring amplicons (see the --super flag).

Fragments which match the expected amplicons have the primer sequences removed from either end (see the --clipping options). 

Any fragments which do not originate from one of these expected amplicons are filtered into the discarded BAM file.

This tool replaced BAMclipper as BAMclipper was found to be clipping more bases than expected and did not take into account super amplicons. 

## What inputs are required for this app to run?
* BAM (file) - A BAM file
* BEDPE (file) - A paired end BED file containing the primer pairs
* clipping (int) - 0=noclip, 1=softclip (default), 2=hardclip
* super (Boolean) - which allows all primer combinations (chimera/superamplicons)
* mask (Boolean) -  which sequence Hardmasking, default False
* primer distance (int) - string maximum distance from nearest possible primer, default = 0
* tolerance (int) - detection end tolerance, default = 0
* maxbuffer (int) - maximum read buffer size, default (when not stated) = 100000. App default is 1000000.
* extratrim (int) - additional trimming after primer, default = 0

## How does this app work?
Downloads reference genome (hs37d5.fasta-index.tar.gz) and the hardcoded version of the docker file from 001.

The parameters/flag inputs described above are concatenated into a string.

The ampliconfilter docker file is loaded and run (the version of ampliconFilter code applied in this app is listed above).

AmpliconFilter produces 2 BAM files, one containing discarded reads and another containing reads which match expected amplicons. These BAM files are sorted and indexed using samtools (v0.1.19-44428cd)


## What does this app output
* clipped_BAM - sorted BAM file (output to /output)
* clipped_BAM_BAI - index of clipped_BAM (output to /output)
* discarded_BAM - BAM file containing discarded reads (output to /output)
* discarded_BAM_BAI - index of discarded_BAM (output to /output)
* metrics - metrics from ampliconFilter (output to /QC )

## This app was made by Viapath Genome Informatics
