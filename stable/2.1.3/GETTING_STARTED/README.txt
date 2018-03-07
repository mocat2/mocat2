##################################
#    MOCAT v2.3.1 README FILE    #
##################################


Content:
- Note
- Release note
- Installation
- Requirements
- Running on SGE or PBS
- Guides and examples
- Quick Guide
- Included software
- CHANGELOG


==========
   Note
==========
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. All external software
used by MCOAT are (c) repesctive author.


==================
   Release note
==================
MOCAT 1 and 2 were developed by Jens Roat Kultima during his time
at EMBL. Unfortunately, he can no longer actively develop MOCAT nor
provide significant help beyond bug fixes. Sincere apologizes for any
inconvenience this may cause and we still hope you can use MOCAT!


==================
   INSTALLATION
==================

=== EASY ===
To setup MOCAT, simply execute these two commands from the directory
you unpacked MOCAT (one up form this).

1. > ./setup.MOCAT2.pl
2. > source ~/.bashrc; source ~/.bash_profile
3. DONE! You can now use MOCAT!
   Check out the section QUICK GUIDE below how to setup a new project!

   Execute MOCAT by running:
   > runMOCAT.sh
   or:
   > MOCAT.pl
   
=== What is the difference between runMOCAT.sh and MOCAT.pl? ===
runMOCAT.sh is a wrapper script for MOCAT.pl that executes, using
custom scripts files, egneral processing pipelines. These pipelines
include generating mOTU and taxonomic profiles. This wrapper also
sends you an email once all the jobs have finished. In short, it's
much easier to use than the original MOCAT.pl script. But of course,
you can still run MOCAT.pl to process samples the way you wish.

=== DETAILS ===
This step will download and install required databases.
The largest file to download is the custom made
human genome 19 database, which has a download size of 6.5GB.

If you wish to download the simulated metagenome, the size
of these files are 3.5GB in total.

After running this step, all you need to do is
start a new session of UNIX (to load needed variables),
and then you should be able to execute the main
script MOCAT.pl from within any direcotry.

ADDITIONAL NOTE FOR USING USEARCH AND METAGENEMARK!
If you want to run gene prediction using MetaGeneMark, or screen
a custom fasta file for contaminants, you have to download
and extract MetaGeneMark and Usearch, respectively, into the
/ext/metagenemark and /ext/usearch directories respectively.

You can download MetaGeneMark from
http://exon.gatech.edu/GeneMark/license_download.cgi

And download Usearch from
Usearch from http://www.drive5.com/usearch/nonprofit_form.html.

NOTE: MOCAT supports Usearch v5 and v6. We recommend v5, because
the memory restrictions for 32bit systems affect input files to
a higher degree in version 6.

This additional step is required because you have to get a custom licence.

IMPORTANT!
To be able to reproduce the results made in the MOCAT article,
you need to download both MetaGeneMark and Usearch.


==================
   Requirements
==================
To run MOCAT you need:
- UNIX 64 bit system
- Perl 5.8.8+

We also recommend:
- SGE, PBS or LSF queuing system

Required CPUs and hard disk space varies
greatly on your needs, but we recommend
a hard disk of 500+ GB and a 8 core CPU,
preferably a UNIX cluster with much more
CPUs and hard disk space.


================================
   Running on SGE, PBS or LSF
================================
FAQ:
Q1. Why doesn't my jobs start?
A1: Try different parallel environment settings, see below.
Also, procesing steps requires 1-8 CPUs depending on the step.
If you have less that 8 CPUs available your job will not start.
You can change this by forcing a lower number of CPUs to be used
with the -p flag when executing MOCAT.

MOCAT has been designed to run on SGE primarily.
The system has also been tested and run on PBS and LSF.

However, the setup of these queuing system is often very
different and it is difficult to design a system that will work for all
permutations of setups of these queuing system.

We have designed the system to have as few user settings
as possible and currently, the only changable options are:
MOCAT_qsub_system        : SGE [SGE,PBS,LSF,none]
MOCAT_SGE_qsub_add_param : [-l mem_free=6G -l h_vmem=6G] 
MOCAT_PBS_qsub_add_param : [-l select=mem=6gb]
MOCAT_LSF_qsub_add_param : [-l select=mem=6gb]
MOCAT_LSF_queue          : []
MOCAT_LSF_memory_limit   : []
MOCAT_SGE_parallell_env  : smp [smp,mpi,make,qstate,-or other setting on your system-]

And these can be changed in the CONFIG file.
Specifically YOU MAY MOST LIKELY HAVE TO CHANGE THE
MOCAT_SGE_parallell_env setting (which only affects SGE systems),
this is done during the installation.

The additional paremeter can be specified if your queuing system has
specific limitations or other constraints that need to be addressed when submitting jobs.

To change number of CPUs required you can use the -p option when running MOCAT.
By default different processing steps requires between 1 and 8 CPUs and memory usage depends
on the size of your metagenomes.

If you have trouble you can also create all the files and submit the jobs manually by
specifying the -x option when running MOCAT.


=========================
   Guides and Examples
=========================
The MOCAT_user_manual.pdf is the user manual in PDF.

If you chose to download the datasets in the article,
they will be located in the /article_datasets folder


=================
   QUICK GUIDE
=================
WE STRONGLY RECOMMEND to execute
a shell script that can run many of the standard
MOCAT jobs for you easy, just execute:
> runMOCAT.sh

Or, if you want to have more fine tuning, use:
> MOCAT.pl

QUICK NOTES!
1. You need to copy the MOCAT.cfg config file to the project folder!
    - Check that this file is correct for you specific needs!
2. You need a sample file containing names of folders, in which you
have the FastQ files

For a complete and extensive manual
how to use MOCAT, execute 'MOCAT.pl -man'

A new project should have the following
structure (if pair-end reads):

CONFIG and SAMPLE file:
| PROJECT_DIR
| PROJECT_DIR/MOCAT.cfg
| PROJECT_DIR/sample_file

FASTQ files:
| PROJECT_DIR/SAMPLE_1
| PROJECT_DIR/SAMPLE_1/sample_file_lane1.1.fq.gz
| PROJECT_DIR/SAMPLE_1/sample_file_lane1.2.fq.gz
| PROJECT_DIR/SAMPLE_1/sample_file_lane2.1.fq.gz
| PROJECT_DIR/SAMPLE_1/sample_file_lane2.2.fq.gz
| ...

IMPORTANT NOTE FOR PAIRED-END READ SAMPLE FILES:
The sample files must end with .1.fq or .1.fq.gz
for the 1st pair of paired end reads, and .2.fq or
.2.fq.gz for the second paired end reads.
You can of course have as many samples, and lanes per
samples as you want. 

IF THE READS ARE SINGLE END:
Th sample files have to end with .fq or .dq.gz.

THE CONFIG FILE:
The MOCAT.cfg file is required in the base project
direcotry. This file can be copied from this folder.

THE SAMPLE FILE:
The sample_file can be any file (anywhere),
containing the names of the samples (subfolders) to analyze,
and is specified when running MOCAT.pl each time.
The structure of the sample file is:

| SAMPLE_1
| SAMPLE_2
| ...


=======================
   INCLUDED SOFTWARE
=======================
2bwt-builder
bwa
faSomeRecords
fastQC
fastx_quality_stats
hmmsearch
msamtools
pigz
prodigal
psort
samtools
soap.coverage
soap1.05
soap1.06
soap2.21
sqlite3
unpigz


=================
    CHANGELOG
=================

=================
  Version 2.1.3
=================
- Major bugfixes
- Support for BWA
- Support for SLURM
- v1.5.3: added support for using pair and single files in raw data input folders
- v1.5.4: fixed major bugs in new RTF using 3 files
- v1.5.6: t18->t19; removed unnecessary code in the NCBI profiling steps
- v1.5.7: Support for pair.1 pair.2 single in screen fasta file
- v2.1.0: Added support for SLURM queue and BWA for mapping
- v2.1.1: Updated version of fastq_trim_filter_v5_EMBL, fixing FastX bug

=================
   Version 1.3
=================
- Major bugfixes
- Updated included software:
  - BWA      version 0.7.5a-r405
  - samtools version 0.1.19-44428cd
- R dependencies removed
- Better LSF support

=================
   Version 1.2
=================
- Major bugfixes
- Supports generating taxonomic and mOTU profiles
- Supports extracting marker genes using fetchMG
- Supports USearch v6 (but we still strongly recommend using v5)
- Some functionality depends on R

=================
   Version 1.1
=================
- Major bugfixes

=================
   Version 1.0
=================
- Initial release



- end of file -



