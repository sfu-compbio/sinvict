# SiNVICT #

SiNVICT is a tool for the detection of SNVs and indels from cfDNA/ctDNA samples obtained by ultra-deep sequencing.

### Prerequisites ###

* g++ version: 4.8.3
* ~~Boost version: 1_53~~ (no longer required)

----

### Getting SiNVICT ###

To install SiNVICT, first you should fetch it from our git repository, or download one of the corresponding compressed zip/tar.gz packages. After downloading, change the current directory to the source directory sinvict, and run make in the terminal. The sinvict binary will be created, which is ready to use.

```
git clone https://github.com/sfu-compbio/sinvict.git
```

To get all submodules included with SiNVICT as well, the following command can be used:

```
git clone --recursive https://github.com/sfu-compbio/sinvict.git
```
----

### Preprocessing Data for SiNVICT ###

SiNVICT requires the readcount file (See [here](https://github.com/genome/bam-readcount) for a detailed description of the format) to detect the SNVs/indels. However, based on the input file you have, you may follow the steps as described below to obtain the readcount file. SiNVICT is pre-packaged with the tools we have tested for different steps of obtaining the readcount file. You may opt to use any other software as you see fit.

#### Obtaining Readcount Files from different  input files ####

* Input file: **FASTQ**
 1. Trimming FASTQ files (optional, fastq->fastq).  
   We recommend trimming your raw fastq file if it has not be qc'ed and contains barcordes and adapters.  Use one of the available FASTQ trimmers (such as [fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)) to remove any remaining adapters, barcodes, etc. Please refer to its [manual](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_trimmer_usage) for details about running this tool.
 2. Mapping (fastq -> bam).  
   Use one of the available short read mappers to map the reads (e.g. [mrFAST](https://github.com/BilkentCompGen/mrfast), [bwa](https://github.com/lh3/bwa)). Please see their manuals for details.
* Input file: **BAM**
 1. Recalibration and error correction (optional, bam -> bam).  
   We recommend using a recalibration and error correction tool on your BAM file if the initial alignments seem noisy. Use one of the available recalibration/error correction tools (such as [ABRA](https://github.com/mozack/abra)) to improve noisy alignments. Please see their manuals for details.
 2. Obtaining mapping statistics per location (bam -> readcount).  
   Use [bam-readcount](https://github.com/genome/bam-readcount) to gather per location statistics from the alignment files. Please see its manuals for details.

----

### Running SiNVICT ###

The simplest way to run SiNVICT is to use the following command:

```
./sinvict -t input-readcount-dir -o .
```

You may run SiNVICT with modified parameters. For example, to require the calls to have a minimum read depth of 70, you may use the following command.
```
./sinvict --minDepth 70 -t input-readcount-dir -o .
```

There are many options that can be added to the SiNVICT command line. Here are the explanation of these parameteres.

 * --error-rate:  Error Rate for the sequencing technology used (e.g. Illumina, Ion Torrent, ...) (Default: 0.01)
 * --min-depth: Minimum required read depth for a call to be considered reliable. (Default: 100)
 * --left-strand-bias and --right-strand-bias: The strand bias values in the range [leftStrandBias, rightStrandBias] will be considered reliable. This [lsb, rsb] interval has to be between [0,1]. (Defaults: 0.3 and 0.7)
 * --read-end-fraction: Despite the trimming step, calls can be marked as low confidence according to the average position of the base on the reads that support a call. This value should be within range 0-1. (Default: 0.01)
 * --qscore-cutoff: The poisson model used by SiNVICT assigns a QScore to every call. Calls with a QScore below the user defined threshold will be considered low confidence. This value should be in range 0-99. (Default: 95)
 * --tumor-directory-path: The path to the directory where the readcount files are located.
 * --output-directory-path: The path to an empty directory where the output files will be generated.  


----

## Output ##

SiNVICT generates a number of output files that are ordered according to their level of confidence/filtering as named below (each filter eliminates some calls from the previous layer):

 1. Poisson model: **calls_level1.sinvict**
 2. Minimum Read Depth filter: **calls_level2.sinvict**
 3. Strand-bias filter: **calls_level3.sinvict**
 4. Average position of called location among all reads supporting the call: **calls_level4.sinvict**
 5. Signal-to-Noise ratio filter: **calls_level5.sinvict**
 6. Homopolymer Regions filter: **calls_level6.sinvict**

Each tab delimited line in an output file corresponds to a call made by SiNVICT with the following fields:

 1. Chromosome Name
 2. Position
 3. Sample Name (readcount file name)
 4. Reference Base
 5. Read Depth
 6. Mutated Base(s)
 7. Number of reads supporting the mutation
 8. Variant Allele Percentage
 9. Number of reads mapped to the + strand (in format "+:`value`")
 10. Number of reads mapped to the - strand (in format "-: `value`")
 11. Average position of the base on reads as a fraction (0.00 - 1.00)
 12. "Somatic" or "Germline", predicted based on variant allele frequency (do not take as the ground truth)


The following is a sample output from SiNVICT:

```
chrX    66943552        22RV1_49C_10_to_1_10ng  A       3709    G       602     16.2308 +:430   -:172   0.53    Somatic
chrX    66943552        22RV1_49C_10_to_1_1ng   A       1979    G       111     5.60889 +:75    -:36    0.54    Somatic
chrX    66943552        22RV1_49C_10_to_1_2.5ng A       6221    G       803     12.9079 +:567   -:236   0.53    Somatic
chrX    66943552        22RV1_49C_10_to_1_5ng   A       2508    G       265     10.5662 +:183   -:82    0.52    Somatic
chrX    66943552        22RV1_49C_20_to_1_10ng  A       5247    G       304     5.79379 +:167   -:137   0.53    Somatic

```

----

### SiNVICT man page ###

To view the full list of SiNVICT options and their descriptions, please run the following:

```
./sinvict -h
```
----

### Contact ###

For any additional questions/comments/suggestions, please send an email to the following address:

ckockan@sfu.ca

