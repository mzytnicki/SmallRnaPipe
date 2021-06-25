# SmallRnaPipe

## Why are we doing this pipeline ?

This pipeline use some programms which aren't use by any other pipeline which analyses Small Rnas. These programms ([srnaMapper]( https://github.com/mzytnicki/srnaMapper) and [mmquant](https://bitbucket.org/mzytnicki/multi-mapping-counter/src/master/)) make that this pipeline is very useful for the reproductibily of the analyses of Small Rnas. 

## Arguments



|  Option  |  Parameter(s)  |  Description  |  Requirement  |
|---   |:-:   |:-:   |--:  |
|  `--reads`  |  `fastq1 ...`  |  Input `fastq` file(s) |  Required  |
|  `--genome`  |  `genome.fa`  | A FA genome file  |  Required  |
|  `--index`  |  `directory/prefix`  |  Input genome index directory built by `bwa`  |  Optional  |
|  `--annotation`  |  `annotation.(gtf/gff)`  |  Input reference annotation file  |  Required  | 


## Workflow


|  Programm  |  Action  |  Inputs  |  Command  |  Outputs  |
|:-:  |:-:  |:-:  |:-:  |:-:  |
|  [FASTQC](https://github.com/s-andrews/FastQC)  |  Control the quality of fasta files in input to have a detailed report about quality  |  Fastq file (reads)  |  `fastqc ech1.fastq ...`  |  HTML with a resume and a ZIP file for the MultiQC  |
|  [trim_galore](https://github.com/FelixKrueger/TrimGalore)  |  



## Usage : a typical command line 

`nextflow run main_script1.nf --reads ech1.fastq --annotation ann1.gff3 --index index1 --genome genome1.fa`

**Be careful :** the index must be a index built by `bwa` so index1 must be the prefix that you have written during the index construction by `bwa`  
