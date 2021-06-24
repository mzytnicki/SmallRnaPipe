# SmallRnaPipe


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
