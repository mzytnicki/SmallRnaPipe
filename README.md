# SmallRnaPipe

## Why are we doing this pipeline ?

This pipeline use some programms which aren't use by any other pipeline which analyses Small Rnas. The existing pipelines analyses only microRNAs, a category of smallRNAs, and they don't consider in a satisfactory way the sequences aligning in a multiple way on the genome. These programms ([srnaMapper]( https://github.com/mzytnicki/srnaMapper) and [mmquant](https://bitbucket.org/mzytnicki/multi-mapping-counter/src/master/)) make that this pipeline is very useful for the reproductibily of the analyses of smallRnas. 

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
|  [trim_galore](https://github.com/FelixKrueger/TrimGalore)  |  Trimmming the adaptators on your reads  |  Fastq file (reads)  |  `trim_galore ech1.fastq ...`  | un `*report.txt` for the multiqc and un fastq trimmed for the next step  |  
|  [prinseq-lite](http://prinseq.sourceforge.net/)  |  Sort the reads with the "good" reads with great complexity and the "bad" reads with a low complexity  |  Fastq trimmed file (reads)  |  `prinseq-lite.pl -fastq ech1.trimmed.fastq ... -lc_method dust -lc_threshold 7`  |  2 fastq : one with the bad reads and un other with the good reads for the next step  |
|  [BWA](https://github.com/lh3/bwa)  | If you haven't a bwa index for your genome, this process will do that.  |  `bwa index -p $prefix $genome.fa`  |  5 five which are your index.  |
|  [srnaMapper](https://github.com/mzytnicki/srnaMapper)  |  Map the Fastq on the reference genome  |  the Index, the clens reads and the prefix  |  `srnaMapper -r $reads_clean -g direction/prefix_of_index -o $prefix.sam`  | A sam file for each reads for the next step  |
| [mmquant](https://bitbucket.org/mzytnicki/multi-mapping-counter/src/master/)  |  Quantifie the expression of the sRnas  |  the annotation file(in GFF or GFF3 pr GTF)and the bam/sam  |  `mmquant -a $annotation_file -r bam/sam -o prefix.tsv`  |  A table where, for each annotation, you have un number which represents the number of times the annotation has been spotted in the BAM  |
|  [MultiQC](https://multiqc.info/)  |  A great representation of all datas and results |  Fastqc report, trimming report  | `multiqc --config $config` (we use a specific config file because we want a specific display order  |  a HTML report  |


## Usage : a typical command line 

`nextflow run main_script1.nf --reads ech1.fastq --annotation ann1.gff3 --index index1 --genome genome1.fa`

**Be careful :** the index must be a index built by `bwa` so index1 must be the prefix that you have written during the index construction by `bwa`  
