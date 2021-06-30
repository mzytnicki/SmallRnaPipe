/
* Parse and check parameters
*/

reads = params.containsKey('reads') ? params.reads : ''
index = params.containsKey('index') ? params.index : ''
genome = params.containsKey('genome') ? params.genome : ''
annotation = params.containsKey('annotation') ? params.annotation : ''



log.info """\


					 ===================================


				                S m a l l R n a P i p e    


				         ===================================
         """



/*
*Control step
*/


error =''

if (!reads) error += "----- No reads provided"


if (!annotation) error += "----- No annotation provided"


if (!genome) error += "----- No genome provided"

if (error) exit 1, error


	
/*
* Channel
*/

Channel
	.fromPath(reads)
	.map { path ->
	filename= path.getName()
	return [path,filename]
}
	.into{ 
	reads_ch_to_fastqc 
	reads_ch_to_trimming
}

Channel.fromPath(genome)
        .map { path ->
        filename= path.getName()
        return [filename,path]
}
        .into{
        genome_to_index
        genome_to_maps
}

Channel.fromPath(annotation)
	.set{ annotation_to_mmquant
	}


/*
* Process to control quality of reads with FASTQC
*/

process control_quality_with_fastqc {

	
	input :
	tuple path(reads), val(prefix)  from reads_ch_to_fastqc

	output:
	path '*_fastqc.zip' into fastqc_to_report
	
	script:
	"""
	fastqc $reads
	"""
}


/*
* Process to trim primers with trim_galore
*/


process trimming {

        publishDir "$baseDir/results/trimming" , mode: 'copy'

        input :
        tuple path(reads), val(prefix)  from reads_ch_to_trimming

        output :
	path "*trimming_report.txt" into trim_to_report
	tuple val(prefix),  path("*_trimmed.fq") into trim_read_to_prinseq

	script:
	"""
	trim_galore $reads --basename $prefix
	"""

}


	

process cleaning_reads_with_low_complexity {

	publishDir "$baseDir/results/prinseq" , mode: 'copy'

	input:

	tuple val(prefix), path(reads_trim) from trim_read_to_prinseq

	output:
	tuple val(prefix), path("*good*.fastq") into read_clean_to_maps
	tuple val(prefix), path("*good*.fastq") into lol
	script:
	"""
	prinseq-lite.pl -fastq $reads_trim -lc_method dust -lc_threshold 7
	"""
}

if (!index){

	process bwa_index_to_create {

		publishDir "$baseDir/results/index_with_bwa" , mode: 'copy'


		input:
		tuple val(prefix), path(reads) from genome_to_index

		output:
		path '*.amb' 
		path '*.sa' 
		path '*.pac' 
		path '*.bwt' 
		path '*.ann' 
		tuple val(prefix) , path("*.amb") into index_not_terminated
	
		script:
		"""
		bwa index -p $prefix $genome
		"""

	}
	
	index_not_terminated.map {path ->
	filepath= path[1]
	filepath = filepath.toString()
	filepath = filepath.substring(0, filepath.length() -4)
	return [filepath]
	}
	.set{
		index_to_maps_sr
	}
	
}

else {

	Channel.value (index)
		.set { index_to_maps_sr }
}


read_clean_to_maps.combine(index_to_maps_sr).set {
	read_clean_to_maps
}


process mapping_with_srnamapper {


        publishDir "$baseDir/results/bam", mode: 'copy'

        input:
        tuple val(prefix), path(reads), val(index)  from read_clean_to_maps


        output:
        path ("*.bam") into bam_map_to_mmquant


        script:
        """
	srnaMapper -r $reads -g $index -o ${prefix}.sam 
	samtools view -b ${prefix}.sam > ${prefix}.bam
	"""

}

process mmquant {

	publishDir "$baseDir/results/mmquant", mode: "copy"


	input:
	path filebam from bam_map_to_mmquant.flatten().collect()
	path annotation_file from annotation_to_mmquant 

	output:
	path "*.tsv" into quantification_to_report

	script:
	"""
	mmquant -a $annotation_file -r $filebam  -o report.tsv
	"""

}  

config_multiqc = Channel.fromPath("$baseDir/multiqc.yaml")

process multiQC {


	publishDir "$baseDir/results/control", mode: "copy"

	input:
	path config from config_multiqc
	path '*' from fastqc_to_report.flatten().collect()
	path '*' from trim_to_report.flatten().collect()
	path '*' from quantification_to_report.flatten().collect()


	output:
	path 'multiqc_report.html' 

	script:
	"""
	multiqc --config $config .
	"""
}
	
	
