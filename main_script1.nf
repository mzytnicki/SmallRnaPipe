/*
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

	publishDir "$baseDir/results/fastqc", mode: 'copy'

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
	prinseq-lite -fastq $reads_trim -lc_method dust -lc_threshold 7
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
		echo $genome
		bwa index -p $prefix $genome
		"""

	}
	
	index_not_terminated.map {path ->
	filepath= path[1]
	filepath = filepath.toString()
	filepath = filepath.substring(0, filepath.length() -4)
	return [path [0], filepath]
	}
	.set{
		index_to_maps_sr
	}
	
}

else {

	Channel.value (index)
		.map { path ->
		filename = path
		path = path
		return [filename, path]
}
		.set { index_to_maps_sr }
}


process mapping_with_srnamapper {


        publishDir "$baseDir/results/sam", mode: 'copy'

        input:
        tuple val(nothing), val(index)  from index_to_maps_sr
        tuple val(prefix), path(reads) from read_clean_to_maps


        output:
        tuple val(prefix), path( "*.sam") into sam_map_to_mmquant


        script:
        """
	
	srnaMapper -r $reads -g $index -o ${prefix}.sam
	"""

}


process mmquant {

	publishDir "$baseDir/results/mmquant", mode: "copy"


	input:
	tuple val(prefix), path(filesam) from sam_map_to_mmquant
	path annotation_file from annotation_to_mmquant

	output:
	path "*.tsv" into quantification_to_report

	script:
	"""
	mmquant -a $annotation_file -r $filesam -o ${prefix}.tsv
	"""

}  

process multiQC {


	publishDir "$baseDir/results/control", mode: "copy"

	input:
	path fastqc from fastqc_to_report
	path trimming from trim_to_report
	path quantify from quantification_to_report


	output:
	path 'multiqc_report.html' 

	script:
	"""
	multiqc $fastqc $trimming $quantify
	"""
}
	
	
