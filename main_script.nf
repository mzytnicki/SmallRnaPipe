#!/bin/bash

/
* Parse and check parameters
*/

reads = params.containsKey('reads') ? params.reads : ''
index = params.containsKey('index') ? params.index : ''
genome = params.containsKey('genome') ? params.genome : ''
annotation = params.containsKey('annotation') ? params.annotation : ''
mature = params.containsKey('mature') ? params.mature : ''
hairpin = params.containsKey('hairpin') ? params.hairpin : ''
species = params.containsKey('species') ? params.species : ''
o_species = params.containsKey('o_species') ? params.o_species : ''

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

Channel.fromPath(reads)
	.map { path ->
	filename= path.getSimpleName()
	return [filename,path]
}
	.set{ 
	reads_ch_to_decompressed
}

Channel.fromPath(genome)
        .map { path ->
        filename= path.getSimpleName()
        return [filename,path]
}
        .into{
        genome_to_index
        genome_to_maps
	genome_to_decompress
}

Channel.fromPath(annotation)
	.set{ annotation_to_mmquant
	}
Channel.fromPath(mature)
	.set{ mature_to_decompress
	}
Channel.fromPath(hairpin)
	.set{ hairpin_to_decompress
	}	

Channel.from(species)
	.into{ 
	species_to_hairpin_extraction
	species_to_mature_extraction
	species_to_mirdeep2
}	
Channel.from(o_species)
	.set{ other_species_to_mature_extraction
}

/*
* Process to decompress FASTQ file 
*/


process decompress_reads {

        input : 
        tuple val(prefix), path(reads) from reads_ch_to_decompressed

        output:
        path '*.fastq' into fastq_decompressed_to_control
        tuple val(prefix), path("*fastq") into fastq_decompressed_to_trimming

        script:
	"""
	if [[ "\$(tail -c 8 <<<$reads)"  == .tar.gz ]]; then
        	tar -xf $reads -O > ${prefix}.fastq
        elif [[ "\$(tail -c 4 <<<$reads)" == .gz ]]; then
        	gzip -c -d $reads > ${prefix}.fastq
        fi
	"""
}

/*
* Process to decompress file for the different steps of mirdeep2
*/

process decompress_to_mirdeep2 {

        input :
        path mature from mature_to_decompress
	path hairpin from hairpin_to_decompress
	tuple val(prefix), path(genome) from genome_to_decompress

        output:
        path 'mature.fa' into mature_to_extraction
	path "mature.fa" into mature_to_extraction_other_species
	path 'hairpin.fa' into hairpin_to_extraction
	tuple val(prefix), path('*_edited.fa') into genome_to_bowtie_index
	path '*_edited.fa' into genome_to_mirdeep2
	
        script:
        """
	MATURE="$mature"
	HAIRPIN="$hairpin"
	if [ \${MATURE: -3} == ".gz" ]; then
        	gunzip -f \$MATURE
        	MATURE=\${MATURE%%.gz}
	fi
	if [ \${HAIRPIN: -3} == ".gz" ]; then
        	gunzip -f \$HAIRPIN
        	HAIRPIN=\${HAIRPIN%%.gz}
	fi
	
	sed '/^[^>]/s/[^ATGCatgc]/N/g' $genome > ${prefix}_edited.fa

	sed -i 's, ,_,g' \$HAIRPIN
	sed -i 's, ,_,g' \$MATURE

        """
}

/*
* Process to extract the sequences of the hairpin RNAs of your species in the hairpin ref
*/

process extract_hairpin_target {
	
	publishDir "$baseDir/results/hairpin_target", mode: 'copy'

	input:
	path hairpin from hairpin_to_extraction
	val species from species_to_hairpin_extraction

	output:

	path 'hairpin_ref.fa' into hairpin_to_mirdeep2 
	
	script:
	"""
	extract_miRNAs.pl $hairpin $species > hairpin_ref.fa
	"""

}

/*
* Process to extract the sequences of the mature RNAs of your species in the mature ref
*/


process extract_mature_target {

        publishDir "$baseDir/results/mature_target", mode: 'copy'

        input:
        path mature from mature_to_extraction
        val species from species_to_mature_extraction

        output:

        path 'mature_ref.fa' into mature_to_mirdeep2

        script:
        """
        extract_miRNAs.pl $mature $species > mature_ref.fa
        """

}

/*
* Process to extract the sequences of the mature RNAs of a other species in the mature ref
*/


process extract_other_mature_target {

        publishDir "$baseDir/results/other_mature", mode: 'copy'

        input:
        path mature from mature_to_extraction_other_species
        val other_species from other_species_to_mature_extraction

        output:

        path 'mature_other.fa' into mature_other_to_mirdeep2

        script:
        """
        extract_miRNAs.pl $mature $other_species > mature_other.fa
        """


}

/*
* Process to index the genome with bowtie. It is for mirdeep2 process 
*/

process index_with_bowtie {
	
	publishDir "$baseDir/results/bowtie_index", mode: 'copy'

	input:
	tuple val(prefix), path(genome) from genome_to_bowtie_index

	output:
	path '*.ebwt' into genome_for_mapper
	script:
	"""
	bowtie-build $genome ${prefix}
	"""

}

/*
* Process to control quality of FASTQ before the cleaning
*/

process control_quality_with_fastqc {
	
	publishDir "$baseDir/results/fastqc" , mode: 'copy'
	
	input :
	path reads from fastq_decompressed_to_control

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
        tuple val(prefix), path(reads) from fastq_decompressed_to_trimming

        output :
	path "*trimming_report.txt" into trim_to_report
	tuple val(prefix),  path("*_trimmed.fq") into trim_read_to_prinseq

	script:
	"""
	trim_galore $reads --basename $prefix
	"""

}

/*
* Process to clean the reads when it has a low complexity
*/
	
process cleaning_reads_with_low_complexity {

	publishDir "$baseDir/results/prinseq" , mode: 'copy'

	input:

	tuple val(prefix), path(reads_trim) from trim_read_to_prinseq

	output:
	tuple val(prefix), path("*good*.fastq") into read_clean_to_maps
	tuple val(prefix), path("*good*.fastq") into read_clean_to_mapper_mirdeep2
	tuple val(prefix), path("*good*.fastq") into control_quality_of_goods_reads
	script:
	"""
	prinseq-lite.pl -fastq $reads_trim -lc_method dust -lc_threshold 7
	"""
}

/*
* Process to control the quality of FASTQ after the cleaning 
*/

process control_quality_good_reads_with_fastqc {
	
	publishDir "$baseDir/results/fastqc_good_reads" , mode: 'copy'

        input :
        tuple val(prefix), path(reads)  from control_quality_of_goods_reads

        output:
        path '*_fastqc.zip' into fastqc_good_reads_to_report

        script:
        """
        fastqc $reads
        """
}

/*
* If you do not have a bwa index, the process will createa bwa index  
*/

if (!index){

	process bwa_index_to_create {

		publishDir "$baseDir/results/index_with_bwa" , mode: 'copy'


		input:
		tuple val(prefix), path(genome) from genome_to_index

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

/*
* Process to map the reads with a specific function on mirdeep2 package
*/

process mapper_for_mirdeep2 {

	input:
	tuple val(prefix), path(reads) from read_clean_to_mapper_mirdeep2
	path genome from genome_for_mapper.collect()

	output:
	path '*_collapsed.fa' into reads_collapsed_to_mirdeep2
	path '*_reads_vs_refdb.arf' into reads_vs_refdb_to_mirdeep2

	script:

	index_base = genome.toString().tokenize(' ')[0].tokenize('.')[0]
	"""
	mapper.pl $reads -e -h -i -j -m -p $index_base -s ${prefix}_collapsed.fa -t ${prefix}_reads_vs_refdb.arf -o 4
	"""
}

reads_collapsed_to_mirdeep2.combine( genome_to_mirdeep2)
	.combine(species_to_mirdeep2)
	.combine(mature_other_to_mirdeep2)
	.combine(mature_to_mirdeep2)
	.combine(hairpin_to_mirdeep2)
	.set{files_to_mirdeep2}

/*
* Process to discover new mRNAs and to annotate the reads with mirdeep2
*/

process mirdeep2 {
	publishDir "$baseDir/results/mirdeep2", mode: 'copy'

	input:
	tuple path(collapsed), path(genome), val(species), path(other_mature), path(mature), path(hairpin) from files_to_mirdeep2
	path vs_refdb from reads_vs_refdb_to_mirdeep2
	output:
	file 'result*.{bed,csv,html}'

	script:
	"""
	perl -ane 's/[ybkmrsw]/N/ig;print;' $hairpin > hairpin_ok.fa
	sed 's/ .*//' $genome | awk '\$1 ~ /^>/ {gsub(/_/,"",\$1); print; next} {print}' > genome_nowhitespace.fa
	miRDeep2.pl $collapsed genome_nowhitespace.fa $vs_refdb $mature $other_mature hairpin_ok.fa -t $species 
	"""
}

/*
* Process to map the reads with srnaMapper
*/

process mapping_with_srnamapper {


        publishDir "$baseDir/results/bam", mode: 'copy'

        input:
        tuple val(prefix), path(reads), val(index)  from read_clean_to_maps


        output:
        path ("*.bam") into bam_to_mmquant
        tuple val (prefix), path ("*.bam") into bam_to_flagstats
        tuple val (prefix), path ("*.bam") into bam_to_idxstats
        tuple val (prefix), path ("*.bam") into bam_to_stats

        script:
        """
	srnaMapper -r $reads -g $index -o ${prefix}.sam 
	samtools sort ${prefix}.sam > ${prefix}.bam
	"""

}

/*
* Process to control the alignements with samtools flagstat
*/


process control_alignments {
	
	publishDir "$baseDir/results/flagstats",mode: "copy"

	input :
	tuple val(prefix), path (bam) from bam_to_flagstats

	output :
	path ('*.flagstat') into control_alignements_to_report

	script:
	"""
	samtools flagstat $bam > ${prefix}.flagstat
	"""
}

/*
* Process to control the metrics data with samtools stats
*/

process control_metrics {

        publishDir "$baseDir/results/stats",mode: "copy"

	input :
	tuple val(prefix), path(bam) from bam_to_stats

	output :
	path ('*.stats') into control_metrics_to_report

	script:
	"""
	samtools stats $bam > ${prefix}.stats
	"""
}

/*
* Process to control the contifs with samtools idxstats
*/

process control_contigs {

        publishDir "$baseDir/results/idxstats",mode: "copy"

	input :
	tuple val(prefix), path(bam) from bam_to_idxstats

	output :
	path ('*.idxstats') into control_contigs_to_report

	script:
	"""
	samtools index $bam
	samtools idxstats $bam > ${prefix}.idxstats
	"""
}

/*
* Process to quantify each annotation on your genome
*/

process mmquant {

	publishDir "$baseDir/results/mmquant", mode: "copy"


	input:
	path filebam from bam_to_mmquant.flatten().collect()
	path annotation_file from annotation_to_mmquant 

	output:
	path "*report.summary" into quantification_to_report
	path "final_tab2.tsv" into quantification_table_to_report 
	script:

	"""
	mmquant -a $annotation_file -r $filebam  -F -O tab.tsv > quantification_report.summary
	sed -e 's/                       /;/g' -e 's/                     //g' -e 's/           //g' -e 's/:          /:/g' -e '2s/          /;/g' -e 's/          /;/g' -e 's/    //g' -e 's/:  /:/g' -e 's/(......)/;/g' -e 's/:  /:/g' -e 's/:/;/g' -e 's/#//g' -e 's/;/\t/g'  tab.tsv > final_tab.tsv
	datamash transpose < final_tab.tsv > final_tab2.tsv
	"""

}  

config_multiqc = Channel.fromPath("$baseDir/multiqc.yaml")

/*
* Process for displaying a multiqc report in an html
*/

process multiQC {


	publishDir "$baseDir/results/control", mode: "copy"

	input:
	path config from config_multiqc
	path '*' from fastqc_to_report.flatten().collect()
	path '*' from trim_to_report.flatten().collect()
	path '*' from quantification_to_report.flatten().collect()
	path '*' from fastqc_good_reads_to_report.flatten().collect()
	path '*' from quantification_table_to_report.flatten().collect()
        path '*' from control_alignements_to_report.flatten().collect()
        path '*' from control_metrics_to_report.flatten().collect()
        path '*' from control_contigs_to_report.flatten().collect()

	output:
	path 'multiqc_report.html' 

	script:
	"""
	multiqc --config $config .
	"""
}
	
	
