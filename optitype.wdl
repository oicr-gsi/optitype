version 1.0

workflow optitype {
	input {
		Array[File] fastqR1
		Array[File] fastqR2
		String outputFileNamePrefix
		Int numChunks = 1
		Int? numReads
		String libtype
	}

	parameter_meta {
		fastqR1: "Array or single of Fastq files for read 1"
		fastqR2: "Array or single of Fastq files for read 2"
		outputFileNamePrefix: "Prefix for output files"
		numChunks: "Number of chunks to split fastq file [1, no splitting]"
		numReads: "Number of reads"
		libtype: "the type of library, which will determine the hla reference to use.  dna|rna"
	}
    
	
	Map[String,File] ref_fasta = {
		"dna":"/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/optitype-1.3.1/ref/hla_reference_dna.fasta",
		"rna":"/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/optitype-1.3.1/ref/hla_reference_rna.fasta"
	}
	
	### concatenating raw files
	call concatRawReads {
  		input:
    		fastqR1 = fastqR1,
    		fastqR2 = fastqR2
	}
	
	### if this is greater than 1, then the fastq files will be split to chunks and
	if (numChunks > 1) {
		call countChunkSize {
			input:
			fastqR1 = concatRawReads.raw_concatenated_fastqR1,
			numChunks = numChunks,
			numReads = numReads
		}
		call slicer as slicerR1{
			input:
			fastqR = concatRawReads.raw_concatenated_fastqR1,
			chunkSize = countChunkSize.chunkSize
		}
		call slicer as slicerR2 {
			input:
			fastqR = concatRawReads.raw_concatenated_fastqR2,
			chunkSize = countChunkSize.chunkSize
		}
	}

	Array[File] fastq1 = select_first([slicerR1.chunkFastq, [concatRawReads.raw_concatenated_fastqR1]])
	Array[File] fastq2 = select_first([slicerR2.chunkFastq, [concatRawReads.raw_concatenated_fastqR2]])
    
	scatter(fq1 in fastq1){
		call HLAReads as HLAReadsR1{
		    input:
			fastq = fq1,
			hlaref = ref_fasta[libtype]
		}
	}	
	scatter(fq2 in fastq2){
		call HLAReads as HLAReadsR2{
		    input:
			fastq = fq2,
			hlaref = ref_fasta[libtype]
		}
	}	
	
	Array[File] hlafastq1 = HLAReadsR1.hlafastq
	Array[File] hlafastq2 = HLAReadsR2.hlafastq
	
	call concatReads as concatR1{
		input:
		fastq=hlafastq1
	}
	call concatReads as concatR2{
		input:
		fastq=hlafastq2
	}
	
	call run_optitype{
		input:
		fastqR1 = concatR1.concatenated_fastq,
		fastqR2 = concatR2.concatenated_fastq,
		prefix = outputFileNamePrefix,
		libtype = libtype
	}
	output {
		File optitypeResults = run_optitype.results
		File optitypePlot = run_optitype.plot
    }
	
	meta {
		author: "Lawrence Heisler, Monica L. Rojas-Pena"
		email: "lheisler@oicr.on.ca, mrojaspena@oicr.on.ca"
		description: "OptiType does HLA genotyping producing 4-digit HLA genotyping predictions from NGS data and selects major/minor HLA Class I alleles. The workflow pre-filters input fastq reads by aligning to and HLA fasta reference based on library type (dna|rna) using RazerS3, as suggested in the tool documentation."
		dependencies: [
			{
				name: "optiType/1.3.1",
				url: "https://github.com/FRED-2/OptiType"
			},
			{
				name: "razers3/3.5.8",
				url: "http://packages.seqan.de/razers3/razers3-3.5.8-Linux-x86_64.tar.xz"
			},
			{
				name: "slicer/0.3.0",
				url: "https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz"
			},
			{ 
				name: "gsi software modules : optitype/1.3.1 slicer/0.3.0",
			    url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
			},
			{
				name: "gsi software module dependencies : singularity/3.9.4 razers/3.5.8 samtools/1.16.1",
				url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
			}
		]
		output_meta: {
			optitypeResults: {
				description: "Results of optitype",
				vidarr_label: "optitypeResults"
			},
			optitypePlot: {
				description: "Plots of optitype",
				vidarr_label: "optitypePlot"
			}
		}
	}
}

# Task to concatenate raw fastq files before chunking
task concatRawReads {
  input {
    Array[File] fastqR1
    Array[File] fastqR2
    Int jobMemory = 16
    Int timeout = 48
  }

  command <<<
    set -euo pipefail
    zcat ~{sep=" " fastqR1} | gzip > raw_concat_R1.fastq.gz
    zcat ~{sep=" " fastqR2} | gzip > raw_concat_R2.fastq.gz
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    timeout: "~{timeout}"
  }

  output {
    File raw_concatenated_fastqR1 = "raw_concat_R1.fastq.gz"
    File raw_concatenated_fastqR2 = "raw_concat_R2.fastq.gz"
  }
}

task countChunkSize{
	input {
		File fastqR1
		Int numChunks
		Int ?numReads
		String modules = "python/3.7"
		Int jobMemory = 16
		Int timeout = 48
	}
	parameter_meta {
		fastqR1: "input fatsq file"
		numChunks: " number of chunks"
		numReads: "number of reads"
		modules: "name and version of modules"
		jobMemory: "Memory allocated for this job"
		timeout: "Hours before task timeout"
	}
	command <<<
		set -euo pipefail
		if [ -z "~{numReads}" ]; then
			totalLines=$(zcat ~{fastqR1} | wc -l)
		else totalLines=$((~{numReads}*4))
		fi
		python3 -c "from math import ceil; print (int(ceil(($totalLines/4.0)/~{numChunks})*4))"
	>>>
	runtime {
		modules: "~{modules}"
		memory: "~{jobMemory} GB"
		timeout: "~{timeout}"
	}
	output {
		String chunkSize =read_string(stdout())
	}
}
task slicer {
	input {
		File fastqR         
		String chunkSize
		String modules = "slicer/0.3.0"
		Int jobMemory = 16
		Int timeout = 48
	}
	parameter_meta {
		fastqR: "Fastq file"
		chunkSize: "Number of lines per chunk"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job"
		timeout: "Hours before task timeout"
	}
	command <<<
		set -euo pipefail
		module load slicer/0.3.0
		slicer -i ~{fastqR} -l ~{chunkSize} --gzip 
	>>>
	runtime {
		modules: "~{modules}"
		memory: "~{jobMemory} GB"
		timeout: "~{timeout}"
	}
	output {
		Array[File] chunkFastq = glob("*.fastq.gz")
	}
	meta {
		output_meta: {
			chunkFastq: "output fastq chunks"
		}
	}
}
task HLAReads{
	input {
		File fastq
		File hlaref
		String modules = "optitype/1.3.1 hla-reference/1.0.0"
		Int jobMemory = 16
		Int timeout = 48
	}
	parameter_meta {
		fastq: "Fastq file"
		hlaref: "hla reference file"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job"
		timeout: "Hours before task timeout"
	}
	command <<<
		set -euo pipefail
		razers3 -i 95 -m 1 -dr 0 -o HLA.bam ~{hlaref} ~{fastq}
		samtools bam2fq HLA.bam > HLA.fastq	  
	>>>
	runtime {
		modules: "~{modules}"
		memory: "~{jobMemory} GB"
		timeout: "~{timeout}"
	}
	output {
		File hlafastq = "HLA.fastq"
	} 
}
task concatReads{
	input{
		Array[File] fastq
		Int jobMemory = 16
		Int timeout = 48
	}
	parameter_meta {
		fastq: "Array of input fastq files"
		jobMemory: "Memory allocated for this job"
		timeout: "Hours before task timeout"
	}
	command <<<
		set -euo pipefail
		cat ~{sep=" " fastq} > hlareads.fastq
	>>>
	runtime {

		memory: "~{jobMemory} GB"
		timeout: "~{timeout}"
	}
	output {
		File concatenated_fastq = "hlareads.fastq"
	}
}
task run_optitype{
	input{
		File fastqR1
		File fastqR2
		String prefix
		String libtype
		String modules = "optitype/1.3.1"
		Int jobMemory = 16
		Int timeout = 48
	}
	parameter_meta {
		fastqR1: "Fastq file 1ead1"
		fastqR2: "Fastq file read2"
		prefix: "Prefix for output files"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job"
		timeout: "Hours before task timeout"
	}
	command <<<
	module load optitype
	optitype -i ~{fastqR1} ~{fastqR2} --~{libtype} -v -o . --prefix ~{prefix}
	>>>
	runtime {
		modules: "~{modules}"
		memory: "~{jobMemory} GB"
		timeout: "~{timeout}"
	}
	output {
		File results = "~{prefix}_result.tsv"
		File plot = "~{prefix}_coverage_plot.pdf"
 	}
}