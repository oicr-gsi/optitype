version 1.0

workflow optitype {
  input {
    File bam
    File bai
    String outputFileNamePrefix
    String libtype
  }

  parameter_meta {
    bam: "Input BAM file containing aligned reads"
    bai: "Index file for BAM"
    outputFileNamePrefix: "Prefix for output files"
    libtype: "Library type determining HLA reference to use (dna|rna)"
  }

  Map[String, File] ref_fasta = {
    "dna": "/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/optitype-1.3.1/ref/hla_reference_dna.fasta",
    "rna": "/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/optitype-1.3.1/ref/hla_reference_rna.fasta"
  }

  call extract_chr6_HLA_region {
    input:
      bam = bam,
      bai = bai,
      libtype = libtype
  }

  call HLAReads {
    input:
      fastqR1 = extract_chr6_HLA_region.fastqR1,
      fastqR2 = extract_chr6_HLA_region.fastqR2,
      hlaref = ref_fasta[libtype]
  }

  call run_optitype {
    input:
      hlafastq_R1 = HLAReads.hlafastq_R1,
      hlafastq_R2 = HLAReads.hlafastq_R2,
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
    description: "OptiType performs 4-digit HLA genotyping predictions from NGS data, selecting major/minor HLA Class I alleles. The workflow pre-filters input fastq reads by aligning to an HLA fasta reference using RazerS3, per tool documentation."
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

task extract_chr6_HLA_region {
  input {
    File bam
    File bai
    String libtype
    String modules = "samtools/1.9"
    Int jobMemory = 16
    Int timeout = 48
  }

  parameter_meta {
    bam: "Input BAM file"
    bai: "BAM index file"
    modules: "Modules to load"
    jobMemory: "Memory allocated (GB)"
    timeout: "Timeout in hours"
  }

  command <<< 
  set -euo pipefail
  module load ~{modules}

  # Extract all contigs starting with chr6 HLA region
  samtools view -b ~{bam} chr6:29677984-33485635 > chr6_region.bam

  # Extract chr6 alt reads
  samtools view -b ~{bam} \
    chr6_GL000250v2_alt chr6_GL000251v2_alt chr6_GL000252v2_alt \
    chr6_GL000253v2_alt chr6_GL000254v2_alt chr6_GL000255v2_alt \
    chr6_GL000256v2_alt chr6_GL383533v1_alt chr6_KB021644v2_alt \
    chr6_KI270758v1_alt chr6_KI270797v1_alt chr6_KI270798v1_alt \
    chr6_KI270799v1_alt chr6_KI270800v1_alt chr6_KI270801v1_alt \
    chr6_KI270802v1_alt chr6_KQ090017v1_alt > chr6_alt.bam

  # Conditionally extract unmapped reads if libtype is 'rna'
  if [ "~{libtype}" == "rna" ]; then
    samtools view -h -b -f 4 ~{bam} > unmapped.bam
    samtools merge chr6_filter.bam chr6_region.bam chr6_alt.bam unmapped.bam
  else
    samtools merge chr6_filter.bam chr6_region.bam chr6_alt.bam
  fi

  # Sort by read name
  samtools sort -n -o chr6_filtered.sorted.bam chr6_filter.bam

  # Convert to paired FASTQ files
  samtools fastq -1 chr6_R1.fastq -2 chr6_R2.fastq -0 /dev/null -s /dev/null -n chr6_filtered.sorted.bam

  # gzip fastq files
  gzip -c chr6_R1.fastq > chr6_R1.fastq.gz
  gzip -c chr6_R2.fastq > chr6_R2.fastq.gz
>>>

  runtime {
    modules: "~{modules}"
    memory: "~{jobMemory} GB"
    timeout: "~{timeout}"
  }

  output {
    File fastqR1 = "chr6_R1.fastq.gz"
    File fastqR2 = "chr6_R2.fastq.gz"
  }
}

task HLAReads {
  input {
    File fastqR1
    File fastqR2
    File hlaref
    String modules = "optitype/1.3.1 hla-reference/1.0.0"
    Int threads = 8
    Int jobMemory = 24
    Int timeout = 72
  }

  parameter_meta {
    fastqR1: "Fastq file read 1"
    fastqR2: "Fastq file read 2"
    hlaref: "hla reference file"
    modules: "Required environment modules"
    threads: "Number of threads to use for razers3"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail
    echo "Starting razers3 at $(date)"
    razers3 -i 95 -tc "~{threads}" -m 1 -dr 0 -o HLA1_R1.bam ~{hlaref} ~{fastqR1}
    razers3 -i 95 -tc "~{threads}" -m 1 -dr 0 -o HLA2_R2.bam ~{hlaref} ~{fastqR2}
    echo "razers3 finished at $(date)"
    samtools bam2fq HLA1_R1.bam > HLA_R1.fastq
    samtools bam2fq HLA2_R2.bam > HLA_R2.fastq
  >>>

  runtime {
    modules: "~{modules}"
    memory: "~{jobMemory} GB"
    cpu: "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File hlafastq_R1 = "HLA_R1.fastq"
    File hlafastq_R2 = "HLA_R2.fastq"
  }
}

task run_optitype {
  input {
    File hlafastq_R1
    File hlafastq_R2
    String prefix
    String libtype
    String modules = "optitype/1.3.1"
    Int jobMemory = 16
    Int timeout = 48
  }

  parameter_meta {
    hlafastq_R1: "HLA Fastq file read1"
    hlafastq_R2: "HLA Fastq file read2"
    prefix: "Prefix for output files"
    libtype: "Library type (dna or rna)"
    modules: "Required environment modules"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }

  command <<<
    module load ~{modules}
    optitype -i ~{hlafastq_R1} ~{hlafastq_R2} --~{libtype} -v -o . --prefix ~{prefix}
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
