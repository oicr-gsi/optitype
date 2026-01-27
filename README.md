# optitype

OptiType performs 4-digit HLA genotyping predictions from NGS data, selecting major/minor HLA Class I alleles. The workflow pre-filters input fastq reads by aligning to an HLA fasta reference using RazerS3, per tool documentation.

## Overview

## Dependencies

* [optiType 1.3.1](https://github.com/FRED-2/OptiType)
* [razers3 3.5.8](http://packages.seqan.de/razers3/razers3-3.5.8-Linux-x86_64.tar.xz)
* [slicer 0.3.0](https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz)
* [gsi software modules : optitype 1.3.1 slicer 0.3.0](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi software module dependencies : singularity 3.9.4 razers 3.5.8 samtools 1.16.1](https://gitlab.oicr.on.ca/ResearchIT/modulator)


## Usage

### Cromwell
```
java -jar cromwell.jar run optitype.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`bam`|Array[File]|One or more BAM files for a single sample
`bai`|Array[File]|BAM index files (same order as bam)
`outputFileNamePrefix`|String|Prefix for output files
`libtype`|String|Library type determining HLA reference to use (dna|rna)


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extract_chr6_HLA_region.modules`|String|"samtools/1.9"|Modules to load
`extract_chr6_HLA_region.jobMemory`|Int|16|Memory allocated (GB)
`extract_chr6_HLA_region.timeout`|Int|48|Timeout in hours
`HLAReads.modules`|String|"optitype/1.3.1 hla-reference/1.0.0"|Required environment modules
`HLAReads.threads`|Int|8|Number of threads to use for razers3
`HLAReads.jobMemory`|Int|24|Memory allocated for this job
`HLAReads.timeout`|Int|72|Hours before task timeout
`run_optitype.modules`|String|"optitype/1.3.1"|Required environment modules
`run_optitype.jobMemory`|Int|16|Memory allocated for this job
`run_optitype.timeout`|Int|48|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`optitypeResults`|File|Results of optitype|vidarr_label: optitypeResults
`optitypePlot`|File|Plots of optitype|vidarr_label: optitypePlot


## Commands
  This section lists command(s) run by Optitype workflow
  
  * Running Optitype
 
  ```
   set -euo pipefail
   module load ~{modules}
 
   # Convert input Array[File] to bash arrays
   bams=(~{sep=' ' bam})
   bais=(~{sep=' ' bai})
 
   # Initialize array to store per-BAM filtered BAMs
   chr6_bams=()
 
   # Loop over each BAM
   for i in $(seq 0 $((${#bams[@]}-1))); do
       # Extract chr6 main region
      samtools view -b "${bams[$i]}" chr6:29677984-33485635 > chr6_region_$i.bam
 
       # Extract chr6 alt contigs
       samtools view -b "${bams[$i]}" \
           chr6_GL000250v2_alt chr6_GL000251v2_alt chr6_GL000252v2_alt \
           chr6_GL000253v2_alt chr6_GL000254v2_alt chr6_GL000255v2_alt \
           chr6_GL000256v2_alt chr6_GL383533v1_alt chr6_KB021644v2_alt \
           chr6_KI270758v1_alt chr6_KI270797v1_alt chr6_KI270798v1_alt \
           chr6_KI270799v1_alt chr6_KI270800v1_alt chr6_KI270801v1_alt \
           chr6_KI270802v1_alt chr6_KQ090017v1_alt > chr6_alt_$i.bam
 
      # Conditionally extract unmapped reads if libtype is 'rna'
 
         # Note:
         # Unmapped reads are included only for RNA-seq inputs to mitigate the
         # systematic loss of HLA reads during initial genome alignment due to extreme polymorphism and multi-mapping.
         # OptiType performs HLA-aware realignment, making these reads informative
         # for RNA but unnecessary and potentially noisy for DNA inputs.
 
       if [ "~{libtype}" == "rna" ]; then
           samtools view -h -b -f 4 "${bams[$i]}" > unmapped_$i.bam
           samtools merge chr6_filtered_$i.bam chr6_region_$i.bam chr6_alt_$i.bam unmapped_$i.bam
       else
           samtools merge chr6_filtered_$i.bam chr6_region_$i.bam chr6_alt_$i.bam
       fi
 
       # Add to array
       chr6_bams+=("chr6_filtered_$i.bam")
   done
 
   # Merge all filtered BAMs if more than one
   if [ ${#chr6_bams[@]} -gt 1 ]; then
       samtools merge -@ 8 chr6_merged.bam "${chr6_bams[@]}"
      samtools index -@ 8 chr6_merged.bam
   else
       cp "${chr6_bams[0]}" chr6_merged.bam
      samtools index -@ 8 chr6_merged.bam
   fi
 
   # Sort by read name and convert to FASTQ
   samtools sort -n -o chr6_merged.sorted.bam chr6_merged.bam
 
   # Convert to paired FASTQ files
   samtools fastq -1 chr6_R1.fastq -2 chr6_R2.fastq -0 /dev/null -s /dev/null -n chr6_merged.sorted.bam
 
   # gzip fastq files
   gzip -c chr6_R1.fastq > chr6_R1.fastq.gz
   gzip -c chr6_R2.fastq > chr6_R2.fastq.gz
 ```
 ```
  	set -euo pipefail
     echo "Starting razers3 at $(date)"
     razers3 -i 95 -tc "~{threads}" -m 1 -dr 0 -o HLA1_R1.bam ~{hlaref} ~{fastqR1}
     razers3 -i 95 -tc "~{threads}" -m 1 -dr 0 -o HLA2_R2.bam ~{hlaref} ~{fastqR2}
     echo "razers3 finished at $(date)"
     samtools bam2fq HLA1_R1.bam > HLA_R1.fastq
     samtools bam2fq HLA2_R2.bam > HLA_R2.fastq 
 ```
 ```
     module load ~{modules}
     optitype -i ~{hlafastq_R1} ~{hlafastq_R2} --~{libtype} -v -o . --prefix ~{prefix}	  
 ``` ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
