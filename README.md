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
`bam`|File|Input BAM file containing aligned reads
`bai`|File|Index file for BAM
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


./commands.txt found, printing out the content...
## Commands
  This section lists command(s) run by Optitype workflow
  
  * Running Optitype
  
  
  ```
                set -euo pipefail
  
                if [ -z "~{numReads}" ]; then
                totalLines=$(zcat ~{fastqR1} | wc -l)
                else totalLines=$((~{numReads}*4))
                fi
                python3 -c "from math import ceil; print (int(ceil(($totalLines/4.0)/~{numChunks})*4))"
 ```
 ```
                set -euo pipefail
                module load slicer/0.3.0
                slicer -i ~{fastqR} -l ~{chunkSize} --gzip 
 ```
 ```
                set -euo pipefail
                razers3 -i 95 -m 1 -dr 0 -o HLA.bam ~{hlaref} ~{fastq}
                samtools bam2fq HLA.bam > HLA.fastq       
 ```
 ```
                set -euo pipefail
                cat ~{sep=" " fastq} > hlareads.fastq
 ```
 ```
                module load optitype
                optitype -i ~{fastqR1} ~{fastqR2} --~{libtype} -v -o . --prefix ~{prefix}
 ``` ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_