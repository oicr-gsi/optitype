# optitype

OptiType does HLA genotyping producing 4-digit HLA genotyping predictions from NGS data and selects major/minor HLA Class I alleles. The workflow pre-filters input fastq reads by aligning to and HLA fasta reference based on library type (dna|rna) using RazerS3, as suggested in the tool documentation.

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
`fastqR1`|File|Fastq file for read 1
`fastqR2`|File|Fastq file for read 2
`outputFileNamePrefix`|String|Prefix for output files
`libtype`|String|the type of library, which will determine the hla reference to use.  dna|rna


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`numChunks`|Int|1|Number of chunks to split fastq file [1, no splitting]
`numReads`|Int?|None|Number of reads


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`countChunkSize.modules`|String|"python/3.7"|name and version of modules
`countChunkSize.jobMemory`|Int|16|Memory allocated for this job
`countChunkSize.timeout`|Int|48|Hours before task timeout
`slicerR1.modules`|String|"slicer/0.3.0"|Required environment modules
`slicerR1.jobMemory`|Int|16|Memory allocated for this job
`slicerR1.timeout`|Int|48|Hours before task timeout
`slicerR2.modules`|String|"slicer/0.3.0"|Required environment modules
`slicerR2.jobMemory`|Int|16|Memory allocated for this job
`slicerR2.timeout`|Int|48|Hours before task timeout
`HLAReadsR1.modules`|String|"optitype/1.3.1"|Required environment modules
`HLAReadsR1.jobMemory`|Int|16|Memory allocated for this job
`HLAReadsR1.timeout`|Int|48|Hours before task timeout
`HLAReadsR2.modules`|String|"optitype/1.3.1"|Required environment modules
`HLAReadsR2.jobMemory`|Int|16|Memory allocated for this job
`HLAReadsR2.timeout`|Int|48|Hours before task timeout
`concatR1.jobMemory`|Int|16|Memory allocated for this job
`concatR1.timeout`|Int|48|Hours before task timeout
`concatR2.jobMemory`|Int|16|Memory allocated for this job
`concatR2.timeout`|Int|48|Hours before task timeout
`run_optitype.modules`|String|"optitype/1.3.1"|Required environment modules
`run_optitype.jobMemory`|Int|16|Memory allocated for this job
`run_optitype.timeout`|Int|48|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`optitypeResults`|File|{'description': 'Results of optitype', 'vidarr_label': 'optitypeResults'}
`optitypePlot`|File|{'description': 'Plots of optitype', 'vidarr_label': 'optitypePlot'}


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
```

 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
