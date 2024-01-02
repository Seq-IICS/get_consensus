# GET_CONSENSUS (AMPLICONS)

[TOC]



Run pipeline via [nextflow](https://www.nextflow.io/). Notes on building and pushing Docker image are [here](https://hub.docker.com/layers/vici0uz/iics/latest/images/sha256-1c11f1fa0b9a1d3f5cb9a7f7e9e6eb97dfc23b00ea3b3fd3e44ff6bd00c8b068?context=repo)

## REQUIREMENTS
- nextflow
- docker


## SETUP
data directory containing
- barcoded | demultiplexed fastq files
- your_samplesheet.tsv
- your_reference.fasta
- a config.json file with the following structure

```
{
"sample": "your_samplesheet.tsv",
    "reference": "your_reference.fasta"
}
```
Execute on linux or wsl terminal
```
nextflow run Seq-IICS/get_consensus --path {myDataDir} -with-docker vici0uz/iics:latest
```
## Tools used
- [Nanoplot](https://github.com/wdecoster/NanoPlot)
- [Minimap2](https://github.com/lh3/minimap2)
- [Samtools](https://github.com/samtools/)
- [iVar](https://github.com/andersen-lab/ivar)