# bam2fastq
Nextflow pipeline for conversion of BAM to FASTQ

Usage:

```
nextflow run ameynert/bam2fastq --input 'dir/*.bam' --outdir /path/to/outdir --tmpdir /path/to/tmpdir

Mandatory arguments:
  --input                Path to input data (must be surrounded with quotes)
  --outdir               The output directory where the results will be saved
  --tmpdir               Directory for large temporary files

Other options:
  --cram                 Flag that input is CRAM files
  --reference            If the input is CRAM files, path to the reference FASTA file
  -name                  Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
```
