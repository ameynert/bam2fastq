#!/usr/bin/env nextflow
/*
========================================================================================
                         bam2fastq
========================================================================================
 https://github.com/ameynert/bam2fastq
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run ameynert/bam2fastq --input 'dir/*.bam' --outdir /path/to/outdir

    Mandatory arguments:
      --input                Path to input data (must be surrounded with quotes)
      --outdir               The output directory where the results will be saved
      --tmpdir               Directory for large temporary files

    Other options:
      --cram                 Flag that input is CRAM files
      --reference            If the input is CRAM files, path to the reference FASTA file
      -name                  Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (!params.input) {
    exit 1, "Input files not specified"
}

if (!params.outdir) {
    exit 1, "Output directory not specified"
}

if (!params.tmpdir) {
    exit 1, "Temporary directory not specified"
}

if (params.cram) {
  if (!params.reference) {
    exit 1, "Reference must be specified for CRAM files"
  }
  reference_ch = Channel.value(file(params.reference))
}

/*
 * Create a channel for input alignment files
 */
Channel
  .fromFilePairs( params.input, size: 1 )
  .ifEmpty { exit 1, "Cannot find any files matching ${params.input}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!" }
  .into {alignment_sort_ch ; alignment_flagstat_ch }


/*
 * STEP 1 - Collate reads and extract FASTQ
 */
if (params.cram) {
process extractFastqCram {

  publishDir params.outdir, mode: 'copy'

  input:
  set val(name), file(alignment) from alignment_sort_ch
  file(reference) from reference_ch

  output:
  set val(name), file('*.fastq.gz') into fastq_ch

  script:
  """
  samtools view -b -T ${reference} ${alignment} |
  samtools collate -Ou - ${params.tmpdir}/${name} | \
  bamtofastq \
    gz=1 \
    F=${name}_R1.fastq.gz \
    F2=${name}_R2.fastq.gz \
    S=${name}_S.fastq.gz \
    O=${name}_U1.fastq.gz \
    O2=${name}_U2.fastq.gz
  for file in *.fastq.gz
  do
    check=`zcat \${file} | head -n 1 | wc -l`
    if [ \$check -ne '1' ]
    then
      rm \$file
    fi
  done
  """
}
} else {
process extractFastqBam {

  publishDir params.outdir, mode: 'copy'

  input:
  set val(name), file(alignment) from alignment_sort_ch

  output:
  set val(name), file('*.fastq.gz') into fastq_ch

  script:
  """
  samtools collate -Ou ${alignment} ${params.tmpdir}/${name} | \
  samtools fastq -1 ${name}_R1.fastq.gz -2 ${name}_R2.fastq.gz -s ${name}_S.fastq.gz
  for file in *.fastq.gz
  do
    check=`zcat \${file} | head -n 1 | wc -l`
    if [ \$check -ne '1' ]
    then
      rm \$file
    fi
  done
  """
}
}

/*
 * STEP 2 - Get read stats from input alignment file
 */
process readStats {

  publishDir params.outdir, mode: 'copy'

  input:
  set val(name), file(alignment) from alignment_flagstat_ch

  output:
  file("*.flagstat") into flagstat_results_ch
  file("*.flagstat.count") into flagstat_counts_ch

  script:
  """
  samtools flagstat ${alignment} > ${name}.flagstat

  flagstat_count1=`grep read1 ${name}.flagstat | awk '{ print \$1 + \$3 }'`
  flagstat_count2=`grep read2 ${name}.flagstat | awk '{ print \$1 + \$3 }'`

  echo \${flagstat_count1} \${flagstat_count2} > ${name}.flagstat.count
  """
}

/*
 * STEP 3 - FastQC on the output FASTQ files.
 */
process fastQC {

  publishDir params.outdir, mode: 'copy'

  input:
  set val(name), file(reads) from fastq_ch

  output:
  file("*.zip") into fastqc_results_ch
  file("*.fastqc.count") into fastqc_counts_ch

  script:
  """
  fastqc -t ${task.cpus} --noextract ${reads}

  fastqc_count1=`unzip -p ${name}_R1_fastqc.zip ${name}_R1_fastqc/fastqc_data.txt | grep 'Total Sequences' | awk '{ print \$3}'`
  fastqc_count2=`unzip -p ${name}_R2_fastqc.zip ${name}_R2_fastqc/fastqc_data.txt | grep 'Total Sequences' | awk '{ print \$3}'`

  echo \${fastqc_count1} \${fastqc_count2} > ${name}.fastqc.count
  """
}

/*
 * STEP 4 - Collect all the counts and check them
 */
process checkCounts {

  publishDir params.outdir, mode: 'copy'

  input:
  file(flagstat_counts) from flagstat_counts_ch.collect()
  file(fastqc_counts) from fastqc_counts_ch.collect()

  output:
  file("read_counts.txt") into collect_counts_ch

  script:
  """
  collect_counts.pl > read_counts.txt
  """
}
