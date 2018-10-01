#!/usr/bin/env nextflow

import CopyHelper
import ChannelUtil
import FastQC
import PathUtil
import ParamsHelper

// ---------------------------------------------------------------------------
// Read preprocessing and alignment for DNA (WES or WGS) reads.
// ---------------------------------------------------------------------------

if (params.verbose)
    echo true

ParamsHelper.checkNonEmptyParam(params.runID, "runID");
ParamsHelper.checkNonEmptyParam(params.runPlatform, "runPlatform");
ParamsHelper.checkNonEmptyParam(params.dataDir, "dataDir");

copyHelper = new CopyHelper(params.dataDir, params.printCopyMsgs)

// Open channel for left and right files and merge it into triples, the
// first entry is the LCS of the file names that can be used as a read
// pair identifier.
readPairs = ChannelUtil.createFilePairChannel(
        params.runID,
        Channel.fromPath([params.dataDir, 'fastq', 'original', '*_{R,}1.fastq.gz'].join(File.separator)),
        Channel.fromPath([params.dataDir, 'fastq', 'original', '*_{R,}2.fastq.gz'].join(File.separator)),
        )

// Genome and index files.
indexFileBWA = file(params.indexBWA)

// Output Directory
output_dir = dir(params.outputDir)

// Duplicate the read pairs into one queue for runFastQCOriginal
// and runTrimming.
(readPairsFastQCOriginal, readPairsRunTrimming) = readPairs.separate(2) { x -> [x, x] }

// --------------------------------------------------------------------------
// Step 1a) Run FastQC
//
// - yields report
// --------------------------------------------------------------------------


fastqc_output_dir = "${output_dir}/fastqc/"
process runFastQCOriginal {
    cpus params.fastqc.cpus
    module 'fastqc/0.11.2'

    input:
    set runID, file(readL), file(readR) from readPairsFastQCOriginal

    output:
    set file('*.zip'), file('*.html') into fastqcOutputOriginal

    script:
    """
    set -x
    fastqc -t params.fastqc.cpus -o . ${readL} ${readR}
    """
}

copyHelper.copyFiles(fastqcOutputOriginal, fastqc_output_dir);

System.exit(0);

// --------------------------------------------------------------------------
// Step 1b) Run adapter trimming
//
// - yields trimmed read, used as downstream input
// --------------------------------------------------------------------------

trimmomatic_output_dir = output_directory + "trimmomatic/"

process runTrimming {
    cpus params.trimmomatic.cpus

    input:
    set runID, file(readL), file(readR) from readPairsRunTrimming

    output:
    set runID, file { "out/${readL}" }, file { "out/${readR}" } into readPairsTrimmed
    set file("*.log") into trimmingLogs

    script:
    """
    print "running trimmomatic"
    set -x
    # call trimmomatic
    NAMEBASE=${readL}
    trimmomatic \\
	PE \\
	-threads num_threads \\
	${readL} ${readR} pe \\
	-z \\
	-t ${params.skewer.cpus} \\
	${readL} \\
	${readR} \\
    \${NAMEBASE%.gz}-trimmed-pair1.fastq.gz
	\${NAMEBASE%.gz}-trimmed-pair2.fastq.gz
	LEADING: min_base_quality \\
	TRAILING: min_base_quality \\
	SLIDINGWINDOW: sliding_window \\
	MINLEN: min_read_length
    """
}

// Duplicate the read pairs into multiple queues for processing / copying out.
(readPairsFastQCTrimmed,
 readPairsRunMapping,
 readPairsTrimmedCopyOut) = readPairsTrimmed.separate(3) { x -> [ x, x, x ] }

// Copy out results from trimming step (map removes the pair).
copyHelper.copyFiles(trimmingLogs, 'reports/trimming');
copyHelper.copyFiles(readPairsTrimmedCopyOut.map { [it[1], it[2]] }, 'fastq/trimmed');

// --------------------------------------------------------------------------
// Step 2a) Run FastQC on trimmed
//
// - yields report
// --------------------------------------------------------------------------

process runFastQCTrimmed {
    cpus params.fastqc.cpus
    module 'fastqc/0.11.2'

    input:
    set runID, file(readL), file(readR) from readPairsFastQCTrimmed

    output:
    set file('*.zip'), file('*.html') into fastqcOutputTrimmed

    script:
    """
    set -x
    fastqc -t ${params.fastqc.cpus} -o . ${readL} ${readR}
    """
}

copyHelper.copyFiles(fastqcOutputTrimmed, 'reports/fastqc-trimmed');

// --------------------------------------------------------------------------
// Step 2b) Align reads using BWA-MEM
//
// - align reads
// - sort
// - mark duplicates
// - yields alignment for downstream processing
// --------------------------------------------------------------------------

// Group trimmed read FASTQ files by runID (the runID is part of the output
// of a previous process).
jointBams = readPairsRunMapping.map{f -> [f[0], f[1], f[2]] }.groupTuple()

// The alignments are written to the temporary files alignment.bam. These
// BAM files are already sorted.
process runReadMapping {
    cpus params.bwa.cpus
 
    input:
    indexFileBWA
    set runID, readL, readR from jointBams

    output:
    file { "${runID}.bam*" } into bamFilesOut
    set runID, file { "${runID}.bam" }, file { "${runID}.bam.bai" } into bamFiles

    script:
    """
    
	print "beginning alignment"
	print fastq_r1, fastq_r2
	cmd = BwaMem[:]
	cmd.append("-M")
	cmd.append("-t")
	cmd.append("6")
	cmd.append("-R")
	cmd.append(read_groups)
	#cmd.append(log)
	cmd.append(bwa_index)
	cmd.append(fastq_r1_location) 
	cmd.append(fastq_r2_location)
	print " ".join(cmd)
	f = open(bwa_sam, "w")
	errcode = subprocess.call(cmd, stdout=f)
	if(errcode == 0):
		print "BWA finished successfully"
	else:
		print "BWA failed !!!!"
		del steps_to_process[:]
	f.close()
	del cmd
		
    set -x
    bwa mem \\
        -R '@RG\tID:${runID}\tSM:${runID}\tPL:${params.runPlatform}' \\
        -t ${params.bwa.cpus} \\
        ${indexFileBWA} \\
        <(zcat ${readL.join(" ")}) \\
        <(zcat ${readR.join(" ")}) \\
        | samblaster \\
        | samtools view -u -Sb - \\
        | samtools sort - ${runID}
    samtools index ${runID}.bam
    """
}

(bamFilesForCoverage,
 bamFilesForQualimap,
 bamFilesForVariantCalling) = bamFiles.separate(3) { x -> [ x, x, x ] }

copyHelper.copyFiles(bamFilesOut, 'bam')
