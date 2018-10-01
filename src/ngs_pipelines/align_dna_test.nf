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
output_dir = file(params.outputDir)

// Duplicate the read pairs into one queue for runFastQCOriginal
// and runTrimming.
(readPairsFastQCOriginal, readPairsRunTrimming) = readPairs.separate(2) { x -> [x, x] }

// --------------------------------------------------------------------------
// Step 1a) Run FastQC
//
// - yields report
// --------------------------------------------------------------------------


fastqc_output_dir = output_dir + "/fastqc/"
println fastqc_output_dir
process runFastQCOriginal {
	
	publishDir fastqc_output_dir
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

