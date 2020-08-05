#!/usr/bin/env nextflow
/*
========================================================================================
                         CLOMP
========================================================================================
 CLinically Okay Metagenomic Pipeline
 #### Homepage / Documentation
 https://github.com/FredHutch/CLOMP
----------------------------------------------------------------------------------------
*/

// Using the Nextflow DSL-2 to account for the logic flow of this workflow
nextflow.preview.dsl=2


def helpMessage() {
    log.info"""
    A entirely open source, fast, accurate, multi-sample and highly configurable end to end metagenomics pipeline for the people.

    Usage:

    An example command for running the pipeline is as follows:

    nextflow run FredHutch/CLOMP \\
        --INPUT_FOLDER reads/ \\
        --PAIRED_END \\
        --HOST_FILTER_FIRST \\
        --OUTDIR output/
        

    Required Arguments:
      --INPUT_FOLDER    Folder containing input files (FASTQ format)
      --ENTREZ_EMAIL
                Email adress requried for interfacing with NCBI servers. Please provide your own email address.
      --OUTDIR
                The output directory where the results will be saved

    Optional Arguments:
      --PAIRED_END      If specified, look for paired-end files with { _R1_ / _R2_ }
      --INPUT_SUFFIX    Suffix of input files (default: .fastq.gz)
      --SNAP_INDEXES
                Comma-delimited list of SNAP indexed references
      --SNAP_OPTIONS
                Options used to run SNAP (must enclose in quotes)
                default: -mrl 65 -d 9 -h 30000 -om 1 -omax 20 -map
      --SNAP_BATCHSIZE
                Number of samples to align in parallel over each SNAP index shard, in a batch
      --HOST_FILTER_FIRST
                If specified, perform host filtering prior to trimming
      --SECOND_PASS
                Should we do a more sensitive second pass of host filtering? (default: false)
      --TRIMMOMATIC_OPTIONS
                Options used to run Trimmomatic (default: ':2:30:10 HEADCROP:10 SLIDINGWINDOW:4:20 CROP:65 MINLEN:65')
      --BBDUK_TRIM_OPTIONS
                Options to be used to run BBDuk trimming (defaults: 'ktrim=r k=27 hdist=1 edist=0 mink=4 qtrim=rl trimq=6 minlength=65 ordered=t qin=33')
      --TRIMMOMATIC_JAR_PATH
                Path to Trimmomatic executable (default: "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/trimmomatic-0.38.jar")
      --TRIMMOMATIC_ADAPTER_PATH
                Path to Trimmomatic adapters (default: "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/adapters.fa")
      --BWT_DB_PREFIX   
                Prefix for the files in the host genome directory (default: hg38)
      --BWT_DB_LOCATION 
                Host genome directory (default: s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/hg38/)
      --BWT_SECOND_PASS_OPTIONS 
                Bowtie2 options used in the second pass analysis (default: '-D 35 -R 5 -N 1 -L 19 -i S,1,0.50')
      --KRAKEN_DB_PATH
                The directory containing the Kraken database (default: s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/kraken_db/)
      --LOGIC
                This controls the different types of tie-breaking logic that are baked in to CLOMP: strict/90/oneoff (default: strict)
      --BLAST_EVAL
                E-value cuttoff for BLAST (default: 0.001)
      --BLAST_CHECK
                Check results with BLAST (default: false)
      --BLAST_CHECK_DB
                Database used to check results with BLAST (default: None)
      --FILTER_LIST
                List of taxids to filter from tiebreaking (default: "[12908,28384,48479]")
      --KRAKEN_DB_PATH
                Path to Kraken database (default: "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/kraken_db/")
      --H_STRICT
                Assign reads with any host alignment to host (default: false)
      --H_TAXID
                Host taxid (default: 9606)
      --INCLUSION_TAXID
                Inclusion TAXID will pull all reads that get classified below this level but not downstream of the exclusion taxid (default: 2759)
      --EXCLUSION_TAXID
                See above (default: 9604)
      --BASE_DELIMITER
                Each sample needs a unique name and everything BEFORE the FIRST instance of this string will be the sample base name (default: '_')
      --MIN_READ_CUTOFF
                The minimum number of reads assigned to a species level taxid that are required to build the sam file (default: 10)
      --SAM_NO_BUILD_LIST
                Taxids that contain these taxids in their lineage will not be built (default: "[2759,77133]")
      --ADD_HOST_FILTERED_TO_REPORT
                Should we add the number of reads called to the host to the final Pavian output? (default: true)
      --HOST_FILTER_TAXID
                The taxid for which we should classify host filtered reads as when adding them to the Pavian report (default: 9606)
      --WRITE_UNIQUES
                Should program write a new folder named after the sample and write all unique reads that mapped to each taxid into this folder? VERY useful for debugging tie-breaking (default: true)
      --EDIT_DISTANCE_OFFSET
                The edit distance offset is the maximum difference in edit distance we will accept for alignments (default: 6)
      --BUILD_SAMS
                Should we build SAM files aligning reads to species level assignments? (default: false)
      --TIEBREAKING_CHUNKS
                The number of chunks to process in parallel for tiebreaking
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Set defaults
// These values are overridden by what the user specifies (e.g. with --R1)
params.INPUT_FOLDER = false
params.INPUT_SUFFIX = ".fastq.gz"
params.PAIRED_END = false
params.OUTDIR = false
params.SNAP_INDEXES = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.00/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.01/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.02/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.03/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.04/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.05/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.06/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.07/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.08/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.09/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.10/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.11/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.12/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.13/" 
params.SNAP_OPTIONS = "-mrl 65 -d 9 -h 30000 -om 1 -omax 20"
params.HOST_FILTER_FIRST = false
params.SECOND_PASS = false
params.TRIMMOMATIC_OPTIONS = ':2:30:10 HEADCROP:10 SLIDINGWINDOW:4:20 CROP:65 MINLEN:65'
params.BBDUK_TRIM_OPTIONS = 'ktrim=r k=27 hdist=1 edist=0 mink=4 qtrim=rl trimq=6 minlength=65 ordered=t qin=33'
params.TRIMMOMATIC_JAR_PATH = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/trimmomatic-0.38.jar"
params.TRIMMOMATIC_ADAPTER_PATH = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/adapters.fa"
params.SEQUENCER = 'ILLUMINACLIP:'
params.BWT_DB_PREFIX = "hg38"
params.BWT_DB_LOCATION = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/hg38/"
params.BWT_SECOND_PASS_OPTIONS = '-D 35 -R 5 -N 1 -L 19 -i S,1,0.50'
params.BLAST_EVAL = 0.001
params.BLAST_CHECK = false
params.DEDUPE = true
params.BLAST_CHECK_DB = false
params.FILTER_LIST = "[12908,28384,48479]"
params.KRAKEN_DB_PATH = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/kraken_db/"
params.H_STRICT = false
params.H_TAXID = 9606
params.LOGIC = "strict"
params.INCLUSION_TAXID = 2759
params.EXCLUSION_TAXID = 9604
params.ENTREZ_EMAIL = "uwvirongs@gmail.com"
params.BASE_DELIMITER = "_"
params.MIN_READ_CUTOFF = 10
params.SAM_NO_BUILD_LIST = "[2759,77133]"
params.ADD_HOST_FILTERED_TO_REPORT = true
params.HOST_FILTER_TAXID = 9606
params.WRITE_UNIQUES = true
params.EDIT_DISTANCE_OFFSET = 6
params.BUILD_SAMS = false
params.SNAP_BATCHSIZE = 20
params.TIEBREAKING_CHUNKS = 2
params.TAXDUMP_NODES = 's3://clomp-reference-data/tool_specific_data/CLOMP/clomp_viz/nodes.dmp'
params.TAXDUMP_MERGED = 's3://clomp-reference-data/tool_specific_data/CLOMP/clomp_viz/merged.dmp'


// Check to make sure that the required parameters have been set
if (!params.INPUT_FOLDER){ exit 1, "Must provide folder containing input files with --INPUT_FOLDER" }
if (!params.OUTDIR){ exit 1, "Must provide output directory with --OUTDIR" }
if (!params.KRAKEN_DB_PATH){ exit 1, "Must provide Kraken database with --KRAKEN_DB_PATH" }

// Make sure the output directory has a trailing slash
if (!params.OUTDIR.endsWith("/")){
    params.OUTDIR = "${params.OUTDIR}/"
}

// Identify some resource files
TRIMMOMATIC_JAR = file(params.TRIMMOMATIC_JAR_PATH)
TRIMMOMATIC_ADAPTER = file(params.TRIMMOMATIC_ADAPTER_PATH)
GENERATE_SUMMARY_SCRIPT = file("modules/summarize_run.r")
SAM_SPLIT = file("${workflow.projectDir}/bin/sam_split.py")
BASH_SPLIT = file("${workflow.projectDir}/bin/split.sh")
TAXDUMP_NODES = file(params.TAXDUMP_NODES)
TAXDUMP_MERGED = file(params.TAXDUMP_MERGED)

if (params.BLAST_CHECK){
    if (!params.BLAST_CHECK_DB){ exit 1, "Must provide BLAST check database with --BLAST_CHECK_DB" }
    BLAST_CHECK_DB = file(params.BLAST_CHECK_DB)
} else {
    BLAST_CHECK_DB = false
}
KRAKEN_DB = [
    file("${params.KRAKEN_DB_PATH}database.idx"),
    file("${params.KRAKEN_DB_PATH}database.kdb"),
    file("${params.KRAKEN_DB_PATH}seqid2taxid.map"),
    file("${params.KRAKEN_DB_PATH}taxDB")
]



/*
 * Import the processes used in this workflow
 */

include collect_snap_results from './modules/clomp_modules'
include validate_single from './modules/clomp_modules'
include validate_paired from './modules/clomp_modules'
include trimmomatic_single from './modules/clomp_modules' params(
    SEQUENCER: params.SEQUENCER, 
    TRIMMOMATIC_OPTIONS: params.TRIMMOMATIC_OPTIONS
)
include trimmomatic_paired from './modules/clomp_modules' params(
    SEQUENCER: params.SEQUENCER, 
    TRIMMOMATIC_OPTIONS: params.TRIMMOMATIC_OPTIONS
)
include filter_human_single from './modules/clomp_modules' params(
    BWT_SECOND_PASS_OPTIONS: params.BWT_SECOND_PASS_OPTIONS, 
    BWT_DB_PREFIX: params.BWT_DB_PREFIX
)
include filter_human_paired from './modules/clomp_modules' params(
    BWT_SECOND_PASS_OPTIONS: params.BWT_SECOND_PASS_OPTIONS, 
    BWT_DB_PREFIX: params.BWT_DB_PREFIX
)
include filter_human_single as filter_human_single_second_pass from './modules/clomp_modules' params(
    BWT_SECOND_PASS_OPTIONS: params.BWT_SECOND_PASS_OPTIONS, 
    BWT_DB_PREFIX: params.BWT_DB_PREFIX
)
include filter_human_paired as filter_human_paired_second_pass from './modules/clomp_modules' params(
    BWT_SECOND_PASS_OPTIONS: params.BWT_SECOND_PASS_OPTIONS, 
    BWT_DB_PREFIX: params.BWT_DB_PREFIX
)
include bbMask_Single from './modules/clomp_modules' params(BBDUK_TRIM_OPTIONS: params.BBDUK_TRIM_OPTIONS)
include deduplicate from './modules/clomp_modules'
include snap_single from './modules/clomp_modules' params(SNAP_OPTIONS: params.SNAP_OPTIONS)
include snap_paired from './modules/clomp_modules' params(SNAP_OPTIONS: params.SNAP_OPTIONS)
include summarize_run from './modules/clomp_modules'
include CLOMP_summary from './modules/clomp_modules' params(
    BLAST_CHECK: params.BLAST_CHECK,
    BLAST_EVAL: params.BLAST_EVAL,
    ADD_HOST_FILTERED_TO_REPORT: params.ADD_HOST_FILTERED_TO_REPORT,
    HOST_FILTER_TAXID: params.HOST_FILTER_TAXID,
    WRITE_UNIQUES: params.WRITE_UNIQUES,
    FILTER_LIST: params.FILTER_LIST,
    H_STRICT: params.H_STRICT,
    H_TAXID: params.H_TAXID,
    LOGIC: params.LOGIC,
    INCLUSION_TAXID: params.INCLUSION_TAXID,
    EXCLUSION_TAXID: params.EXCLUSION_TAXID,
    ENTREZ_EMAIL: params.ENTREZ_EMAIL,
    BASE_DELIMITER: params.BASE_DELIMITER,
    MIN_READ_CUTOFF: params.MIN_READ_CUTOFF,
    SAM_NO_BUILD_LIST: params.SAM_NO_BUILD_LIST,
    EDIT_DISTANCE_OFFSET: params.EDIT_DISTANCE_OFFSET,
    BUILD_SAMS: params.BUILD_SAMS,
    SECOND_PASS: params.SECOND_PASS,
)
include generate_report from './modules/clomp_modules'

// Run the CLOMP workflow
workflow {

    BWT_FILES = Channel
        .fromPath("${params.BWT_DB_LOCATION}/${params.BWT_DB_PREFIX}*")
        .ifEmpty { error "Cannot find any files in ${params.BWT_DB_LOCATION} starting with ${params.BWT_DB_PREFIX}" }
        .map{it -> file(it)}
        .collect()

    SNAP_INDEXES_CH = Channel
        .from(params.SNAP_INDEXES.split(","))
        .map { it -> file(it) }
        .ifEmpty { error "Please specify SNAP indexes with --SNAP_INDEXES" }

    if ( params.PAIRED_END ){
        input_read_ch = Channel
            .fromFilePairs("${params.INPUT_FOLDER}*_R{1,2}*${params.INPUT_SUFFIX}")
            .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT_FOLDER} ending with ${params.INPUT_SUFFIX}" }
            .map { it -> [it[0], it[1][0], it[1][1]]}

        // Validate that the inputs are paired-end gzip-compressed FASTQ
        // This will also enforce that all read pairs are named ${sample_name}.R[12].fastq.gz
        validate_paired(
            input_read_ch
        )

        if ( params.HOST_FILTER_FIRST ){
            filter_human_paired(
                validate_paired.out,
                BWT_FILES
            )
            trimmomatic_paired(
                filter_human_paired.out[0],
                TRIMMOMATIC_JAR,
                TRIMMOMATIC_ADAPTER
            )
            if ( params.SECOND_PASS ){
                filter_human_paired_second_pass(
                    trimmomatic_paired.out,
                    BWT_FILES
                )
                snap_paired(
                    filter_human_paired_second_pass.out.collate(params.SNAP_BATCHSIZE),
                    SNAP_INDEXES_CH
                )
            } else {
                snap_paired(
                    trimmomatic_paired.out.collate(params.SNAP_BATCHSIZE),
                    SNAP_INDEXES_CH
                )
            }
        } else {
            trimmomatic_paired(
                validate_paired.out,
                TRIMMOMATIC_JAR,
                TRIMMOMATIC_ADAPTER
            )
            filter_human_paired(
                trimmomatic_paired.out,
                BWT_FILES
            )
            snap_paired(
                filter_human_paired.out[0].collate(params.SNAP_BATCHSIZE),
                SNAP_INDEXES_CH
            )
        }

        collect_snap_results(
            snap_paired.out.flatten().map{
                it -> [it.name.split("__")[0], it]
            }.groupTuple()
        )

        CLOMP_summary(
            collect_snap_results.out.join(
                filter_human_paired.out[1].map{
                    it -> [it.name.split(".log")[0], it]
                }
            ),
            BLAST_CHECK_DB,
            KRAKEN_DB
        )
    } else {
        input_read_ch = Channel
            .fromPath("${params.INPUT_FOLDER}*${params.INPUT_SUFFIX}")
            .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT_FOLDER} ending with ${params.INPUT_SUFFIX}" }
            .map { it -> [it.name.replace(/${params.INPUT_SUFFIX}/, ""), file(it)]}

        // Validate that the inputs are single-end gzip-compressed FASTQ
        // This will also enforce that all read pairs are named ${sample_name}.R1.fastq.gz
        validate_single(
            input_read_ch
        )

        if ( params.HOST_FILTER_FIRST ){
            filter_human_single(
                validate_single.out,
                BWT_FILES
            )
            trimmomatic_single(
                filter_human_single.out[0],
                TRIMMOMATIC_JAR,
                TRIMMOMATIC_ADAPTER
            )
            if ( params.SECOND_PASS ){
                filter_human_single_second_pass(
                    trimmomatic_single.out,
                    BWT_FILES
                )
                snap_single(
                    filter_human_single_second_pass.out[0].collate(params.SNAP_BATCHSIZE),
                    SNAP_INDEXES_CH
                )
            } else {
                snap_single(
                    trimmomatic_single.out.collate(params.SNAP_BATCHSIZE),
                    SNAP_INDEXES_CH
                )
            }
        } else {
            //Trimmomatic depracated. Using BBDuk for quality filtering, adapter trimming, and entropy filtering 

            trimmomatic_single(
              validate_single.out,
              TRIMMOMATIC_JAR,
              TRIMMOMATIC_ADAPTER
            )
            // bbMask_Single(
            //     validate_single.out,
            //     TRIMMOMATIC_ADAPTER
            //     )
            if ( params.DEDUPE){ 
            deduplicate(
                //bbMask_Single.out
                trimmomatic_single.out
                )
            filter_human_single(
                deduplicate.out,
                BWT_FILES
                )
            }
            else { 
            filter_human_single(
                //bbMask_Single.out,
                trimmomatic_single.out,
                BWT_FILES
            )

            }
            snap_single(
                filter_human_single.out[0].toSortedList().flatten().collate(params.SNAP_BATCHSIZE),
                SNAP_INDEXES_CH
                
            )
        }

        collect_snap_results(
            snap_single.out.flatten().map{
                it -> [it.name.split("__")[0], it]
            }.groupTuple(),
            SAM_SPLIT,
            BASH_SPLIT
        )

        CLOMP_summary(
            collect_snap_results.out.transpose(),
            BLAST_CHECK_DB,
            KRAKEN_DB
        )
        generate_report(
            CLOMP_summary.out[0].groupTuple(
            ).join(
                CLOMP_summary.out[1].groupTuple()
            ).join(
                CLOMP_summary.out[2].groupTuple()
            ),
            BLAST_CHECK_DB,
            KRAKEN_DB,
            deduplicate.out,
            TAXDUMP_NODES,
            TAXDUMP_MERGED
        )
        // summarize_run( 
        //     generate_report.out[0].toList(), 
        //         generate_report.out[1].toList(), 
        //         generate_report.out[2].toList(),
        //         GENERATE_SUMMARY_SCRIPT
        // )
    }    
    publish:
    generate_report.out to: "${params.OUTDIR}"
        //summarize_run.out to: "${params.OUTDIR}"
        //filter_human_single.out[1] to: "${params.OUTDIR}/logs/"
}
