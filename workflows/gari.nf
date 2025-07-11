/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowGari.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULES
//
include { DOWNLOAD_KRAKENDB } from '../modules/local/downloadKrakenDB.nf'
include { SKESA } from '../modules/local/skesa'
include { STAT_SUMMARY } from '../modules/local/createMetaReport'
include { STAT_SUMMARY_QC } from '../modules/local/createMetaReport_QC'
include { BBMAP_ALIGN } from '../modules/local/bbmap_cov'
include { BBMAP_RENAME } from '../modules/local/bbmap_rename'
include { CREATE_REPORT } from '../modules/local/summarizeReports'
include { CHECKM_LINEAGEWF } from '../modules/local/checkM/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_READ } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_ASM } from '../modules/nf-core/kraken2/kraken2/main'
include { SPADES } from '../modules/nf-core/spades/main'
include { SHOVILL } from '../modules/nf-core/shovill/main'
include { ASSEMBLYSCAN } from '../modules/nf-core/assemblyscan/main'
include { SKANI_SEARCH } from '../modules/nf-core/skani/search/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow GARI {
    
    ch_versions = Channel.empty()

    threshold_file = file(params.thresholds)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)  

    if (params.kraken_db==null){
        ch_kraken_db = DOWNLOAD_KRAKENDB().krakenDB
        ch_versions = ch_versions.mix(DOWNLOAD_KRAKENDB.out.versions)  
    }
    else{
        ch_kraken_db = params.kraken_db
    }

    if (params.skani_db==null){
        //placeholder for now
        error "Error: no database for skani specified."
    }
    else{
        ch_skani_db = params.skani_db
    }

    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    if (params.qc_mode==false) {
        FASTP (
            INPUT_CHECK.out.reads, 
            [], 
            false, 
            false
        )
        ch_versions = ch_versions.mix(FASTP.out.versions)     

        KRAKEN2_READ (
            FASTP.out.reads, 
            ch_kraken_db, 
            false, 
            false,
            true
        )
        ch_versions = ch_versions.mix(KRAKEN2_READ.out.versions)

        // modify input stream for spades e.g. add empty fields for pb and ont
        fastp_adjust = FASTP.out.reads.map { meta, fastq -> [ meta, fastq, [], [] ] }
        
        if (params.assembler=="shovill") {
            SHOVILL (
                FASTP.out.reads
            )
            ch_versions = ch_versions.mix(SHOVILL.out.versions)
            asm_tmp = SHOVILL.out.contigs
        }
        else if (params.assembler=="skesa"){
            SKESA (
                FASTP.out.reads
            )
            ch_versions = ch_versions.mix(SKESA.out.versions)
            asm_tmp = SKESA.out.contigs
        }
        else {
            SPADES (
                fastp_adjust, 
                [], 
                []
            )
            ch_versions = ch_versions.mix(SPADES.out.versions)
            asm_tmp = SPADES.out.scaffolds
        }
        BBMAP_RENAME (
            asm_tmp,
            params.preset,
            params.min_size
        )
        ch_versions = ch_versions.mix(BBMAP_RENAME.out.versions)

        asm = BBMAP_RENAME.out.rename
        asm_adjust = BBMAP_RENAME.out.rename.map {[ [id: it[0].id + '-ASM', single_end:true, species: it[0].species], it[1] ] }
        asm_adjust2 = BBMAP_RENAME.out.rename.map {[ [id: it[0].id, single_end: it[0].single_end, species: it[0].species, batch_input:true], it[1] ] }        

        bbmapIn_ch = FASTP.out.reads.join(asm)
        BBMAP_ALIGN (
            bbmapIn_ch
        )    
        ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)
    } else {
        BBMAP_RENAME (
            INPUT_CHECK.out.reads,
            params.preset,
            params.min_size
        )
        asm = BBMAP_RENAME.out.rename
        asm_adjust = INPUT_CHECK.out.reads.map {[ [id: it[0].id + '-ASM', single_end:true, species: it[0].species], it[1] ] }
        asm_adjust2 = INPUT_CHECK.out.reads.map {[ [id: it[0].id, single_end: it[0].single_end, species: it[0].species, batch_input:true], it[1] ] }
    }

    ASSEMBLYSCAN (
        asm
    )
    ch_versions = ch_versions.mix(ASSEMBLYSCAN.out.versions)
        
    KRAKEN2_ASM (
        asm_adjust, 
        ch_kraken_db, 
        false, 
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_ASM.out.versions)
        

    SKANI_SEARCH (
        asm_adjust2,
        params.skani_db
    )
    ch_versions = ch_versions.mix(SKANI_SEARCH.out.versions)

    CHECKM_LINEAGEWF (
        asm,
        "fasta",
        params.checkm_db
    )
    ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)

    if (params.qc_mode==false) {
        skani_ch = SKANI_SEARCH.out.search.map{[ [id: it[0].id, single_end:false, species: it[0].species], it[1] ] }
        asmscan_ch = ASSEMBLYSCAN.out.json.map{[ [id: it[0].id, single_end:false, species: it[0].species], it[1]] }
        fastp_ch = FASTP.out.json.map{[ [id: it[0].id, single_end:false, species: it[0].species], it[1] ] }
        krakenR_ch = KRAKEN2_READ.out.report.map{[ [id: it[0].id, single_end:false, species: it[0].species], it[1]] }
        krakenA_ch = KRAKEN2_ASM.out.report.map{[ [id: it[0].id.minus("-ASM"), single_end:false, species: it[0].species], it[1]] }
        bbmap_ch = BBMAP_ALIGN.out.log.map{[ [id: it[0].id, single_end:false, species: it[0].species], it[1]] }
        checkm_ch = CHECKM_LINEAGEWF.out.checkm_tsv.map{[ [id: it[0].id, single_end:false, species: it[0].species], it[1]] }

        concat_ch = asm.join(skani_ch).join(asmscan_ch).join(fastp_ch).join(krakenR_ch).join(krakenA_ch).join(bbmap_ch).join(checkm_ch)


        STAT_SUMMARY (
            concat_ch,
            threshold_file,
            ch_kraken_db,
            ch_skani_db
        )
        ch_versions = ch_versions.mix(STAT_SUMMARY.out.versions)

        CREATE_REPORT (
            STAT_SUMMARY.out.json.collect(),
            params.outdir
        )
    }
    else {
        // needs to be set to single_end since the assemblies are parsed as single_end
        skani_ch = SKANI_SEARCH.out.search.map{[ [id: it[0].id, single_end:true, species: it[0].species], it[1] ] }
        asmscan_ch = ASSEMBLYSCAN.out.json.map{[ [id: it[0].id, single_end:true, species: it[0].species], it[1]] }
        krakenA_ch = KRAKEN2_ASM.out.report.map{[ [id: it[0].id.minus("-ASM"), single_end:true, species: it[0].species], it[1]] }
        checkm_ch = CHECKM_LINEAGEWF.out.checkm_tsv.map{[ [id: it[0].id, single_end:true, species: it[0].species], it[1]] }

        concat_ch = asm.join(skani_ch).join(asmscan_ch).join(krakenA_ch).join(checkm_ch)

        STAT_SUMMARY_QC (
            concat_ch,
            threshold_file,
            ch_kraken_db,
            ch_skani_db
        )
        ch_versions = ch_versions.mix(STAT_SUMMARY_QC.out.versions)        

        CREATE_REPORT(
            STAT_SUMMARY_QC.out.json.collect(), 
            params.outdir
        )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
