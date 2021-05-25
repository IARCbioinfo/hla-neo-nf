#!/usr/bin/env nextflow


//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run iarcbioinfo/purple-nf -singularity [OPTIONS]

    Mandatory arguments:
      --tn_file		         [file] File containing list of T/N bam/cram files to be processed
      --cohort_dir         [dir]  directory where the BAM or CRAM  file are stored
      --ref               [file] fasta file of chr6 of reference genome [chr6-hg38.fa], should be indexed [chr6-hg38.fa.fai]
    Optional arguments:
      --output_folder       [string] name of output folder
      --cpu                 [Integer]  Number of CPUs[def:2]
      --mem 		            [Integer] Max memory [def:8Gb]
      """.stripIndent()
}





//we display help information
if (params.help){ show_help(); exit 0;}
//we display the header of the tool
log.info IARC_Header()
log.info tool_header()
//Check mandatory parameters
assert (params.ref != null) : "please specify --ref chr6-hg38.fasta"
assert (params.tn_file != null ) : "please specify --tn_file"
assert (params.cohort_dir != null ) : "please specify --cohort_dir"

//function that read the tumors to process from a tn_file
if(params.tn_file){
  def cram = params.bam ? false:true
 tn_pairs = parse_tn_file(params.tn_file,params.cohort_dir,cram)
 //we duplicate the tn_pairs channel
 tn_pairs.into {xHLA_input}
}

//chanel for reference genome
ref_fasta = Channel.value(file(params.ref)).ifEmpty{exit 1, "reference file not found: ${params.ref}"}
ref_fai = Channel.value(file(params.ref+'.fai')).ifEmpty{exit 1, "index file not found: ${params.ref}.fai"}
//BWA indexes for re-mapping MHC reads
ref_sa  = file(params.ref+'.sa')
ref_bwt =  file(params.ref+'.bwt')
ref_ann =  file(params.ref+'.ann')
ref_amb =  file(params.ref+'.amb')
ref_pac =  file(params.ref+'.pac')

print_params()

//PATHS in the container for databases
// /hmftools/hg38

// /hmftools/hg38/DiploidRegions.38.bed
// /hmftools/hg38/GC_profile.1000bp.38.cnp
// /hmftools/hg38/GermlineHetPon.38.vcf

process xHLA {

 cpus params.cpu
 memory params.mem+'G'

  publishDir params.output_folder+'/xHLA/', mode: 'copy'
  input:
  set val(tumor_id), file(tumor), file(tumor_index), file(normal), file(normal_index) from xHLA_input
  file(ref) from ref_fasta
  file(fai) from ref_fai
  file ref_sa
  file ref_bwt
  file ref_ann
  file ref_amb
  file ref_pac

  output:
  set val(tumor_id), path("${tumor_id}_xHLA") into xHLA
  script:
       """
      #we run for the normal CRAM it create prefix.mhc.bam
      echo perl ${baseDir}/scripts/extract_mhc_reads_hg38alt.pl -a ${baseDir}/db/hla_regions.lst -b ${normal} -r ${params.ref} -p ${tumor_id}
      #we run xHLA
      echo run.py  --sample_id ${tumor_id} --input_bam_path ${tumor_id}.mhc.bam --output_path ${tumor_id}_mhc
      #we run for the tumor CRAM
      #perl ${baseDir}/scripts/extract_mhc_reads_hg38alt.pl -a ${baseDir}/db/hla_regions.lst -b ${normal} -r ${ref_fasta}
      mkdir ${tumor_id}_xHLA
      """
}

// process AMBER {
//
//  cpus params.cpu
//  memory params.mem+'G'
//
//
//
//   publishDir params.output_folder+'/AMBER/', mode: 'copy'
//   input:
//   set val(tumor_id), file(tumor), file(tumor_index), file(normal), file(normal_index) from tn_pairs_amber
//   file(ref) from ref_fasta
//   file(fai) from ref_fai
//   output:
//   set val(tumor_id), path("${tumor_id}_AMBER") into amber
//   script:
//      if(params.tumor_only){
//        """
//       AMBER  -Xms1g -Xmx${params.mem}g  -loci /hmftools/hg38/GermlineHetPon.38.vcf -ref_genome ${ref} -tumor_only \\
//               -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_AMBER -threads ${params.cpu}
//       """
//      }else{
//        """
//        AMBER   -Xms1g -Xmx${params.mem}g -loci /hmftools/hg38/GermlineHetPon.38.vcf -ref_genome ${ref} \\
//                -reference ${tumor_id}_N -reference_bam ${normal}  \\
//                -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_AMBER -threads ${params.cpu}
//         """
//      }
//
// }
// //we merge previous results from amber and cobalt
// amber_cobalt=amber.join(cobalt, remainder: true)
//
// process PURPLE {
//
//  cpus params.cpu
//  memory params.mem+'G'
//
//   publishDir params.output_folder+'/PURPLE/', mode: 'copy'
//
//   input:
//   //set val(tumor_id), file(tumor), file(tumor_index), file(normal), file(normal_index) from tn_pairs_amber
//   set val(tumor_id), path(amber_dir), path(cobalt_dir) from amber_cobalt
//   file(ref) from ref_fasta
//   file(fai) from ref_fai
//   file(dict) from ref_dict
//   output:
//   set val(tumor_id), path("${tumor_id}_PURPLE") into purple
//   //MESO_071_T_T.purple.purity.tsv
//   //set val(tumor_id), file("${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.tsv") into stats_purple
//   file("${tumor_id}_T.purple.purity.sample.tsv") into stats_purple
//
//   script:
//      if(params.tumor_only){
//        """
//        PURPLE  -Xms1g -Xmx${params.mem}g -tumor_only  -tumor ${tumor_id}_T \\
//                -no_charts \\
//                -output_dir ${tumor_id}_PURPLE \\
//                -amber ${amber_dir} \\
//                -cobalt ${cobalt_dir} \\
//                -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \\
//                -threads ${params.cpu} \\
//                -ref_genome ${ref}
//
//         awk -v tumor=${tumor_id} '{print tumor"\t"\$0}' ${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.tsv > ${tumor_id}_T.purple.purity.sample.tsv
//
//        """
//      }else{
//        """
//         PURPLE  -Xms1g -Xmx${params.mem}g -reference ${tumor_id}_N  -tumor ${tumor_id}_T \\
//                -no_charts \\
//                -output_dir ${tumor_id}_PURPLE \\
//                -amber ${amber_dir} \\
//                -cobalt ${cobalt_dir} \\
//                -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \\
//                -threads ${params.cpu} \\
//                -ref_genome ${ref}
//
//          awk -v tumor=${tumor_id} '{print tumor"\t"\$0}' ${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.tsv > ${tumor_id}_T.purple.purity.sample.tsv
//        """
//      }
//
// }
//
// stats_purple.collectFile(name: 'purple_summary.txt', storeDir: params.output_folder, seed: 'tumor_id\tpurity\tnormFactor\tscore\tdiploidProportion\tploidy\tgender\tstatus\tpolyclonalProportion\tminPurity\tmaxPurity\tminPloidy\tmaxPloidy\tminDiploidProportion\tmaxDiploidProportion\tversion\tsomaticPenalty\twholeGenomeDuplication\tmsIndelsPerMb\tmsStatus\ttml\ttmlStatus\ttmbPerMb\ttmbStatus\tsvTumorMutationalBurden\n', newLine: false, skip: 1)


/*
*
* Functions to create channels from TSV or directories containing BAM/CRAM
*
*/

//we read the pairs from tn_file
def parse_tn_file (tn_file,path,cram){
	    // FOR INPUT AS A TAB DELIMITED FILE
			def file_ext = cram ? '.crai':'.bai'
			//[sample t[.bam,cram] t[.bai,crai] n[.bam,.cram] n[.bai,.crai]]
    def tn_pairs=Channel.fromPath(tn_file)
      .splitCsv(header: true, sep: '\t', strip: true)
      .map{row -> [ row.tumor_id,
               file(path + "/" + row.tumor),
               file(path + "/" + row.tumor+file_ext),
               file(path + "/" + row.normal),
               file(path + "/" + row.normal+file_ext)]}
      .ifEmpty{exit 1, "${tn_file} was empty - no tumor/normal supplied" }
	//we return the channel
  return tn_pairs
}

// print the calling parameter to the log and a log file
def print_params () {
  //software versions for v2.0
  def software_versions = ['xhla' : '0.0.0',
                          'pVactools'   : 'x.y']
  //we print the parameters
  log.info "\n"
  log.info "-\033[2m------------------Calling PARAMETERS--------------------\033[0m-"
  log.info params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"
  log.info "-\033[2m------------------Software versions--------------------\033[0m-"
  log.info software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"


  //we print the parameters to a log file
   def output_d = new File("${params.output_folder}/nf-pipeline_info/")
   if (!output_d.exists()) {
       output_d.mkdirs()
   }
   def output_tf = new File(output_d, "run_parameters_report.txt")
   def  report_params="------------------Calling PARAMETERS--------------------\n"
        report_params+= params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="\n------------------NEXTFLOW Metadata--------------------\n"
        report_params+="nextflow version : "+nextflow.version+"\n"
        report_params+="nextflow build   : "+nextflow.build+"\n"
        report_params+="Command line     : \n"+workflow.commandLine.split(" ").join(" \\\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="-----------------Software versions--------------------\n"
        report_params+=software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"

   output_tf.withWriter { w -> w << report_params}
}


//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        HLA-NEO: Pipeline to predict neoantigens from WGS data (${workflow.manifest.version})
        """
}

//header for the IARC tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}