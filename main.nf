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
      --cram_dir         [dir]  directory where the BAM or CRAM  file are stored
      --vcf_dir         [dir]  directory where the VCF files are stored
      --vep_dir         [dir] directory containing VEP database for annotation [hg38, GENCODE 33]
      --ref               [file] fasta file of chr6 of reference genome [chr6-hg38.fa], should be indexed [chr6-hg38.fa.fai]
    Optional arguments:
      --output_folder       [string] name of output folder
      --cpu                 [Integer]  Number of CPUs[def:2]
      --mem 		            [Integer] Max memory [def:8Gb]
      --pvactools_predictors [string] predicttions tools to compute neoantigens [def:all_class_i,all_class_ii or NetMHCpan,NetMHCIIpan]
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
assert (params.cram_dir != null ) : "please specify --cram_dir"
assert (params.vcf_dir != null ) : "please specify --vcf_dir"
assert (params.vep_dir != null ) : "please specify --vep_dir"

//function that read the tumors to process from a tn_file
if(params.tn_file){
  def cram = params.bam ? false:true
 tn_pairs = parse_tn_file(params.tn_file,params.vcf_dir,params.cram_dir,cram)
 //we duplicate the tn_pairs channel
 tn_pairs.into {xHLA_input; xVEP_input; pvactools_input}
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
vep_dir_path = file(params.vep_dir)
print_params()


process xHLA {

 cpus params.cpu
 memory params.mem+'G'

  publishDir params.output_folder+'/xHLA/', mode: 'copy'
  input:
  set val(tumor_id), file(vcf), file(normal), file(normal_index), val(normal_id), val(tumor_id_name) from xHLA_input
  file(ref) from ref_fasta
  file(fai) from ref_fai
  file ref_sa
  file ref_bwt
  file ref_ann
  file ref_amb
  file ref_pac

  output:
  set val(tumor_id), file("report-${tumor_id}-hla.json") into xHLA_out
  script:
       """
      # we get the mhc reads for the normal CRAM/BAM it create prefix.mhc.bam
      perl ${baseDir}/scripts/extract_mhc_reads_hg38alt.pl -a ${baseDir}/db/hla_regions.lst -b ${normal} -r ${params.ref} -p ${tumor_id}
      # we run xHLA
      run.py  --sample_id ${tumor_id} --input_bam_path ${tumor_id}.mhc.bam --output_path ${tumor_id}_mhc
      mv ${tumor_id}_mhc/report-${tumor_id}-hla.json .
      """
}

// Annnot the VCF with VEP tools
//create a local VEP database (gencode 33) ~ 16Gb size
//vep_install -a cf -s homo_sapiens -y GRCh38 -c vep-db-99 --CONVERT
//https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/vep.html
process VEP {
 cpus params.cpu
 memory params.mem+'G'

  publishDir params.output_folder+'/VEP/', mode: 'copy'
  input:
  set val(tumor_id), file(vcf), file(normal), file(normal_index), val(normal_id), val(tumor_id_name) from xVEP_input
  file (vep_dir_path)
  output:
    set val(tumor_id), file("${tumor_id}.vep.vcf") into xVEP_out
  script:
       """
        vep -i ${vcf} \\
        -o ${tumor_id}.vep.vcf \\
        --cache --offline \\
        --dir_cache ${vep_dir_path} \\
        --format vcf \\
        --vcf \\
        --symbol  \\
        --terms SO \\
        --tsl \\
        --hgvs \\
        --fasta ${vep_dir_path}/homo_sapiens/99_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \\
        --plugin Frameshift \\
        --plugin Wildtype \\
        --dir_plugins ${baseDir}/VEP_plugins \\
        --pick  --transcript_version
	# we remove the VAF from the VCF
	bcftools annotate -x  FORMAT/AF ${tumor_id}.vep.vcf > ${tumor_id}.vep.noAF.vcf
	mv ${tumor_id}.vep.noAF.vcf ${tumor_id}.vep.vcf
       """
}

//we have to sync the xHLA, VEP ouputs and pvac input
xHLA_xVEP=xHLA_out.join(xVEP_out, remainder: true)
pvac_hla_vep=pvactools_input.join(xHLA_xVEP,remainder: true)

process pVactools {
 cpus params.cpu
 memory params.mem+'G'

  publishDir params.output_folder+'/pVACTOOLS/', mode: 'copy'
  input:
  set val(tumor_id), file(vcf), file(normal), file(normal_index), val(normal_id), val(tumor_id_name), file(hla_dir_out),file(vcf_vep) from pvac_hla_vep

  output:
    set val(tumor_id), path("${tumor_id}*_pvactools") into pVACTOOLS_out
    file("${tumor_id}.pvactools.log")
  script:
       """
         perl ${baseDir}/scripts/pbactools_wrapper.pl -a ${hla_dir_out} \\
            -b ${baseDir}/db/xHLA2PVAC_alleles.txt -c ${normal_id}   -d ${vcf_vep} -t ${tumor_id_name} -p ${tumor_id} \\
            -e ${params.pvactools_predictors} > ${tumor_id}.pvactools.log
       """
}


/*
*
* Functions to create channels from TSV or directories containing BAM/CRAM
*
*/

//we read the pairs from tn_file
def parse_tn_file (tn_file,path_vcf,path_cram,cram){
	    // FOR INPUT AS A TAB DELIMITED FILE
			def file_ext = cram ? '.crai':'.bai'
			//[sample t[.bam,cram] t[.bai,crai] n[.bam,.cram] n[.bai,.crai]]
      //id      vcf     normal_cram     normal_id       tumor_id


    def tn_pairs=Channel.fromPath(tn_file)
      .splitCsv(header: true, sep: '\t', strip: true)
      .map{row -> [ row.id,
               file(path_vcf + "/" + row.vcf),
               file(path_cram + "/" + row.normal_cram),
               file(path_cram + "/" + row.normal_cram+file_ext),
               row.normal_id,
               row.tumor_id]}
      .ifEmpty{exit 1, "${tn_file} was empty - no tumor/normal supplied" }
	//we return the channel
  return tn_pairs
}

// print the calling parameter to the log and a log file
def print_params () {
  //software versions for v2.0
  def software_versions = ['xhla' : '0.0.0',
                           'vep' : '99.2',
                          'pVactools'   : '2.0.2']
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
