# nf-hla-neo
Pipeline to predict neoantigens  from WGS of T/N pairs


## Usage
  ```
  #using a tn_pairs file
  nextflow run iarcbioinfo/nf-hla-neo -r v1.0  \
  -profile singularity --ref chr6.mhc.fa \
  --tn_file cohort_neoantigen.tsv --cram_dir cram \
  --vcf_dir vcfs --vep_dir vep-db-99 --output_folder results_hla_neo
  
  ```

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
	- [xHLA](https://github.com/humanlongevity/HLA)
	- [VEP](https://github.com/Ensembl/ensembl-vep)
	- [pvactools](https://github.com/griffithlab/pVACtools)
	
You can avoid installing all the external software by only installing Docker or singularity.
See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.


## Input (mandatory)

  | Type      | Description   |
  |-----------|---------------|
  |--tn_file		|         [file] File containing list of T/N bam/cram files to be processed |
  |    --cram_dir|         [dir]  directory where the BAM or CRAM  file are stored |
  |    --vcf_dir |         [dir]  directory where the VCF files are stored |
  |--vep_dir      |   [dir] directory containing VEP database for annotation [hg38, GENCODE 33] |
  |    --ref      |         [file] fasta file of chr6 of reference genome [chr6-hg38.fa], shold be indexed with BWA-mem|


### Example of cohort_neoantigen.tsv file (--tn_file)
A text file tabular separated, with the following header:

```
id      vcf     normal_cram     normal_id       tumor_id
sample1	sample1.vcf.gz	sample1_N.cram	NORMAL1	TUMOR1
sample2	sample2.vcf.gz	sample2_N.cram	sample2_N.cram	NORMAL2	TUMOR2
sample3	sample3.vcf.gz	sample3_N.cram	sample3_N.cram	NORMAL3	TUMOR3_1,TUMOR3_2
``` 

### Optional parameters

| Name      | type | Description     |
|-----------|---------------|-----------------|
| --pvactools_predictors | [string]| predictions tools to compute neoantigens [def:all_class_i,all_class_ii or NetMHCpan,NetMHCIIpan]|
|      --bam     |       [flag] |active bam mode [def:cram]|
|     --output_folder |  [string] |name of output folder |
|      --cpu          |[Integer] | Number of CPUs[def:2] |
|      --mem |        [Integer] | Max memory [def:8Gb] |  



## Output

```
results
└── xHLA                        # HLA-typing predictions
│   ├── report-MESO_001-hla.json
│   ├── report-MESO_002-hla.json
│   ├── report-MESO_003-hla.json
│   ├── report-MESO_004-hla.json
│   ├── report-MESO_005-hla.json
│   ├── report-MESO_006-hla.json
│   ├── report-MESO_007-hla.json
│   ├── report-MESO_008-hla.json
│   ├── report-MESO_009-hla.json
├── VEP							# Annotation of variant impact
│   ├── MESO_001.vep.vcf
│   ├── MESO_002.vep.vcf
│   ├── MESO_003.vep.vcf
│   ├── MESO_004.vep.vcf
│   ├── MESO_005.vep.vcf
│   ├── MESO_006.vep.vcf
|
├── pVACTOOLS									#pvactools predictions
│   ├── MESO_002.pvactools.log			#pvactools log
│   ├── MESO_002_T1_pvactools				
│   │   ├── combined						#neoatigens for class I and class II
│   │   │   ├── B00JALW.all_epitopes.aggregated.tsv    # raw predictions aggregated
│   │   │   ├── B00JALW.all_epitopes.tsv               # raw predictions 
│   │   │   └── B00JALW.filtered.tsv	                  # filtered predictions
│   │   ├── MHC_Class_I						#neoatigens for class I
│   │   │   ├── B00JALW.all_epitopes.aggregated.tsv
│   │   │   ├── B00JALW.all_epitopes.tsv
│   │   │   ├── B00JALW.filtered.tsv
│   │   │   └── log
│   │   │       └── inputs.yml
│   │   └── MHC_Class_II				#neoatigens for class II
│   │       ├── B00JALW.all_epitopes.aggregated.tsv
│   │       ├── B00JALW.all_epitopes.tsv
│   │       ├── B00JALW.fasta
│   │       ├── B00JALW.filtered.tsv
│   │       └── log
|	 				└── inputs.yml
├── nf-pipeline_info # Nextflow information directory
    ├── hla-neo_report.html
    ├── hla-neo_timeline.html
    ├── hla-neo_trace.txt
    └── run_parameters_report.txt # Custom file providing info for software versions and calling parameters
```


## Common errors

### vep-db

The first time is necesary to get a local copy of the vep database, you can achieve this by running the following command within the vep singularity container:

```
#get the singularity container
singularity pull docker://docker.io/iarcbioinfo/ensembl-vep:v1.0
#open a shell and run the folloing command
singularity shell ensembl-vep_v1.0.sif
#get a local copy of the vep cache database (gencode v33)
vep_install -a cf -s homo_sapiens -y GRCh38 -c vep-db-99 --CONVERT

```

### chr6 BWA index
To perform the HLA-typing this pipeline extract reads from the CRAM/BAM file and remap them to the chr6 of GRCh38 using BWA-MEM (0.7.15-r1140), hence the chr6.fa of GRCh38 should be indexed (BWA index) to perform the mapping  of the reads.


### Singularity
The first time that the container is built from the docker image, the TMPDIR  should be defined in a non parallel file-system, you can set this like:

```
export TMPDIR=/tmp
```

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Matthieu Foll*    |            follm@iarc.fr | Developer to contact for support (link to specific gitter chatroom) |
  | Alex Di Genova | digenovaa@fellows.iarc.fr| Developer |



