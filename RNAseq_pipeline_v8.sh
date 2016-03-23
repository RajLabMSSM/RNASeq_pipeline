#! /bin/bash
# this script has only been tested with bash and may not work with other shells

# script exits if return value of a command is not zero
set -e
# this forces all variables to be defined
set -u
# for debugging prints out every line before executing it
set -x

# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { error "$*"; exit 1; }
try() { "$@" || stop "cannot $*"; }
function files_exist() { for file in $*; do if [ ! -e $file ]; then stop "File $file does not exist" ; fi; done }

#computer=vanHeel
computer=
superLong=

if [[ "$computer" == "vanHeel" ]]
then
    software=/data_n2/vplagnol/Software
    pythonbin=/software/additional/epd-7.3.1/bin/python
    Rbin=/data_n2/vplagnol/Software/R-3.0.2/bin/R
    Rscript=/data_n2/vplagnol/Software/R-3.0.2/bin/Rscript
    dexseqCount=/data_n2/vplagnol/Rlibs/installed/DEXSeq/python_scripts/dexseq_count.py
    javaTemp2="/data_n1/vanheel_singlecellgenomics/tmp"
    javaTemp="TMP_DIR=${javaTemp2}"
else
    software=/cluster/project8/vyp/vincent/Software
    pythonbin=/share/apps/python-2.7.1/bin/python2.7
    if [ ! -e $pythonbin ]; then pythonbin=/share/apps/python-2.7.8/bin/python2.7; fi
    ##Rbin=/cluster/project8/vyp/vincent/Software/R-3.1.2/bin/R
    misoRunEvents=/cluster/project8/vyp/vincent/Software/misopy-0.4.9/misopy/run_events_analysis.py
    runMiso=/cluster/project8/vyp/vincent/Software/misopy-0.4.9/misopy/run_miso.py
    javaTemp2="/scratch2/vyp-scratch2/vincent/java_temp"
    if [ ! -e $javaTemp2 ]; then javaTemp2="/cluster/scratch3/vyp-scratch2/vincent/java_temp/"; fi
    javaTemp="TMP_DIR=${javaTemp2}"
    java=/share/apps/jdk1.7.0_45/bin/java
    if [ ! -e $java ]; then java=/share/apps/jdk1.8.0_25/bin/java; fi
    dexseqCount="/cluster/project8/vyp/vincent/libraries/R/installed/DEXSeq/python_scripts/dexseq_count.py"
    bigFilesBundleFolder=/scratch2/vyp-scratch2/reference_datasets
    if [ ! -e $bigFilesBundleFolder ]
    then
        bigFilesBundleFolder=/cluster/scratch3/vyp-scratch2/reference_datasets
    fi
    #which R should we be using?
    #optparse available but just run this to install:
    # > install.packages('getopt')
    # wget --no-check-certificate https://cran.r-project.org/src/contrib/optparse_1.3.2.tar.gz
    # $Rbin CMD install optparse_1.3.2.tar.gz
    Rbin=/cluster/project8/vyp/vincent/Software/R-3.2.2/bin/R
    Rscript=/cluster/project8/vyp/vincent/Software/R-3.2.2/bin/Rscript
    # data.table cannot be installed with this version of R
    #Rbin=/share/apps/R/bin/R
    #Rscript=/share/apps/R/bin/Rscript
fi


#RNASEQPIPBASE=/Users/pontikos/bin/RNASeq_pipeline/
#RNASEQPIPBASE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline

echo "Base of RNA-Seq pipeline is located here: $RNASEQPIPBASE"
countPrepareR=${RNASEQPIPBASE}/counts_prepare_pipeline.R
files_exist $countPrepareR
dexseqFinalProcessR=${RNASEQPIPBASE}/dexseq_pipeline_v2.R
deseqFinalProcessR=${RNASEQPIPBASE}/deseq2_pipeline.R
pathwayGOAnalysisR=${RNASEQPIPBASE}/pathwayGO_pipeline.R
topGOAnalysisR=${RNASEQPIPBASE}/topGO_pipeline.R
novosort=${software}/novocraft3/novosort
trim_galore=${software}/trim_galore/trim_galore
cutadapt=/share/apps/python-2.7.8/bin/cutadapt
#for the old cluster
if [ ! -e $cutadapt ];then cutadapt=/share/apps/python-2.7.6/bin/cutadapt;fi

starexec=/cluster/project8/vyp/vincent/Software/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR
samtools=${software}/samtools-1.2/samtools
rseqQCscripts=${software}/RSeQC-2.3.3/scripts

picardDup=${software}/picard-tools-1.100/MarkDuplicates.jar
picardStats=${software}/picard-tools-1.100/BamIndexStats.jar
picardMetrics=${software}/picard-tools-1.100/CalculateHsMetrics.jar
picardReorder=${software}/picard-tools-1.100/ReorderSam.jar

## should sex chromosomes be kept in the differential expression analysis?
## Niko: if sex matches then yes?
keepSex=FALSE  
force=yes
species=mouse
segmentLength=25

mart=ensembl
db=mmusculus_gene_ensembl
summary=no
prepareCounts=no
Rdexseq=no
Rdeseq=no
RpathwayGO=no
RtopGO=no
oFolder=temp ##default output

#mkdir -p creates folders only if they do not already exist
mkdir -p $oFolder

submit=no
# If QC is wanted, the pipeline will send each fastq or pair of fastqs through FastQC. If adapters are present then the offending fastq files will be trimmed with Trim Galore! and these trimmmed fastqs will be the ones aligned with STAR. 
QC=no
starStep1a=no
starStep1b=no
starStep2=no

stranded=no
libstrand=fr-unstranded
keepDups=FALSE

code=""
dataframe=none
iFolder=""
misoindex=NA
oFolder=RNAseq_processed
stem=""

until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--fastqFiles)
	    shift
	    i=0
	    for fileloc in $@; do 
		fastqFiles[ $i ]=$fileloc
		((i=i+1))
	    done;;
	--stranded)
	    shift
	    stranded=yes
	    libstrand=$1;;
	--species)
	    shift
	    species=$1;;
	--force)
	    shift
	    force=$1;;
	--starStep1a)
	    shift
	    starStep1a=$1;;
	--starStep1b)
	    shift
	    starStep1b=$1;;
	--starStep2)
	    shift
	    starStep2=$1;;
	--keep.sex)
	    shift
	    keepSex=$1;;
	--misoindex)
	    shift
	    misoindex=$1;;
	--summary)
	    shift
	    summary=$1;;
	--miso)
	    shift
	    miso=$1;;
	--prepareCounts)
	    shift
	    prepareCounts=$1;;
	--Rdexseq)
	    shift
	    Rdexseq=$1;;
	--Rdeseq)
	    shift
	    Rdeseq=$1;;
	--RpathwayGO)
	    shift
	    RpathwayGO=$1;;
	--RtopGO)
	    shift
	    RtopGO=$1;;
	--iFolder)
	    shift
	    iFolder=$1;;
	--oFolder)
        shift
        oFolder=$1;;
	--dataframe)
        shift
        dataframe=$1;;
	--code)
        shift
        code=$1;;
	--submit)
	    shift
        submit=$1;;
    --stem)
        shift
        stem=$1;;
	--keepDups)
	    shift
	    keepDups=TRUE;;
	--QC)
	    shift
	    QC=$1;;
	-* )
	    stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done


echo "Strand information $stranded $libstrand"

########################## estimate the duration of the jobs

countStrand=no
if [[ "$libstrand" == "fr-firststrand" ]]
then
  countStrand=yes
  countStrandReverse=reverse
elif [[ "$libstrand" == "fr-secondstrand" ]]
then
  countStrand=reverse
  countStrandReverse=yes
else
    echo unknown libstrand $libstrand
fi

if [[ "$stem" == "" ]]
then
    stem=$code
fi

if [[ "$superLong" == "yes" ]]
then
    ((nhours=nhours+nhours))
fi

## create the output folders
clusterFolder=${oFolder}/cluster

mkdir -p ${oFolder} ${clusterFolder} ${clusterFolder}/out ${clusterFolder}/error ${clusterFolder}/R  ${clusterFolder}/submission

files_exist $dataframe 

cp "$dataframe" "${oFolder}/"


case "$species" in
    zebrafish)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Danio_rerio/NCBI/Zv9
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf    
        fasta=${refFolder}/Sequence/WholeGenomeFasta/genome.fa
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
        gffFile=${RNASEQBUNDLE}/zebrafish/GTF/zebrafish_iGenomes_Zv9_with_ensembl.gff
        annotationFile=${RNASEQBUNDLE}/zebrafish/biomart/biomart_annotations_zebrafish.tab
        #geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/zebraFish_refSeqTable_zv9_nochr.bed
        #geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/zebraFish_refSeqTable_zv9_chr1.bed
        ;;
    DvH_sc_human)
        refFolder=/data_n1/vanheel_singlecellgenomics/support/Homo_sapiens/NCBI/build37.2
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf    
        fasta=${refFolder}/Sequence/WholeGenomeFasta/genome.fa
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/human_b37_with_spikes
        ;;
    human_hg38)
        fasta=${bigFilesBundleFolder}/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
        gtfFile=${bigFilesBundleFolder}/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gtf
        gffFile=${bigFilesBundleFolder}/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gff
        STARdir=${bigFilesBundleFolder}/RNASeq/Human_hg38/STAR 
	annotationFile=${bigFilesBundleFolder}/RNASeq/Human_hg38/biomart_annotations_human.tab
        #annotationFile=/scratch2/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab
	;;
    humanmuscle)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Homo_sapiens/NCBI/build37.2
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf	
        gffFile=/cluster/project8/vyp/vincent/data/reference_genomes/gff/humanmuscle_iGenomes_NCBI37_with_ensembl.gff
        #annotationFile=${bigFilesBundleFolder}/human/biomart/biomart_annotations_human.tab
        annotationFile=/scratch2/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab
        geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/homoSapiens_geneTable_hg19_nochr.bed
        geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/homoSapiens_geneTable_hg19_chr1.bed
        db=hsapiens_gene_ensembl
        ;;
    Dict_Disc_masked)
        IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Dict/dicty_masked_ERCC92
        gtfFile=${bigFilesBundleFolder}/RNASeq/Dict/dict_no_spike.gtf
        gffFile=MISSING
        annotationFile=not_done_yet
        ;;
    Dict_Disc)
        IndexBowtie2=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/sequence/Dictyostelium_discoideum.dictybase.01.23.dna.genome
        gtfFile=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/GTF/Dictyostelium_discoideum.dictybase.01.23.gtf
        gffFile=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/GTF/Dictyostelium_discoideum.dictybase.01.23.gff
        annotationFile=not_done_yet
        ;;
    pig)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Sus_scrofa/NCBI/Sscrofa10.2
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf
        gffFile=${RNASEQBUNDLE}/pig/GTF/pig_iGenomes_NCBI_10_2_with_ensembl.gff
        annotationFile=${RNASEQBUNDLE}/pig/biomart/biomart_annotations_pig.tab
        ;;
    chicken)
        IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Chicken/Gallus_gallus.Galgal4.dna.toplevel
        gtfFile=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/chicken/GTF/Gallus_gallus.Galgal4.78.gtf
        gffFile=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/chicken/GTF/Gallus_gallus.Galgal4.78.gff
        annotationFile=${RNASEQBUNDLE}/chicken/biomart/biomart_annotations_chicken.tab
        ;;
    rat)
        IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Rat/Rattus_norvegicus.Rnor_5.0.dna_rm.toplevel
        gtfFile=${bigFilesBundleFolder}/RNASeq/Rat/Rattus_norvegicus.Rnor_5.0.79.gtf
        gffFile=${bigFilesBundleFolder}/RNASeq/Rat/Rattus_norvegicus.Rnor_5.0.79.gff
        annotationFile=${RNASEQBUNDLE}/rat/biomart/biomart_annotations_rat.tab
        ;;
    sheep)
        IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Sheep/Ovis_aries.Oar_v3.1.dna_rm.toplevel
        gtfFile=${bigFilesBundleFolder}/RNASeq/Sheep/Ovis_aries.Oar_v3.1.80.gtf
        gffFile=${bigFilesBundleFolder}/RNASeq/Sheep/Ovis_aries.Oar_v3.1.80.gff
        annotationFile=${RNASEQBUNDLE}/sheep/biomart/biomart_annotations_sheep.tab
        ;;
    drosophila)
        refFolder=${bigFilesBundleFolder}/RNASeq/Drosophila/Drosophila_melanogaster/NCBI/build5.41
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome	
        gtfFile=${RNASEQBUNDLE}/drosophila/GTF/Drosophila_melanogaster.BDGP5.75.gtf
        gffFile=${RNASEQBUNDLE}/drosophila/GTF/Drosophila_melanogaster.BDGP5.75.gff
        annotationFile=${RNASEQBUNDLE}/drosophila/biomart/biomart_annotations_drosophila.tab
        ;;
    dog)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Canis_familiaris/NCBI/build3.1
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome	
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf
        gffFile=${RNASEQBUNDLE}/dog/GTF/dog_iGenomes_NCBI_3_1_with_ensembl.gff
        annotationFile=${RNASEQBUNDLE}/dog/biomart/biomart_annotations_dog.tab
        ;;
    mouse)
        STARdir=${bigFilesBundleFolder}/RNASeq/Mouse/STAR
        annotationFile=${bigFilesBundleFolder}/RNASeq/Mouse/biomart_annotations_mouse.tab
        gffFile=${bigFilesBundleFolder}/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gff
        gtfFile=${bigFilesBundleFolder}/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf
        ;;
    tc1_mouse)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/Zanda_AD_Tc1J20_RNASeq/Zanda_Tc1_reference/build1 
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
        gtfFile=${RNASEQBUNDLE}/Tc1_mouse/GTF/Tc1.gtf #${refFolder}/Annotation/Genes/genes.gtf
        gffFile=${RNASEQBUNDLE}/Tc1_mouse/GTF/Tc1.gff  #${refFolder}/gff/tc1.gff
        cleanGtfFile=${RNASEQBUNDLE}/Tc1_mouse/GTF/Tc1.gtf #${refFolder}/Annotation/Genes/genes.gtf
        annotationFile=${RNASEQBUNDLE}/Tc1_mouse/tc1_annotations.tab #${refFolder}/annotations/biomart/tc1.tab
        geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/mm9_NCBI37_Ensembl_chr1.bed
        geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/mm9_NCBI37_Ensembl_nochr.bed
        ##db=mmusculus_gene_ensembl
        ;;
    *)
        stop "unknown species $species"
esac

files_exist $gtfFile $annotationFile $gffFile $STARdir

############### checking the input dataframe
h1=`awk '{if (NR == 1) print $1}'  $dataframe` 
h2=`awk '{if (NR == 1) print $2}'  $dataframe` 
h3=`awk '{if (NR == 1) print $3}'  $dataframe` 
if [[ "$h1" != "sample" ]]; then echo "header 1 must be sample"; exit; fi
if [[ "$h2" != "f1" ]]; then echo "header 2 must be f1 for fastq1"; exit; fi
if [[ "$h3" != "f2" ]]; then echo "header 3 must be f2 for fastq2"; exit; fi



hold=""
SCRATCH_DIR=/scratch0/RNASeq_${code}
JAVA_DIR=${SCRATCH_DIR}/javastar
 
# alignment
function starSubmissionStep1a {
    starSubmissionStep1a=${oFolder}/cluster/submission/starSubmissionStep1a.sh
    echo "
#$ -S /bin/bash
#$ -l h_vmem=15G,tmem=15G
#$ -l h_rt=72:00:00
#$ -pe smp 4
#$ -R y
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -N step1a_${code}
#$ -wd ${oFolder}
echo \$HOSTNAME >&2
date >&2
mkdir -p $JAVA_DIR
" > $starSubmissionStep1a
    tail -n +2  $dataframe | while read sample f1 f2 condition
    do
        if [[ "$f2" == "NA" ]]; then paired=no;  else paired=yes; fi;
        echo "Sample $sample"
        finalOFolder=${oFolder}/${sample}
        dexseqfolder=${oFolder}/${sample}/dexseq
	mkdir -p ${finalOFolder} ${dexseqfolder}
        # go no further
        if [[  -e ${finalOFolder}/${sample}_unique.bam.bai ]]
        then
            echo ${finalOFolder}/${sample}_unique.bam.bai exists will not align again
            return 0
        fi
        files_exist ${iFolder}/$f1
        echo Aligning with STAR
        if [[ $paired == "yes" ]]
        then
            files_exist ${iFolder}/$f2
            #QC paireend
            if [[ "$QC" == "yes" ]];then
# use the currently redundant flag $summary in case the fastqs are trimmed but the alignment had failed etc.
	    	echo $summary		
		if [[ ! "$summary" == "trimmed_exist" ]];then
                	echo "
$trim_galore --gzip -o $iFolder --path_to_cutadapt $cutadapt --paired ${iFolder}/$f1 ${iFolder}/$f2
"
		 >>  $starSubmissionStep1a
		fi

			#the trimmed files have a slightly different output
## For some stupid reason Trim_Galore in paired end mode appends file names differently than in single end mode. Who'd have thought? 		
# some fastqs have the naming scheme sample.a1.fastq.gz so I need to specifically split off the .fastq.gz.
		if [[ "$f2" == "NA" ]]; then
			f1=`echo $f1 | sed 's/.fastq.gz/_trimmed.fq.gz/g'`
		elif [[ ! "$f2" == "NA" ]]; then 
			f1=`echo $f1 | sed 's/.fastq.gz/_val_1.fq.gz/g'`
			f2=`echo $f2 | sed 's/.fastq.gz/_val_2.fq.gz/g'`
		fi
                #check that the trimming has happened. If not then exit
                echo "
if [ ! -e ${iFolder}/$f1 ]; then exit;fi
                " >> $starSubmissionStep1a
            fi
            #if QC step is wanted and ran successfully then the trimmed fastqs should be aligned.
            echo "
${starexec} --readFilesIn ${iFolder}/$f1 ${iFolder}/$f2 --readFilesCommand zcat --genomeLoad LoadAndKeep --genomeDir ${STARdir} --runThreadN  4 --outFileNamePrefix ${SCRATCH_DIR}/${sample} --outSAMtype BAM Unsorted --outSAMunmapped Within --outSAMheaderHD ID:${sample} PL:Illumina
date >&2
# sort reads and mark duplicates all in one step
$novosort --md --xs -f -t /scratch0/ -6 -c 4 -m 50G ${SCRATCH_DIR}/${sample}Aligned.out.bam -o ${finalOFolder}/${sample}_unique.bam
date >&2
mv ${SCRATCH_DIR}/${sample}Log* ${finalOFolder}/
rm ${SCRATCH_DIR}/${sample}Aligned.out.bam

" >> $starSubmissionStep1a
        #if single ended
        else
            if [[ "$QC" == "yes" ]];then
		echo $summary           
                if [[ ! "$summary" == "trimmed_exist" ]];then
                echo "
$trim_galore --gzip -o $iFolder --path_to_cutadapt $cutadapt ${iFolder}/$f1
" >> $starSubmissionStep1a
                fi
		#the trimmed files have a slightly different output
                f1=`echo $f1 | sed 's/.fastq.gz/_trimmed.fq.gz/g'`
                #check that the trimming has happened. If not then exit
                echo "  if [ ! -e ${iFolder}/$f1 ]; then exit;fi " >> $starSubmissionStep1a
            fi
            #if trimming has occurred then the trimmed fastq will be aligned
        # STAR    
	echo "
${starexec} --readFilesIn ${iFolder}/$f1 --readFilesCommand zcat --genomeLoad LoadAndKeep --genomeDir ${STARdir} --runThreadN  4 --outFileNamePrefix ${SCRATCH_DIR}/${sample} --outSAMtype BAM Unsorted --outSAMunmapped Within --outSAMheaderHD ID:${sample} PL:Illumina
date >&2
# sort reads and mark duplicates
$novosort --md --xs -f -t /scratch0/ -6 -c 4 -m 50G  ${SCRATCH_DIR}/${sample}Aligned.out.bam -o ${finalOFolder}/${sample}_unique.bam
date >&2
mv ${SCRATCH_DIR}/${sample}Log* ${finalOFolder}/
rm ${SCRATCH_DIR}/${sample}Aligned.out.bam

" >> $starSubmissionStep1a
        fi
    done
    echo "
rm -rf $JAVA_DIR
" >> $starSubmissionStep1a
    echo "
${starexec} --genomeLoad Remove --genomeDir ${STARdir}

rm -rf ${SCRATCH_DIR}

" >> $starSubmissionStep1a
}


# sorting and duplication removal
function starSubmissionStep1b {
# per sample
# if the master table has been made before then remove it. Otherwise every time the submission script is run then more lines are appended to it.
  starMasterTableStep1b=${oFolder}/cluster/submission/starMasterTableStep1b.tab
  if [ -e $starMasterTableStep1b ]; then rm $starMasterTableStep1b;fi
  tail -n +2  $dataframe | while read sample f1 f2 condition
  do
	finalOFolder=${oFolder}/${sample}
   echo "
${samtools} index ${finalOFolder}/${sample}_unique.bam
${samtools} flagstat ${finalOFolder}/${sample}_unique.bam > ${finalOFolder}/${sample}_mappingStats.tab
" > ${oFolder}/cluster/submission/star_step1b_${sample}.sh
      echo "${oFolder}/cluster/submission/star_step1b_${sample}.sh" >> $starMasterTableStep1b
  done
  njobs1b=`wc -l $starMasterTableStep1b | awk '{print $1}'`
# submission script
  starSubmissionStep1b=${oFolder}/cluster/submission/starSubmissionStep1b.sh
  echo "
#$ -S /bin/bash
#$ -l h_vmem=9.5G
#$ -l tmem=9.5G
#$ -l h_rt=24:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -N step1b_${code}
#$ -wd ${oFolder}
#$ -t 2-${njobs1b}
#$ -tc 20
echo \$HOSTNAME >&2
script=\`awk '{if (NR == '\$SGE_TASK_ID') print}' $starMasterTableStep1b\`
sh \$script
" > $starSubmissionStep1b
}

# dexseqCount
function starSubmissionStep2 {
# per sample
    starMasterTableStep2=${oFolder}/cluster/submission/starMasterTableStep2.tab
    echo "scripts" > $starMasterTableStep2

    #echo $PWD/${dataframe}; exit
    #echo ${oFolder}/${dataframe}; exit
    
    ##echo $dataframe; exit
    tail -n +2  $dataframe | while read sample f1 f2 condition; do
        if [[ "$f2" == "NA" ]]; then paired=no;  else paired=yes; fi;	
        dexseqfolder=${oFolder}/${sample}/dexseq

	finalOFolder=${oFolder}/${sample}
	echo "
$samtools view -F 0x0400 ${finalOFolder}/${sample}_unique.bam |  ${pythonbin} ${dexseqCount} --order=pos --paired=${paired} --stranded=${countStrand}  ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts.txt
$samtools view ${finalOFolder}/${sample}_unique.bam |  ${pythonbin} ${dexseqCount} --order=pos --paired=${paired} --stranded=${countStrand}  ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts_keep_dups.txt
" > ${oFolder}/cluster/submission/star_step2_${sample}.sh

        echo "${oFolder}/cluster/submission/star_step2_${sample}.sh" >> $starMasterTableStep2

    done
    #((nscripts=nscripts+1))
    njobs2=`wc -l $starMasterTableStep2 | awk '{print $1}'`
# submission script
    starSubmissionStep2=${oFolder}/cluster/submission/starSubmissionStep2.sh
    echo "
#$ -S /bin/bash
#$ -l h_vmem=5.9G
#$ -l tmem=5.9G
#$ -l h_rt=12:00:00
#$ -N step2_${code}
#$ -R y
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -wd ${oFolder}
#$ -t 2-${njobs2}
#$ -tc 20
echo \$HOSTNAME >&2
script=\`awk '{if (NR == '\$SGE_TASK_ID') print}' $starMasterTableStep2\`
sh \$script
" > $starSubmissionStep2
}




################################################# Now the scripts that take all samples together    

function starSubmissionStep3 {
# analyse all samples together
    starSubmissionStep3=${oFolder}/cluster/submission/starSubmissionStep3.sh    
    ncores=1
    nhours=0
    nminutes=0
    mem=0
    if [[ "$prepareCounts" == "yes" ]]; then mem=13.9; ((nminutes=nminutes+50)); fi
    if [[ "$Rdeseq" == "yes" ]]; then ((nhours=nhours+3)); mem=13.9; fi
#scale with the number of samples
    if [[ "$Rdexseq" == "yes" ]]; then ((nhours=nhours+18)); ncores=4;mem=5.9; fi
    #if [[ "$RpathwayGO" == "yes" ]]; then ((nhours=nhours+3)); mem=6; fi
    #if [[ "$RtopGO" == "yes" ]]; then ((nhours=nhours+3)); mem=6; fi
    echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -cwd
#$ -pe smp $ncores
#$ -l tmem=${mem}G,h_vmem=${mem}G
#$ -V
#$ -R y
#$ -l h_rt=${nhours}:${nminutes}:00
echo \$HOSTNAME >&2
" > $starSubmissionStep3
    if [[ "$prepareCounts" == "yes" ]]
    then
        files_exist ${annotationFile} ${dataframe} ${oFolder} ${countPrepareR}
    echo "
${Rscript} ${countPrepareR} --gff ${gffFile} --annotation.file ${annotationFile} --keep.dups ${keepDups} --support.frame ${dataframe} --code ${code} --iFolder ${oFolder} > ${clusterFolder}/R/count_prepare.out
    " >> $starSubmissionStep3
    fi
    ##############
    if [[ "$Rdeseq" == "yes" ]]
    then
        files_exist $deseqFinalProcessR $dataframe
    echo "
${Rscript} ${deseqFinalProcessR} --keep.sex ${keepSex} --support.frame ${dataframe} --keep.dups ${keepDups} --code ${code} --annotation.file ${annotationFile} --iFolder ${oFolder} > ${clusterFolder}/R/deseq_${stem}.out 
    " >> $starSubmissionStep3
    fi
    ##############
    if [[ "$Rdexseq" == "yes" ]]
    then
        files_exist $dexseqFinalProcessR $dataframe
    echo "
${Rscript} ${dexseqFinalProcessR} --gff ${gffFile} --keep.sex ${keepSex} --keep.dups ${keepDups} --support.frame ${dataframe} --code ${code} --annotation.file ${annotationFile} --iFolder ${oFolder} > ${clusterFolder}/R/dexseq_${stem}.out
    " >> $starSubmissionStep3
    fi
    ##############
    if [[ "$RpathwayGO" == "yes" ]]
    then
        files_exist $pathwayGOAnalysisR $dataframe
    echo "
${Rscript} ${pathwayGOAnalysisR} --support.frame ${dataframe} --code ${code} --mart ${mart} --db ${db} --iFolder ${oFolder} > ${clusterFolder}/R/pathwayGO_${stem}.out 
    " >> $starSubmissionStep3
    fi
    ##############
    if [[ "$RtopGO" == "yes" ]]
    then
        files_exist $topGOAnalysisR $dataframe
    echo "
${Rscript} ${topGOAnalysisR} --support.frame ${dataframe} --code ${code} --mart ${mart} --db ${db} --iFolder ${oFolder} > ${clusterFolder}/R/topGO_${stem}.out 
    " >> $starSubmissionStep3
    fi
    #############
}




if [[ "$starStep1a" == "yes" ]]
then
   echo step1a: align
   starSubmissionStep1a
   files_exist $starSubmissionStep1a 
   if [[ "$submit" == "yes" ]]
   then
       qsub $hold $starSubmissionStep1a
       if [[ "$hold" == "" ]]; then hold="-hold_jid step1a_${code}"; else hold="$hold,step1b_${code}"; fi
   fi
fi

if [[ "$starStep1b" == "yes" ]]
then
   echo step1b: sorting and duplication removal
   starSubmissionStep1b
   files_exist $starSubmissionStep1b
   if [[ "$submit" == "yes" ]]
   then
       qsub $hold $starSubmissionStep1b
       if [[ "$hold" == "" ]]; then hold="-hold_jid step1b_${code}"; else hold="$hold,step1b_${code}"; fi
   fi
fi

if [[ "$starStep2" == "yes" ]]
then
   echo step2: dexseq count
   starSubmissionStep2
   files_exist $starSubmissionStep2
   if [[ "$submit" == "yes" ]]
   then
       qsub $hold $starSubmissionStep2
       if [[ "$hold" == "" ]]; then hold="-hold_jid step2_${code}"; else hold="$hold,step2_${code}"; fi
   fi
fi

if [[ "$prepareCounts" == "yes" || "$Rdeseq" == "yes" || "$Rdexseq" == "yes" || "$RpathwayGO" == "yes" || "$RtopGO" == "yes" ]]
then
    echo "Taking all samples together"
    starSubmissionStep3
    ls -ltrh $starSubmissionStep3
    if [[ "$submit" == "yes" ]]
    then
        qsub $hold $starSubmissionStep3
    fi
fi


