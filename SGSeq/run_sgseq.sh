#!/bin/bash 
memPerCore=2.8G # not recommended by Tristan


until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
        --support)
            shift
            support=$1;;
        --pipelineBaseDir)
            shift
            pipelineBaseDir=$1;;
        --code) 
            shift 
            code=$1;;
        --submit)
            shift
            submit=$1;; 
        --outputDir) 
            shift 
            outputDir=$1;;
        --step) 
            shift 
            step=$1;;
        --species) 
            shift 
            species=$1;;
        --submit)
            shift
            submit=$1;;
        -* )
            stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done


step0=${pipelineBaseDir}/sgseq_step0.R
step1a=${pipelineBaseDir}/sgseq_step1a.R
step1b=${pipelineBaseDir}/sgseq_step1b.R
step2=${pipelineBaseDir}/sgseq_step2.R

for script in $step0 $step1a $step1b $step2 $support; do
    if [ ! -e $script ]; then
        stop "$script not found"
    fi
done

script_step0=${outputDir}/cluster/submission/sgseq_step0.sh
script_step1a=${outputDir}/cluster/submission/sgseq_step1a.sh 
script_step1b=${outputDir}/cluster/submission/sgseq_step1b.sh 
script_step2a=${outputDir}/cluster/submission/sgseq_step2a.sh 
script_step2b=${outputDir}/cluster/submission/sgseq_step2b.sh 
Rscript=/share/apps/R-3.3.2/bin/Rscript

Rdir=${outputDir}/cluster/R
out_step0=${Rdir}/sgseq_step0.out
out_step1a=${Rdir}/sgseq_step1a.out
out_step1b=${Rdir}/sgseq_step1b.out
out_step2a=${Rdir}/sgseq_step2a.out
out_step2b=${Rdir}/sgseq_step2b.out
# make directories
for dir in $outputDir ${outputDir}/cluster/out ${outputDir}/cluster/error ${outputDir}/cluster/submission ${Rdir}; do
	if [ ! -e $dir ];then
		mkdir -p $dir
	fi
done

# find species-specific annotations
refFolder=/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/
if [ ! -e $refFolder ];then
    stop "reference folder $refFolder is missing"
fi

case "$species" in  
	worm)
	   gtf=${refFolder}/Worm/Caenorhabditis_elegans.WBcel235.89.gtf
	   annotation=${refFolder}/Worm/biomart_annotations_worm.tab 
     sgseqAnno=${refFolder}/Worm/Caenorhabditis_elegans.WBcel235.89.sgseqAnno.Rdata
	;;
	fly)
     gtf=${refFolder}/Fly/Drosophila_melanogaster.BDGP6.89.gtf
     annotation=${refFolder}/Fly/biomart_annotations_fly.tab
     sgseqAnno=${refFolder}/Fly/Drosophila_melanogaster.BDGP6.89.sgseqAnno.Rdata
  ;;
	mouse)
	   gtf=${refFolder}/Mouse/Mus_musculus.GRCm38.82_fixed.gtf
	   annotation=${refFolder}/Mouse/biomart_annotations_mouse.tab
     sgseqAnno=${refFolder}/Mouse/Mus_musculus.GRCm38.82_fixed.sgseqAnno.Rdata
	;;
	human)
	   gtf=${refFolder}/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gtf
	   annotation=${refFolder}/Human_hg38/biomart_annotations_human.tab
     sgseqAnno=${refFolder}/Human_hg38/Homo_sapiens.GRCh38.82_fixed.sgseqAnno.Rdata 
	;;
	*)
        stop "unknown species $species"
esac

if [ ! -e $gtf ];then
    stop "GTF file $gtf is missing"
else
    echo "GTF file exists"
fi

if [ ! -e $annotation ];then
    stop "annotation file $annotation is missing"
else
  echo "annotation file exists"
fi

function step0 {
    echo "
#$ -S /bin/bash
#$ -l h_vmem=4G,tmem=4G
#$ -l h_rt=72:00:00
#$ -pe smp 1 
#$ -N SGSeq_${code}_step0 
#$ -R y
#$ -o ${outputDir}/cluster/out
#$ -e ${outputDir}/cluster/error
#$ -cwd 

$Rscript $step0 --gtf $gtf --sgseq.anno $sgseqAnno > $out_step0
" > $script_step0
echo "step 0 - create SGSeq transcript objects"
echo $script_step0
}

function step1a {
echo "  
#$ -S /bin/bash
#$ -l h_vmem=${memPerCore},tmem=${memPerCore}
#$ -l h_rt=72:00:00
#$ -pe smp 8  
#$ -R y
#$ -N SGSeq_${code}_step1a
#$ -o ${outputDir}/cluster/out
#$ -e ${outputDir}/cluster/error
#$ -cwd 
}

$Rscript --vanilla ${step1a} --support.tab ${support} \
  --code ${code} --output.dir ${outputDir} --gtf ${gtf} \
  --sgseq.anno ${sgseqAnno} > ${out_step1a}
" > $script_step1a

    echo "step1a - find all annotated splicing events"
    echo $script_step1a

}

function step1b {
#if [ $step = "step1b" ]; then  
    echo "  
#$ -S /bin/bash
#$ -l h_vmem=${memPerCore},tmem=${memPerCore}
#$ -l h_rt=72:00:00
#$ -pe smp 8
#$ -N SGSeq_${code}_step1b  
#$ -R y
#$ -o ${outputDir}/cluster/out
#$ -e ${outputDir}/cluster/error
#$ -cwd 

$Rscript --vanilla ${step1b} --support.tab ${support} \
  --code ${code} --output.dir ${outputDir} --gtf ${gtf} \
  --sgseq.anno ${sgseqAnno}  > ${out_step1b}
" > $script_step1b
    echo "step1b - find all annotated AND novel splicing events"
    echo $script_step1b
}

function step2a {
    echo "  
#$ -S /bin/bash
#$ -l h_vmem=${memPerCore},tmem=${memPerCore}
#$ -l h_rt=72:00:00
#$ -pe smp 1  
#$ -N SGSeq_${code}_step2a
#$ -R y
#$ -o ${outputDir}/cluster/out
#$ -e ${outputDir}/cluster/error
#$ -cwd 

$Rscript --vanilla ${step2} --step step2a --support.tab ${support} \
  --code ${code} --output.dir ${outputDir} --annotation ${annotation} > ${out_step2a}
" > $script_step2a
    echo "step2a - run DEXSeq using the annotated event counts from SGSeq"
    echo $script_step2a
}

#if [ $step = "step2b" ]; then  

function step2b {
    echo "  
#$ -S /bin/bash
#$ -l h_vmem=${memPerCore},tmem=${memPerCore}
#$ -l h_rt=72:00:00
#$ -pe smp 1
#$ -N SGSeq_${code}_step2b  
#$ -R y
#$ -o ${outputDir}/cluster/out
#$ -e ${outputDir}/cluster/error
#$ -cwd 

$Rscript --vanilla ${step2} --step step2b --support.tab ${support} \
  --code ${code} --output.dir ${outputDir} --annotation ${annotation} > ${out_step2b}
" > $script_step2b
    echo "step2b - run DEXSeq using the annotated AND novel event counts from SGSeq"
    echo $script_step2b
}

# make SGSeq_${code} transcript annotation data if doesn't exist already
if [ ! -e $sgseqAnno ];then
    step0
    if [[ "$submit" == "yes" ]];then
        qsub $script_step0
        hold="-hold_jid SGSeq_${code}_step0"
    fi
fi


if [[ "$step" == "step1a" ]]; then  
    step1a
    step2a
    if [[ "$submit" == "yes" ]];then
        qsub $hold $script_step1a
        qsub -hold_jid SGSeq_${code}_step1a $script_step2a
    fi
fi

if [ $step == "step1b" ]; then
    step1b
    step2b
    if [[ "$submit" == "yes" ]];then
        qsub $hold $script_step1b
        qsub -hold_jid SGSeq_${code}_step1b $script_step2b
    fi
fi

