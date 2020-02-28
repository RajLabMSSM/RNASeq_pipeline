import pandas as pd
import os

# SGSeq-pipeline
shell.prefix("ml samtools; ml R/3.6.0") # change if you're on a different HPC

code = config["code"]
outputDir = config["outputDir"]
support = config["support"]

gtf = "" # GENCODE V30
annotation = "biomart_annotations_human.tab" 

# build SGSeq reference
rule step0:
    input:
        gtf = gtf
    output:
        "indexes/GENCODE_v30_hg38.sgseqAnno.Rdata"
    params:
        script = "scripts/sgseq_step0.R"
    shell:
        "Rscript $step0 --gtf {input.gtf} --sgseq.anno {output}"


metadata = pd.read_csv(support, sep = "\t")

bam_files = metadata["file_bam"]

# get all condition columns in metadata

# for each one generate a condition string
condition_col_names = [cond for cond in list(metadata.columns) if "condition" in cond]

condition_strings = []
for cond in condition_col_names:
    cond_values = metadata[cond]
    cond_values = [x for x in cond_values if not pd.isna(x) ]
    cond_string = "_".join(sorted(list(set(cond_values))))
    condition_strings.append(cond_string)

region = "chr??" # FUS gene

rule all:
    input:
        expand(outputDir + "{comparison}/" + dataCode + "{comparison}_sgseq_variant_type_table.tab", comparison = condition_strings)

# extract region from BAM files and index
rule extractFUS:
    input:
        inFolder + "{sample}.bam"
    output:
        inFolder + "{sample}_FUS.bam",
        inFolder + "{sample}_FUS.bam.bai"
    shell:
        "samtools view -bh {input} {region} > {output}"
        "samtools index {output}" 


# run SGSeq to find novel and annotated isoforms
rule step1b:
	input:
        bams = expand( inFolder + "{sample}_FUS.bam" , sample = samples )
        support,
        sgseqAnno = "indexes/GENCODE_v30_hg38.sgseqAnno.Rdata"
    params:
        script = "scripts/sgseq_step1b.R",
        code = dataCode,
        outputDir = outputDir
    output:
       outputDir + dataCode + "_txf_novel.RData" 
    shell:
        "Rscript --vanilla {params.script} --support.tab {input.support}"
        " --code {params.code} --output.dir {params.outputDir} --gtf ${gtf} "
        " --sgseq.anno ${sgseqAnno}  "

# run DEXSeq to test
rule step2b:
    input:
        support = support
        # output of step1 
    output:
        outputDir + comparisonFolder + "/" + dataCode + comparisonFolder + "_res_clean_novel.tab"
    params:
        script = "scripts/sgseq_step2.R",
        code = dataCode
    shell:
        "Rscript --vanilla {params.script} --step step2b --support.tab {input.support} "
        " --code {params.code} --output.dir {outputDir} --annotation ${annotation}"

# create output tables
rule step3:
    input:
        outputDir + "{comparison}/" + dataCode + "{comparison}_sgseq_variant_type_table.tab"
    output:
        outputDir + "{comparison}/" + dataCode + "{comparison}_sgseq_variant_type_table.tab"
    params:
        createVarTable = "scripts/createVariantTable.R",
        findCentralExons = "scripts/findCentralExons.R"
    shell:
        "Rscript --vanilla {params.createVarTable} --step step2b --support.tab {support} "
        "--code {dataCode} --output.dir {outputDir} ;"
        "Rscript --vanilla $findCentralExons --step step2b --support.tab {support} "
        "--code {dataCode} --output.dir {outputDir}"




