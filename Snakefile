import pysam
import pandas as pd
import os

# SGSeq-pipeline
shell.prefix("ml samtools; ml R/3.6.0; ") # change if you're on a different HPC

gtf = "/sc/orga/projects/ad-omics/data/references//hg38_reference/GENCODE/gencode.v30.annotation.gtf" # GENCODE V30 or whatever you have lying around
annotation = "data/gencode_v30.ensemblid_gene_name.tsv.gz"  # will have to reverse engineer - I think just matches EnsemblIDs  to gene names


dataCode = config["dataCode"]
inFolder = config["inFolder"]
outFolder = config["outFolder"]
support = config["metadata"]
region = config["region"]

metadata = pd.read_csv(support, sep = "\t")

# get out bam files and sample names
samples = metadata["sample_name"]
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

print( " * SGSeq pipeline")
print(" * %s comparisons to test" % len(condition_strings) )

rule all:
    input:
        expand(outFolder + "{comparison}/" + dataCode + "_{comparison}_sgseq_variant_type_table.tab", comparison = condition_strings)

# symlink bams to rename by sample name
rule symlinkBAMs:
    output:
        expand( inFolder + "{sample}.bam", sample = samples)
    run:
        for i in range(len(samples)):
            os.symlink(os.path.abspath(bam_files[i]), output[i], target_is_directory = False, dir_fd = None)

# index BAMs if not already
rule indexBAMs:
    input:
        inFolder + "{sample}.bam",
    output:
        inFolder + "{sample}.bam.bai"
    shell:
        "samtools index {input}"

# extract region from BAM files and index
rule extractRegion:
    input:
        bam = inFolder + "{sample}.bam",
        bai = inFolder + "{sample}.bam.bai"
    output:
        bam = inFolder + "{sample}_subset.bam",
        bai = inFolder + "{sample}_subset.bam.bai"
    shell:
        "samtools view -bh {input.bam} {region} > {output.bam};"
        "samtools index {output.bam}" 

# create new support using the region bams
rule createSGSeq_support:
    input:
        expand(inFolder + "{sample}_subset.bam", sample = samples)
    output:
        outFolder + dataCode + "_sgseq_region_support.tsv"
    run:
        lib_sizes = []
        for bam in metadata["file_bam"]:
            ls = 0
            samfile = pysam.AlignmentFile(bam, "rb")
            for read in samfile:
                ls += 1
            lib_sizes.append(ls)
        metadata["file_bam"] = input[0]
        metadata["lib_size"] = lib_sizes
        metadata.to_csv(output[0], sep = "\t", na_rep = "NA", index = False)    

# build SGSeq reference
rule step0:
    input:
        gtf = gtf
    output:
        "indexes/hg38.sgseqAnno.Rdata"
    params:
        script = "scripts/sgseq_step0.R"
    shell:
        "Rscript {params.script} --gtf {input.gtf} --sgseq.anno {output}"


# run SGSeq to find novel and annotated isoforms
rule step1b:
    input:
        sgseq_support = outFolder + dataCode + "_sgseq_region_support.tsv",
        bams = expand( inFolder + "{sample}_subset.bam" , sample = samples ),
        sgseqAnno = "indexes/hg38.sgseqAnno.Rdata"
    params:
        script = "scripts/sgseq_step1b.R"
    output:
       outFolder + dataCode + "_txf_novel.RData" 
    shell:
        "Rscript --vanilla {params.script} --support.tab {input.sgseq_support}"
        " --code {dataCode} --output.dir {outFolder} --gtf {gtf} "
        " --sgseq.anno {input.sgseqAnno}  "

# run DEXSeq to test
rule step2b:
    input:
        sgseq_support = outFolder + dataCode + "_sgseq_region_support.tsv", 
        txf_novel = outFolder + dataCode + "_txf_novel.RData"
    output:
        outFolder + "{comparison}/" + dataCode + "_{comparison}_res_clean_novel.tab"
    params:
        script = "scripts/sgseq_step2.R",
        code = dataCode
    shell:
        "Rscript --vanilla {params.script} --step step2b --support.tab {input.sgseq_support} "
        " --code {params.code} --output.dir {outFolder} --annotation {annotation}"

# create output tables
rule step3:
    input:
        sgseq_support = outFolder + dataCode + "_sgseq_region_support.tsv",
        step1b_output = outFolder + "{comparison}/" + dataCode + "_{comparison}_res_clean_novel.tab"
    output:
        outFolder + "{comparison}/" + dataCode + "_{comparison}_sgseq_variant_type_table.tab"
    params:
        createVarTable = "scripts/createVariantTable.R",
        findCentralExons = "scripts/findCentralExons.R"
    shell:
        "Rscript --vanilla {params.createVarTable} --step step2b --support.tab {input.sgseq_support} "
        "--code {dataCode} --output.dir {outFolder} ;"
        "Rscript --vanilla {params.findCentralExons} --step step2b --support.tab {input.sgseq_support} "
        "--code {dataCode} --output.dir {outFolder}"




