
# SGSeq-pipeline
gtf = ""

# build SGSeq reference
rule step0:
    output:
        "GENCODE_v30_hg38.sgseqAnno.Rdata"
    shell:
        "Rscript $step0 --gtf {gtf} --sgseq.anno {output}"
