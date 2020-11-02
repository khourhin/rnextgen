library(DESeq2)
library(biomaRt)
library(BiocParallel)

## DEFAULTS
## Normalized counts using counts function in DESEQ (could use rlog and vst as well)

## USING 10 threads by default
## Used only for DESeq wrapper ()
register(MulticoreParam(10))

check_counts_meta_tables = function(counts, meta){
    # Verify rownames in meta are the same as colnames in counts
    if (!isTRUE(all.equal(rownames(meta),colnames(counts)))){
    warning("Metadata doesn't seem to fit count matrix. Check both inputs")
    }
    else {
        message("OK: Count and meta tables seems to correspond")
    }
}

make_DR = function(counts, meta, design){

    check_counts_meta_tables(counts, meta)

    dds = DESeqDataSetFromMatrix(countData=counts, colData=meta, design=design)
    dds = DESeq(dds, parallel=TRUE)

    return(dds)
}

pairwise_comparison = function(dds, comps, meta_col_name){
    ## Perform all comparison specified in comps (a list of size 2
    ## lists). Using the groups specified in metadata table column 'meta_col_name'
    
    compa_res = list()

    for (comp in comps){
        compa_name = paste(c(comp[[1]], comp[[2]]), collapse="_VS_")
        message(compa_name)
        res_i = results(dds, contrast=c(meta_col_name, comp[[1]], comp[[2]]), alpha=0.05, parallel=TRUE, independentFiltering=TRUE, cooksCutoff=TRUE)
        summary(res_i)
        
        compa_res[[compa_name]] = res_i
    }
    return(compa_res)
}


export_counts = function(dds, prefix=''){

    # Export both raw and VST transformed counts

    ## rlog is accounting for lib size and apparently vst too (see deseq2 manual)
    ## rlog = rlog(dds, blind=FALSE)

    counts = counts(dds)
    vst_counts = assay(vst(dds, blind=FALSE))

    write.csv(counts, paste0(prefix, 'counts_raw.csv'))
    write.csv(vst_counts, paste0(prefix, 'counts_vst_norm.csv'))
}

export_results = function(res, prefix='', orderCol='padj'){
    ## Export the DE results table

    # Sort by the column specified
    if (orderCol != ''){
        res = res[order(res[,orderCol]), ]
    }
    
    write.csv(res, file=paste(prefix, 'degs_DESeq2.csv', sep="_"), row.names=TRUE)
    return(res)
}

export_pairwise = function(res, orderCol='padj'){
    ## For pairwise comparisons
    for (compa in names(res)){
        export_results(res[[compa]], prefix=compa, orderCol=orderCol)
    }
}

