library(compiler)
library(data.table)
library(ggplot2)
library(ggpubr)
library(doParallel)
library(psych)

selgenes_file <- function (input_file) {
  
  # Extract, transform and load gene fold changes from experiments.
  data = etl_experiments(input_file)
  xgene = data$xgene
  
  # Select genes by using all available methods.
  res = selgenes(data, xgene)
  
  # Write results into the output file.
  fwrite(data.table(res$direct), paste0("output_genes/", tools::file_path_sans_ext(basename(input_file)),"_selgenes_direct.csv"), row.names = F, col.names = F)
  fwrite(data.table(res$inverse), paste0("output_genes/", tools::file_path_sans_ext(basename(input_file)),"_selgenes_inverse.csv"), row.names = F, col.names = F)
  fwrite(data.table(res$direct_scores), paste0("output_genes/", tools::file_path_sans_ext(basename(input_file)),"_selgenes_direct_scores.csv"), row.names = F, col.names = T)
  fwrite(data.table(res$inverse_scores), paste0("output_genes/", tools::file_path_sans_ext(basename(input_file)),"_selgenes_inverse_scores.csv"), row.names = F, col.names = T)
  
}



etl_experiments <- cmpfun(function (filename, genetypes_filename = "data/biomart_human_type.tsv") {
  # Read fold change values from a experiments file.
  data = fread(filename, header=T)
  
  # Ignore the three first rows. 
  data = data[4:.N]
  
  # Convert data columns to numeric.
  data[,3:ncol(data)] = data[, lapply(.SD, as.numeric), .SDcols = 3:ncol(data)]
  
  # Determine the xgene based upon file name.
  xgene = tools::file_path_sans_ext(basename(filename))
  if (!startsWith(xgene, "ENSG")) {
    xgene = data[V2 == xgene, V1]  
  }
  
  # Rename the column of the gene name.
  setnames(data, "V1", "Gene")
  data[, V2 := NULL]
  
  # Filter genes with type "protein_coding".
  types = fread(genetypes_filename)
  types = types[, `Gene type`:=as.factor(`Gene type`)]
  coding_genenames = types[`Gene type` == "protein_coding"]$`Gene stable ID`
  data = data[Gene %in% coding_genenames]

  # Determination of duplicated genes.
  dup = data[duplicated(data[,.(Gene)]),Gene]
  genes = unique(data[,.(Gene)])
  
  # Averaging fold change values for repeated genes.
  data = data[, lapply(.SD, mean), by=Gene, .SDcols = 3:ncol(data)]
  
  # Transpose the data table (experiments -x- genes).
  datat = as.data.table(t(data[,-1]))
  setnames(datat, names(datat), t(genes))
  
  return(list(exp2genes = datat, genes2exp = data, duplicated = dup, different_genes = genes, xgene = xgene))
  
})

selgenes <- function(data, xgene) {
  res = selgenes_method_SFCS(data, xgene)
  return(res)
}

metric_SFCS <- cmpfun(function (fc_xgene, fc_other) {
  r = fc_xgene * fc_other
  p = r[r[[1]] > 0,.N] / r[,.N]
  n = r[r[[1]] < 0,.N] / r[,.N]
  z = 1 - p - n
  return(data.table(p,n,z))
})

selgenes_method_SFCS <- function(data_pack, xgene, method = "pvalue", sd_factor = 2) {
  
  # Extract only the used variables from the data_pack.
  data = data_pack$genes2exp
  datat = data_pack$exp2genes

  # Assess the score for each gene.
  score <- foreach(col=1:ncol(datat), .combine=rbind) %dopar% {
    metric_SFCS(datat[,xgene,with=F], datat[,col,with=F])
  }
  
  # Build the result table.
  ppv = 1 - pnorm(scale(score$p))
  npv = 1 - pnorm(scale(score$n))
  genes_table = data.table(gene = names(datat), score = score, ppv = ppv, npv = npv)
  setnames(genes_table, c("ppv.V1","npv.V1"), c("ppv","npv"))
  genes_table = genes_table[score.p != 1]
  
  # Selection of relevant related genes to xgene.
  if (method == "pvalue") {
    # (by p-value).
    sel_genes_direct = genes_table[ppv <= 0.05]
    sel_genes_inverse = genes_table[npv <= 0.05]
  } else {
    # (by mean + sd * sd_factor).
    threshold_direct = mean(score$p) + sd_factor * sd(score$p)
    threshold_inverse = mean(score$n) + sd_factor * sd(score$n)
    sel_genes_direct = genes_table[score.p >= threshold_direct]
    sel_genes_inverse = genes_table[score.n >= threshold_inverse]
  }

  # Sort by their scores.
  sel_genes_direct = sel_genes_direct[order(-score.p)]
  sel_genes_inverse = sel_genes_inverse[order(-score.n)]
  
  return(list(direct = sel_genes_direct[,gene], inverse = sel_genes_inverse[,gene], scores = genes_table, 
              direct_scores = sel_genes_direct, inverse_scores = sel_genes_inverse))
  
}
