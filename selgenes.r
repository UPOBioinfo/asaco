#!/usr/bin/env Rscript

# Libraries and imports.
library(data.table)
source("selgenes_functions.R")

# Get input file names.
input_genes = list.files("input_genes", full.names = T)

# Prepare cluster for parallel execution.
if (Sys.info()['sysname'] == "Windows") {
  cl<<-makeCluster(spec = as.numeric(detectCores()), type="PSOCK")
  clusterExport(cl, c('data.table') )
} else {
  cl<<-makeCluster(spec = as.numeric(detectCores()), type="FORK")
}
registerDoParallel(cl)

# For each input file.
for (file in input_genes) {
  selgenes_file(file)
}

# Stop cluster.
stopCluster(cl)
