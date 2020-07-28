# Correlations for searching drugs again SARS-CoV-2
# AJPerez, March 2020
# Updated, July 2020

library(ggplot2)
library(reshape2)
library(topGO)
library(scales)
library(RColorBrewer)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(biomaRt)
library(rWikiPathways)
library(drugbankR)
library(tidyverse)
library(clipr)

setwd("/home/ajperez/Nextcloud/coronavirus")

# Parameters
TARGET_GENE <- 0 # 0 for all  
CORR_SIGN   <- 2 # 1=+, 2=-, 0=all
UPDOWN_SIGN <- 3 # 1=up, 2=down, 0=both, 3=off
FC          <- 1 # Fold change threshold
EXPRESSIONS <- "expressions/seeds/" # tr_itc
OUT_FOLDER  <- "results_def/"

print(list.files(path = EXPRESSIONS)) # expression matrices

# Databases' constants
dv <- "drugbank_5.1.5.db"
dids <- queryDB(type = "getIDs", db_path = dv)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")#, host = "uswest.ensembl.org")

# Genes and correlations
NEXP <- 0.9 # %ncols(expr_table)
GUPO_array <- list.files(path = EXPRESSIONS)
#GUPO_array <- GUPO_array[31:length(GUPO_array)]
#GUPO_array <- GUPO_array[c(19)] # <------------- REMAININGs ----------------
  if (TARGET_GENE != 0) {
  GUPO_array <- GUPO_array[TARGET_GENE]
}
if (CORR_SIGN == 0) {
  CORR_array <- c(1, -1)
} else if (CORR_SIGN == 1) {
  CORR_array <- 1
} else {
  CORR_array <- -1
}
if (UPDOWN_SIGN == 0) {
  UPDOWN_array <- c("up", "down")
} else if (UPDOWN_SIGN == 1) {
  UPDOWN_array <- "up"
} else if (UPDOWN_SIGN == 2) {
  UPDOWN_array <- "down"
} else {
  UPDOWN_array <- "off"
}
  
# Iterate genes
for (GUPO in GUPO_array) { 
  for (CORR in CORR_array) {
    for (UPDOWN in UPDOWN_array) { 
    
      GENE <- strsplit(GUPO, '\\(')[[1]][[1]]
      print(paste(GENE, CORR, UPDOWN, FC))

      # Organism variables and folder
      ORG  <- strsplit(GUPO, '_')[[1]][[3]]
      ORG1 <- "human"
      ORG2 <- "hsa"
      ORG3 <- "hsapiens_gene_ensembl"
      if (ORG == "mus musculus") {
        ORG1 <- "mouse"
        ORG2 <- "mmu"
        ORG3 <- "mmusculus_gene_ensembl"
      }

      ############
      # Files
      FILE <- paste0(EXPRESSIONS, GUPO, "/resultados.tsv")
      header <- read.csv(FILE, sep = "\t", nrows = 2, header = FALSE)
      header <- rbind(header, t(c("", "", paste0("exp", seq(1, ncol(header)-2)))))
      expr_table <- read.csv(FILE, sep = "\t", skip = 2, header = FALSE)
      names(expr_table) <- c("Ensembl", "Genename", paste0("exp", seq(1, ncol(expr_table)-2)))
      expr_table <- expr_table %>% filter(Ensembl != "") # empty gene (controls from   array)
      
      # Add gene type and experiment ID
      type_table <- read.csv(paste0("../galiciame/biomart_", ORG1, "_type.tsv"), sep = "\t")
      expr_table <- merge(expr_table, type_table, by.x = "Ensembl", by.y = "Gene.stable.ID")

      # Remove rows with one or more NA
      #expr_table <- na.omit(expr_table)
      #############
      
      # Reference gene
      ref_gene_index <- which(expr_table$Genename == GENE)
      expr_table <- expr_table[,order(as.numeric(expr_table[ref_gene_index,]), decreasing = T)]
      
      # Filter by activators/inhibitors (half of correlation)
      if (UPDOWN == "up") {
        temp <- expr_table[ref_gene_index, 4:ncol(expr_table)]
        temp <- temp[1, ] < FC
        expr_table <- expr_table %>% select(- names(temp[,temp]))
      } else if (UPDOWN == "down") {
        temp <- expr_table[ref_gene_index, 4:ncol(expr_table)]
        temp <- temp[1, ] > -FC
        expr_table <- expr_table %>% select(- names(temp[,temp]))
      } else {
        temp <- expr_table[ref_gene_index, 4:ncol(expr_table)]
        temp <- temp[1, ] > -FC & temp[1, ] < FC
        expr_table <- expr_table %>% select(- names(temp[,temp]))
      }
        
      # Follow with reference gene
      ref_gene <- expr_table %>% filter(Genename == GENE) # GENE="SMN2" | "SMN1" if smn1+smn2
      ref_expr <- as.numeric(ref_gene[,4:ncol(ref_gene)]) # -2 if smn1+smn2 and 4 -> 6
      
      # Correlation and number of NA
      fcor <- function(x) {
        c <- cor(as.numeric(x), ref_expr, use = "pairwise.complete.obs")
        return(c)
      }
      #expr_table2 <- expr_table
      expr_table$Correlation <- apply(expr_table %>% select(4:ncol(expr_table)), 1, fcor) # -2 if smn1+smn2 and 4 -> 6 (add smn2 to figure)
      expr_table$nexp_count <- apply(expr_table %>% select(5:ncol(expr_table)-1), 1, function(x) sum(!is.na(x))) # -2 if smn1+smn2 and 4 -> 6 (add smn2 to figure)
      
      # Correlation distribution
      corrs <- as.numeric(expr_table  %>% filter(nexp_count >= NEXP*length(ref_expr)) %>% pull(Correlation))
      mean_corr <- mean(corrs)
      if (CORR > 0) {
        CORR <- as.numeric(formatC(as.vector(quantile(corrs, 0.75) + IQR(corrs)), digits = 3, format = "f"))
        CORR <- 0.7
      } else {
        CORR <- as.numeric(formatC(as.vector(quantile(corrs, 0.25) - IQR(corrs)), digits = 3, format = "f"))
        CORR <- -0.7
      }
      
      # Filer by correlation and nexp
      filtered_expr_table <- expr_table  %>% filter((CORR < 0 & Correlation <= CORR &
                                                      nexp_count >= NEXP*length(ref_expr)) |
                                                      (CORR > 0 & Correlation >= CORR &
                                                      nexp_count >= NEXP*length(ref_expr)) |
                                                      Genename == GENE)
      
      filtered_expr_table <- expr_table %>% filter(Ensembl %in% readLines(paste0("gualberto/v3/output_genes/",
                                                   GENE, "_selgenes_inverse_scores.csv_0.01")) | Genename == GENE) # Gualberto (direct/inverse)
      #filtered_expr_table <- expr_table %>% filter(Ensembl %in% readLines(paste0("../galiciame/gualberto/smns/",
      #                                             GENE, "_selgenes_direct.csv")) | Genename == GENE) # Gualberto (direct/inverse)
      
      filtered_expr_table_corr <- filtered_expr_table
      filtered_expr_table %>% select(- c(Correlation, nexp_count))
      
      # Figure
      fig_expr_table <- melt(filtered_expr_table, id = c("Ensembl","Genename", "Gene.type"))
      fig_expr_table <- merge(fig_expr_table, t(header[c(1,3),3:ncol(header)]), by.x = "variable", by.y = "3")
      names(fig_expr_table)[6] <- "Atlas"
      fig_ref <- fig_expr_table %>% filter(Genename == GENE)
      
      p <- ggplot() +
        geom_hline(yintercept = 0, color = "blue") +
        geom_line(data = fig_expr_table, aes(variable, value, group = Genename, color = Gene.type), alpha = .8) +
        xlab("Experiments") +
        ylab("Log2-fold change") +
        theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.text = element_text(size=10.5),
              axis.text.x = element_text(angle = 60, hjust = 1), legend.position="none", #"top" 
                                         legend.title = element_blank()) +
        scale_y_continuous(breaks = seq(-10, 10, 2)) +
        geom_tile(data=fig_expr_table, aes(x = variable, y = -10, fill = Atlas)) +
        geom_line(data=fig_ref, aes(variable, value, group = Genename), size = 1.2) +
        guides(colour = guide_legend(nrow = 3))
      print(p)
      
      # Folder  
      FOLDER <- paste0(OUT_FOLDER, GENE, "_", CORR, "_", NEXP, "_", UPDOWN , "_", FC, "_", ORG1)
      if (!dir.exists(FOLDER)){
        dir.create(FOLDER)
      }
      
      # Save figure
      pdf(paste0(FOLDER, "/../", GENE, "_correlation.pdf"), width=16, height=8, paper='special')
      print(p)
      dev.off()
      
      # Correlation distribution
      h1 <- ggplot(expr_table %>% filter(nexp_count >= NEXP*length(ref_expr)), aes(Correlation)) + 
        #geom_histogram(binwidth = 20, color="grey20", aes(fill=ind), position = "dodge") +
        xlab("Correlation value") +
        ylab("Number of genes") +
        #scale_x_continuous(limits =c(-1, 1)) +
        geom_histogram(aes(y = ..count..), binwidth=0.01, alpha = .5) +
        #geom_density() +
        geom_vline(xintercept = mean_corr, color = "blue", linetype = "dashed", size = 1.25) +
        scale_x_continuous(breaks = seq(-1, 1, by = 0.1)) +
        theme_grey() +
        theme(legend.title = element_blank(), legend.position=c(0.82, 0.87), text = element_text(size=16))
      print(h1)
      
      pdf(paste0(FOLDER, "/corr_distribution.pdf"), width=16, height=8, paper='special')
      print(h1)
      dev.off()
      
      # Save file and count results
      filtered_expr_table_wo_SMN <- filtered_expr_table %>% filter(Genename != GENE)
      
      if (nrow(filtered_expr_table_wo_SMN) == 0) { next }
      
      write.table(filtered_expr_table_corr %>% filter(Gene.type == "protein_coding") %>% select(Ensembl, Genename, Correlation), 
                  file = paste0(FOLDER, "/genes_correlations.tsv"), quote = F, col.names = F, row.names = F, sep = "\t")
      
      write.table(filtered_expr_table_wo_SMN %>% filter(Gene.type == "protein_coding") %>% select(Ensembl), 
                  file = paste0(FOLDER, "/genes.ensembl"), quote = F, col.names = F, row.names = F)
      
      write.table(filtered_expr_table_wo_SMN %>% filter(Gene.type == "protein_coding") %>% select(Genename), 
                  file = paste0(FOLDER, "/genes.gn"), quote = F, col.names = F, row.names = F)
      
      write.table(filtered_expr_table_wo_SMN %>% group_by(Gene.type) %>%  summarise(count=n()), 
                  file = paste0(FOLDER, "/types.tsv"), quote = F, col.names = F, row.names = F, sep = "\t")
      genes <- as.vector(unlist(filtered_expr_table_wo_SMN %>% filter(Gene.type == "protein_coding") %>% select(Ensembl)))
      
      ##############  
      # Enrichment #
      ##############
next      
      # Files
      file_bg <- paste0("../galiciame/biomart_", ORG1, "_go2_protein.tsv")
      Nodes <- 50 # number of processes to show
      Ontology <- "GOES" #GO.P.ID : PFC (BP MF CC)
      
      #Create temp file
      data <- read.csv(file_bg, sep = "\t", header = TRUE, row.names = NULL)[,(c('ID', Ontology))]
      file_temp <- paste0(file_bg,"2")
      write.table(data, file = file_temp, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
      
      # Get background annotation
      GOesByID <- readMappings(file = file_temp)
      bg_genes <- names(GOesByID)
      
      compared_genes <- factor(as.integer(bg_genes %in% genes))
      names(compared_genes) <- bg_genes
      
      # Iterate through the two ontologies
      for (ONT in c("BP", "CC")) {
      next
        # Create topGO object
        GOdata <- new("topGOdata", ontology = ONT, allGenes = compared_genes,
                    annot = annFUN.gene2GO, gene2GO = GOesByID)
        asd <- unlist(Term(GOTERM))
      
        # Run Fisher test
        resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      
        # Create and print table with enrichment result
        allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = Nodes)
      
        # Different palettes
        palette <- c("#F52A2A", "#D561EA", "#61B0EA", "green", "#E89B57", "#E4EA61", "white") # alternative palette
        myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
      
        # Figure
        ########
        layout(t(1:2), widths=c(8,3))
        par(mar=c(4, .5, .7, .7), oma=c(3, 15, 3, 4), las=1)
      
        pvalue <- as.numeric(gsub("<", "", allRes$classicFisher)) # remove '<' symbols
        allRes$classicFisher <- pvalue
        max_value <- as.integer(max(-log(pvalue)))+1
        pv_range <- exp(-seq(max_value, 0, -1))
        allRes <- mutate(allRes, plot_id = paste(GO.ID, Term, sep = " - "))
      
        mylabels <- paste (allRes$GO.ID, "-",  asd[allRes$GO.ID])
        mybreaks <- 10^-(0:30)
      
        p <- ggplot(data=allRes, aes(x=reorder(plot_id, Significant), y=Significant)) +
          geom_bar(stat="identity", color="black", aes(fill=as.numeric(log(classicFisher))), size = 0.3)+
          geom_text(aes(label=mylabels), position=position_fill(vjust=0), hjust=0, fontface="bold", size = 5)+
          coord_flip() +
          theme(panel.background = element_blank(), panel.grid.major.x = element_line(colour = "darkgrey", size=0.75),
                panel.grid.minor.x = element_line(colour = "grey",size=0.75), axis.title.y=element_blank(), 
                axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
                axis.ticks.x =element_blank(), axis.line.y=element_blank(), axis.text=element_text(size=12)) +
          ylab("Number of genes") +
          guides(fill = guide_colourbar(barheight = 25, reverse=T)) +
          scale_fill_gradientn(name = "p-value", colours = myPalette(4), breaks = log(mybreaks), 
                               guide = guide_colourbar(reverse = TRUE), labels=mybreaks) +
          scale_y_continuous(breaks = seq(0, max(allRes$Significant), by = 20))
        print(p)
      
        # Save results
        if (exists("allRes") & nrow(allRes) > 4) {
          # Save table and figure 
          write.table(allRes, file = paste0(FOLDER, "/enrichment_", ONT, ".tsv"), quote = F, row.names = F, sep = "\t")
          pdf(paste0(FOLDER, "/enrichment_", ONT, ".pdf"), width=16, height=8, paper='special')
          print(p)
          dev.off()
      
          # Save list of genes by enriched GO
          #GOnames <- as.vector(allRes$GO.ID)
          #allGenes <- genesInTerm(GOdata, GOnames)
          #significantGenes <- list() 
          #for(x in 1:Nodes){
          #  significantGenes[[x]] <- allGenes[[x]][allGenes[[x]] %in% as.vector(genes)]
          #}
          #names(significantGenes) <- allRes$Term
          #new_folder <- paste0(FOLDER, "/", ONT, "/")
          #dir.create(new_folder)
          #for(x in 1:Nodes){
          #  write.table(significantGenes[x], quote = FALSE, sep = "\t", col.names = F, row.names = F,
          #              file = paste(new_folder, x, "-", str_replace_all(names(significantGenes)[x], "/", "_"), ".tsv", sep=""))
          #}
        }
      }
      
      ###################
      # ClusterProfiler #
      ###################
      
      entrez <- getBM(attributes="entrezgene_id", filters = 'ensembl_gene_id', values = genes, mart = ensembl)
      entrez <- as.vector(unlist(entrez['entrezgene_id']))
      
      ############
      # Reactome #
      ############
      xr = ""
      xr <- enrichPathway(gene = entrez, pvalueCutoff = 0.05, readable = T, organism = ORG1, pAdjustMethod = "fdr")
      
      if (exists("xr") & nrow(xr) > 1) {
        r <- dotplot(xr, showCategory = 15, font.size = 14)
        #scale_fill_gradientn(name = "p_adjust", colours = myPalette(4), guide = guide_colourbar(reverse = TRUE))
        #print(r)
      
        #pdf(paste0(FOLDER, "/../", GENE, ".pdf"), width=16, height=8, paper='special')
        png(paste0(FOLDER, "/../", GENE, ".png"), width = 800, height = 600)
        print(r)
        dev.off()
      
        #write.table(data.frame(xr$Description), quote = FALSE, col.names = F, row.names = F, file = paste0(FOLDER, "/../", GENE, ".txt"))
        write.table(entrez, quote = FALSE, col.names = F, row.names = F, file = paste0(FOLDER, "/../", GENE, ".txt"))
        
        # Write data
        #new_folder <- paste0(FOLDER, "/REACT/")
        #dir.create(new_folder)
        #for (i in 1:nrow(xr)) {
        #  xrgeneid <- unlist(strsplit(xr$geneID[i], split="/"))
        #  xrensembl <- getBM(attributes="ensembl_gene_id", filters = 'entrezgene_id', values = xrgeneid, mart = ensembl)
        #  xrensembl <- as.vector(unlist(xrensembl['ensembl_gene_id']))
        #  write.table(xrensembl, quote = FALSE, sep = "\t", col.names = F, row.names = F,
        #              file = paste(new_folder, i, "-", str_replace_all(xr$Description[i], "/", "_"), ".tsv", sep=""))
        #}
        write.table(as.data.frame(xr), file = paste0(FOLDER, "/enrichment_REACTOME.tsv"), quote = F, row.names = F, sep = "\t")
      }
      #next
      ################
      # WikiPathways #
      ################
      
      xw = ""
      ORG_WIKI <- ORG
      substr(ORG_WIKI, 1, 1) <- toupper(substr(ORG_WIKI, 1, 1))
      wp.gmt <- rWikiPathways::downloadPathwayArchive(organism = ORG_WIKI, format = "gmt")
      #listOrganisms()
      wp2gene <- clusterProfiler::read.gmt(wp.gmt)
      wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
      wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
      wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
      xw <- clusterProfiler::enricher(entrez, pAdjustMethod = "fdr",  pvalueCutoff = 0.1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

      if (exists("xw") & nrow(xw) > 4) {
        #as.data.frame(xw)
        b <- barplot(xw, showCategory=8, colorEdge = TRUE) + 
          ylab("Number of genes")
        #  scale_fill_gradientn(name = "p_adjust", colours = myPalette(4), guide = guide_colourbar(reverse = TRUE))
        print(b)
      
        pdf(paste0(FOLDER, "/wiki_plot.pdf"), width=16, height=8, paper='special')
        print(b)
        dev.off()
      
        # Write data
        #new_folder <- paste0(FOLDER, "/WIKI/")
        #dir.create(new_folder)
        #for (i in 1:nrow(xw)) {
        #  xwgeneid <- unlist(strsplit(xw$geneID[i], split="/"))
        #  xwensembl <- getBM(attributes="ensembl_gene_id", filters = 'entrezgene_id', values = xwgeneid, mart = ensembl)
        #  xwensembl <- as.vector(unlist(xwensembl['ensembl_gene_id']))
        #  write.table(xwensembl, quote = FALSE, sep = "\t", col.names = F, row.names = F,
        #              file = paste(new_folder, i, "-", str_replace_all(xw$Description[i], "/", "_"), ".tsv", sep=""))
        #}
        write.table(as.data.frame(xw), file = paste0(FOLDER, "/enrichment_WIKI.tsv"), quote = F, row.names = F, sep = "\t")
      }
      
      ############
      # DrugBank #
      ############
      #install.packages("remotes")
      #remotes::install_github("yduan004/drugbankR")
      
      # genes
      #filtered_expr_table_wo_SMN %>% filter(Gene.type == "protein_coding") %>% select(Genename)
      cgenes <- expr_table %>% filter(Ensembl %in% genes) %>% select(Genename) %>% pull()
      
      # KEGG
      #cgenes <- as.vector(unlist(expr_table %>% filter(Ensembl %in% xensembl) %>% select(Genename)))
      
      dtable <- queryDB(ids = dids, type = "getTargets", db_path = dv) %>% filter (t_gn_sym %in% cgenes)
      
      if (nrow(dtable) > 0) {
        # Write files
        #write_clip(dtable$t_Uni_id)
        write.table(dtable, file = paste0(FOLDER, "/drugbank.tsv"), quote = T, row.names = F)
        write.table(expr_table %>% filter(Genename %in% dtable$t_gn_sym) %>% select(Ensembl), 
                    file = paste0(FOLDER, "/genes_drugs.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
        write.table(expr_table %>% filter(Genename %in% dtable$t_gn_sym) %>% select(Ensembl, Correlation, nexp_count), 
                    file = paste0(FOLDER, "/genes_drugs_corr.tsv"), quote = F, row.names = F, col.names = F, sep = "\t")
      
        #drugs <- dtable %>% filter(t_Uni_id %in% c(dtable$t_Uni_id))
        #print(drugs)
        #drugs <- dtable %>% filter(t_Uni_id %in% "Q14376")
        drugs <- dtable$q_db_id
        fda <- queryDB(ids = drugs, type = "whichFDA", db_path = dv)
        write.table(fda, file = paste0(FOLDER, "/drugbank_fda.tsv"), quote = F, row.names = F, sep = "\t")
        next;
        # Match with pathways
        command <- paste0("grep -f ", FOLDER, "/", "genes_drugs.txt ", FOLDER, "/*/*")
        rgrep <- system(command, intern = TRUE)
        rgrep <- read.table(text = rgrep, sep = "/")
        rgrep <- separate(data = rgrep, col = V4, into = c("a", "b"), sep = ":")
        rgrep <- rgrep[,-c(1,2)]
        names(rgrep) <- c("Source", "Annotation", "Ensembl")
        write.table(rgrep, file = paste0(FOLDER, "/drugables_vs_annots.tsv"), quote = F, row.names = F, sep = "\t")
      }
    }
  }
}

