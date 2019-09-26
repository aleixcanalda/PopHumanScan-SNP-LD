install.packages("gsubfn")
library(gsubfn)
install.packages("disgenet2r")
library(disgenet2r)

source("https://bioconductor.org/biocLite.R")
biocLite("KEGGprofile")
library(KEGGprofile)

source("http://bioconductor.org/biocLite.R")
biocLite( c("pathview", "gage", "gageData") )
library(pathview)
library(gage)
library(gageData)

library(biomaRt)

#to start off we want each gene in a different row
gene.table <- tibble()

for (gene in 1:nrow(useful.table)) {
  print(paste("-->", gene, sep=" "))
  more.genes <- unlist(strsplit(as.character(useful.table[gene,1]),","))
  if (length(more.genes)>1) {
    for (hi in more.genes) {
      print(paste("-->", hi, sep=" "))
      GeneID <- unlist(strapplyc(more.genes, hi, simplify = TRUE))
      GeneID <- cbind(GeneID, useful.table[gene,2:6])
      gene.table <- rbind(gene.table,GeneID)
    }
  }
  else {
    gene.table <- rbind(gene.table,useful.table[gene,])
  }
}
#now we need only the types that are genes
c.gene.table <- tibble()
for (houhou in 1:nrow(gene.table)) {
  if (gene.table[houhou,2] == "Gene") {
    c.gene.table <- rbind(c.gene.table,gene.table[houhou,])
  }
  else {next}
} 

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
temp <- getBM(attributes = c('entrezgene_id','external_gene_name'), 
              filters = c('external_gene_name'), 
              values = list(unlist(c.gene.table[,1])),
              mart = mart)
#from 478 to 385 rows because it eliminated all noncoding genes and other gene symbols that don't seem to have an entrez id

#now we want to analyze the genes and look which pathway they belong to
kegg <- find_enriched_pathway(unlist(temp[,1]),species = "hsa")


write.table(gene.table, file="/home/aleixcanalda/Desktop/AleixCanalda/VEP/genetable.txt",
            append = FALSE, quote = F, sep = "\t",
            col.names = TRUE, row.names = F)
         