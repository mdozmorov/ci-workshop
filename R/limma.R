library(limma)
library(GEOquery)
library(gplots)
library(preprocessCore)
library(topGO)
library(hgu133a.db)

gds <- getGEO("GDS858")
pheno <- Columns(gds)
genes <- Table(gds)$IDENTIFIER

expr <- {
    tbl <- Table(gds)
    expr <- data.matrix(tbl[,3:ncol(tbl)])
    expr <- normalize.quantiles(expr)
    rownames(expr) <- tbl$ID_REF
    colnames(expr) <- colnames(tbl)[3:ncol(tbl)]
    t(scale(t(log2(expr))))
}

heatmap(expr, col=redgreen(75))

infection <- pheno$infection
genotype <- pheno$"genotype/variation"

design <- cbind(WT=1, 
                Infected=infection!="uninfected")
fit <- lmFit(expr, design)
fit <- eBayes(fit)
result <- topTable(fit)

geneList <- as.character(result$P.Value)
names(geneList) <- as.character(result$ID)
#genes <- result$ID[:30]
go.data <- new("topGOdata", ontology="BP",
               allGenes=rownames(expr),
               geneSel=geneList,
               nodeSize=10, annot=hgu133a.db)
