
Differential expression analysis
===
Limma
---
SAM
---
Biomart
---
Exercises
---

Exercise 1
===
Using the ALL dataset, do the following:

1. Choose a comparison that interests you - you'll make contrasts for it. 
For more information about all the pData fields, you can see [http://www.bioconductor.org/packages/release/data/experiment/manuals/ALL/man/ALL.pdf](http://www.bioconductor.org/packages/release/data/experiment/manuals/ALL/man/ALL.pdf)
2. Use genefilter to select only probes which map to gene symbols.
3. Encode your contrast of interest into a limma design and contrast matrix. Use limma alongside the filtered expression matrix to find the differentially expressed probes. 
4. Perform Gene Ontology enrichment analysis on the differentially expressed probes. Try using the molecular function ("MF") or cellular component ("CC") ontologies. Try performing GOEA using only upregulated or only downregulated genes. (Recall that this information can be found using the "decideTests" function.)
5. Go back to the gene filtering stage, and use genefilter to filter out probes with kOverA(10, 5.0). Re-run the rest of the pipeline. How do the results change?


