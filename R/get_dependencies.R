# Install required packages

required.packages = c()
required.bioC.packages = c("affy", "GEOquery",
    "gplots", "preprocessCore", "topGO", "hgu133a.db",
    "hgu95a2probe", "Biostrings", "Biobase",
    "gcrma", "limma", "annotate", "hgu95av2.db",
    "pamr", "e1071", "ShortRead", "Rsamtools",
    "GenomicRanges", "IRanges", "DESeq", "edgeR", "goseq",
    "SRAdb", "GEOmetadb")

for (pkg in required.packages) {
    if (!require(pkg, character.only=T))
        install.packages(pkg)
}
source("http://bioconductor.org/biocLite.R")
for (pkg in required.bioC.packages) {
    if (!require(pkg, character.only=T))
        biocLite(pkg)
}

# Create data directories and download data

ensure.downloaded = function(uri, subdir="") {
    path = paste("data", subdir, basename(uri), sep="/")
    if (!file.exists(path))
        download.file(uri,path)
}

ensure.dir = function(path) {
    if (!file.exists(path)) dir.create(path)
}

ensure.dir("data")
ensure.dir("data/cel")

ensure.downloaded("http://www-stat.stanford.edu/~tibs/PAM/Rdist/khan.txt")
ensure.downloaded("http://gbnci.abcc.ncifcrf.gov/geo/GEOmetadb.sqlite.gz")
ensure.downloaded("http://watson.nci.nih.gov/~zhujack/SRAmetadb.sqlite.gz")

su.base = "http://cals.arizona.edu/~anling/MCB516/data_lecture13/"
su.cels = c("Brain_1.CEL", "Brain_2.CEL", 
            "Fetal_brain_1.CEL", "Fetal_brain_2.CEL",
            "Fetal_liver_1.CEL", "Fetal_liver_2.CEL",
            "Liver_1.CEL", "Liver_2.CEL")
for (cel in su.cels)
    ensure.downloaded(paste(su.base,cel,sep=""), "cel")
