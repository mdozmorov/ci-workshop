========================
Sequence analysis with R
========================

RNAseq is rapidly superseding microarrays
=========================================

(TODO Figure showing growth in datasets)

- Currently, microarrays are still cheaper
- ... but RNAseq provides far more information

RNAseq technical protocol
=========================

(TODO figure showing procedure here)

A typical RNAseq analysis pipeline
==================================

(TODO another figure)

Preprocessing reads with ShortRead
==================================

.. code-block:: r

    > fq <- readFastq("data/fq/")

Using Bowtie for alignment
==========================

- First install Bowtie2 and the appropriate reference index (here: UCSC hg19) from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
- And samtools from http://samtools.sourceforge.net

.. code-block:: bash

    $ bowtie2 -p 12 -x ~/ref/hg19 -U SRR452389.fastq.gz \
        | samtools view -bS - > SRR452389.bam

Postprocessing alignments with Rsamtools
========================================

.. code-block:: r

    > library(Rsamtools)
    > bam <- "data/bam/SRR452389.bam"
    > sorted.bam <- "data/bam/SRR452389.sorted.bam"
    > sortBam(bam, sorted.bam)
    > indexBam(sorted.bam)

- This can also be achieved in bash:

.. code-block:: bash

    $ samtools sort data/bam/SRR452389.bam data/bam/SRR452389.sorted
    $ samtools index data/bam/SRR452389.sorted.bam
