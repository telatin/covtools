# covtools

_Work in progress_

* **covtocounts** will count the _number of alignments_ in a BAM file per feature of a target BED or GFF file (basically, adds GFF support to `read-count` found in [nim-hts-tools](https://github.com/brentp/hts-nim-tools))
* **covtotarget** will count the _total nucleotide coverage_ per feature in a BED or GFF file using as input the output of [covtobed](https://github.com/telatin/covtobed) (also from STDIIN).


## References
 * Brent S Pedersen,  Aaron R Quinlan, [hts-nim: scripting high-performance genomic analyses](https://academic.oup.com/bioinformatics/article/34/19/3387/4990493) (Bioinformatics)
 * Giovanni Birolo, Andrea Telatin, [covtobed: a simple and fast tool to extract coverage tracks from BAM files](https://joss.theoj.org/papers/10.21105/joss.02119)
