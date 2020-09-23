# covtools

_Work in progress_


### covtotarget
will count the _total nucleotide coverage_ per feature in a BED or GFF file using as input the output of [covtobed](https://github.com/telatin/covtobed) (also from STDIN).

```
covToTarget 0.1.0

  Usage: covtotarget [options] <Target> [<covtobed-output>]

Arguments:                                                                                                                                                 

  <Target>           the BED (or GFF) file containing regions in which to count reads
  <covtobed-output>  covtobed output, or STDIN if not provided

Options:

  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]
  -h, --help                   Show help
```

Example, can be used in a stream from the BAM emitter to covtobed:
```bash
cat input/mini.bam | covtobed | covtotarget input/mini.gff
```

Example, can be used on the bed output _produced by covtobed_:
```bash
covtotarget input/mini.gff input/mini.cov.bed 
```

### covtocounts
will count the _number of alignments_ in a BAM file per feature of a target BED or GFF file (basically, adds GFF support to `read-count` found in [nim-hts-tools](https://github.com/brentp/hts-nim-tools))
```
covToCounts 0.4.1

  Usage: covtocounts [options] <Target> <BAM-or-CRAM>

Arguments:                                                                                                                                                 

  <Target>       the BED (or GFF) file containing regions in which to count reads
  <BAM-or-CRAM>  the alignment file for which to calculate depth

Options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: ].
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]
  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]
  -n, --rpkm                   Add a RPKM column
  -l, --norm-len               Add a counts/length column (after RPKM when both used)
  --header                     Print header
  --debug                      Enable diagnostics    
  -h, --help                   Show help


## References
 * Brent Pedersen,  Aaron Quinlan, [hts-nim: scripting high-performance genomic analyses](https://academic.oup.com/bioinformatics/article/34/19/3387/4990493) (Bioinformatics)
 * Giovanni Birolo, Andrea Telatin, [covtobed: a simple and fast tool to extract coverage tracks from BAM files](https://joss.theoj.org/papers/10.21105/joss.02119) (JOSS)
