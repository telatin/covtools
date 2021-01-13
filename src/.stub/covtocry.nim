import os
import hts
import docopt
import heapqueue
import strutils
import tables
import algorithm
import itertools
import iterutils

const
  version = "0.0.1"

# THIS SHOW HOW TO "GET NEXT ALIGNMENT"

template initClosure(id:untyped,iter:untyped) =
  let id = iterator():auto {.closure.} =
    for x in iter:
      yield x
 

proc covtobed(bam:Bam, mapq:uint8, eflag:uint16) =
  var
    coverage_ends = initHeapQueue[covEnd]()
    next_change   = 0
    chrSize       = initTable[string, int]()
    aln           : Record
    cov           = newCov()
    
  initClosure(nextaln,bam.items())
  while true:
    aln = nextaln()
    
    if aln.isNil:
      echo type(aln)
      echo "Bye!"
      quit()
    echo aln.chrom

proc main(argv: var seq[string]): int =
  let doc = format("""
  covToBed $version

  Usage: covtobed [options] [<BAM>]

Arguments:                                                                                                                                                 

  <BAM>          the alignment file for which to calculate depth

Options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
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
  """ % ["version", version])

  let args = docopt(doc, version=version, argv=argv)
  let mapq = parse_int($args["--mapq"])


  

  var
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])
    bam:Bam
    

  if $args["<BAM>"] != "nil":
    try:
      open(bam, $args["<BAM>"], threads=threads)
    except:
      stderr.writeLine("FATAL ERROR: Unable to read input file: ", $args["<BAM>"]) 
      quit(1)
  else:
    open(bam, "-", threads=threads)

  covtobed(bam, uint8(mapq), eflag)
  return 0

when isMainModule:
  var args = commandLineParams()
  discard main(args)