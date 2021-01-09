import os
import hts
import docopt
import heapqueue
import strutils
import tables
import algorithm

const
  version = "0.0.1"
#[
fixedStep  chrom=chrN start=position  step=stepInterval [span=windowSize]
  dataValue1
  dataValue2
  ... etc ...
]#

var
  verbose = false



 


proc bam_to_wig(bam:Bam, mapq:uint8, eflag:uint16, step:int) =
  var
    chromosomes = newSeq[string]()
    lastpos:int
    nextpos:int
    coverage:int

  for align in bam.items():
    if align.chrom notin chromosomes:
      # new chromosome!
      echo "fixedStep chrom=", align.chrom, " start=", align.start, " step=", step
      chromosomes.add(align.chrom)


    echo '#', align.start, "-", align.stop


proc main(argv: var seq[string]): int =
  let doc = format("""
  covToWig $version

  Usage: covtowig [options] [<BAM>]

Arguments:                                                                                                                                                 
  <BAM>          the alignment file for which to calculate depth

Options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]
  --verbose                    Print a lot of info                      
  -h, --help                   Show help
  """ % ["version", version])

  let args = docopt(doc, version=version, argv=argv)
  let mapq = parse_int($args["--mapq"])
  verbose  = args["--verbose"]
   
  var
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])
    bam:Bam
    

  if $args["<BAM>"] != "nil":
    try:
      open(bam, $args["<BAM>"], threads=threads)
      if verbose:
        stderr.writeLine( "Reading BAM: ", $args["<BAM>"])
    except:
      stderr.writeLine("FATAL ERROR: Unable to read input file: ", $args["<BAM>"]) 
      if verbose:
        stderr.writeLine( "Reading STDIN... ")
      quit(1)
  else:
    open(bam, "-", threads=threads)

  

  bam_to_wig(bam, uint8(mapq), eflag, 50)
  return 0

when isMainModule:
  var args = commandLineParams()
  discard main(args)