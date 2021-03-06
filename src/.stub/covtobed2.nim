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
  
]# 
template initClosure(id,iter:untyped) =
  let id = iterator():auto{.closure.} =
    for x in iter:
      let strand = if x.flag.reverse: '+' 
      else: '-'

      stderr.writeLine( " +---> ", x.qname, "\t", x.chrom , ":", x.start, ":", x.stop, " ", strand)
      yield x

type
  coverage = ref object
    forward: int
    reverse: int

proc doAssert(condition: bool, message: string) =
  if condition == false:
    stderr.writeLine("ERROR: ", message)
    quit(1)

proc newCov(f = 0, r = 0): coverage =
  coverage(forward: f, reverse: r)

proc inc(c: coverage, reverse=false) =
  if reverse == false:
    c.forward += 1
  else:
    c.reverse += 1

proc dec(c: coverage, reverse=false) =
  if reverse == false:
    c.forward -= 1
  else:
    c.reverse -= 1

proc tot(c: coverage): int =
  c.forward + c.reverse

proc topStop(q: HeapQueue): int =
  return q[0].stop


proc topReverse(q: HeapQueue): bool =
  return not q[0].reverse

proc empty(q: HeapQueue): bool =
  if len(q) == 0:
    return true
  else:
    return false

# Class that stores info about the end position of alignments, used in the alignment queue
type 
  covEnd = ref object
    stop: int
    reverse: bool



#TODO: fix this hack to make it returning the maximum item
proc `<`(a, b: covEnd): bool = a.stop > b.stop

#proc `>`(a, b: covEnd): bool = a.stop < b.stop
proc writeCoverage(chrName: string, last_pos, next_change: int, c: coverage) =
  echo "[", chrName, "]\t", last_pos, "\t", next_change, "\t", c.tot(), "x\t", c.forward, ',', c.reverse 
  discard


proc covtobed(bam:Bam, mapq:uint8, eflag:uint16) = 
  var
    next_change   = 0
    chrSize       = initTable[string, int]()
    aln           : Record

  initClosure(nextAlignment,bam.items())

  aln = nextAlignment()
  for reference in bam.hdr.targets():
    stderr.writeLine("===", reference.name, '\t', reference.length)
    chrSize[reference.name] = int(reference.length)
    
 
    var 
      more_alignments = not aln.isNil
      last_pos = 0
      more_alignments_for_ref = more_alignments and aln.chrom == reference.name
      coverage_ends = initHeapQueue[covEnd]()
      cov           = newCov()

    while true:

      # calculate next change
      if more_alignments_for_ref:
        if coverage_ends.empty():
          next_change = int(aln.start)
        else:
          next_change = min(int(aln.start), coverage_ends.topStop())
      else:
        if coverage_ends.empty():
          next_change = int(reference.length)
        else:
          next_change = coverage_ends.topStop()

      # output coverage ...
      doAssert(coverage_ends.len() == cov.tot(), "coverage not equal to queue size")
      writeCoverage(reference.name, last_pos, next_change, cov)
      
      stderr.writeLine("-+-  next=", next_change, "\tMore.aln=", more_alignments, "|", more_alignments_for_ref,";Cov=", cov.tot(), ";Size=", len(coverage_ends))
      if more_alignments:
        stderr.writeLine( " +-> more aln >pos=", aln.start)
      
      
      # increment coverage with aln starting here
      while more_alignments_for_ref and (next_change == aln.start):
        if "physical_coverage" == "not":
          echo "to be defined"
          #[						
            if (alignment.InsertSize > 0) {
						        debug cerr << "   [phy] pos:" << alignment.Position << " size:" << alignment.InsertSize << endl;
							coverage_ends.push({alignment.Position + alignment.InsertSize, alignment.IsReverseStrand()});
							coverage.inc(alignment.IsReverseStrand());
						}
            ]#
        else:
          coverage_ends.push( covEnd(stop: aln.stop, reverse: true) )
          cov.inc(aln.flag.reverse)

        aln = nextAlignment()
        more_alignments = not aln.isNil
        more_alignments_for_ref = more_alignments and aln.chrom == reference.name

      #decrement coverage with alignments that end here
      while ( not coverage_ends.empty()  and next_change == coverage_ends.topStop()):
        cov.dec(coverage_ends.topReverse())
        discard coverage_ends.pop()
      
      last_pos = next_change
      # End chromosome loop
      if last_pos == int(reference.length) or not more_alignments:
        stderr.writeLine("Lastpos=", last_pos, "==", int(reference.length), " OR NOT more_aln=", more_alignments,)
        doAssert(cov.tot()==0, "coverage not null at the end of chromosome " & reference.name & ": " & $cov.tot() & " = " & $cov.forward & "+" & $cov.reverse) 
        doAssert(coverage_ends.len() == 0, "coverage queue not null at the end of chromosome "  & reference.name & ": " & $coverage_ends.len())
        break
    
    if not coverage_ends.len() == 0:
      stderr.writeLine("Coverage not zero when expected. Try samtools fixmate.")
      raise
  
  #if more_alignments:
  #  stderr.writeLine("Is the BAM sorted?")
  #ß  raise


proc main(argv: var seq[string]): int =
  let doc = format("""
  covToBed $version

  Usage: covtobed [options] [<BAM>]

Arguments:                                                                                                                                                 

  <BAM>          the alignment file for which to calculate depth

Options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -p, --physical               Calculate physical coverage
  -s, --stranded               Report coverage separate by strand
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