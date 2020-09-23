import os
import system
import hts
import docopt
import lapper
import strutils
import tables
import algorithm

const
  version = "0.1.0"
#[
  **covToTarget**, part of MAGENTA Flow
  based on count-reads in the "hts-nim-tools" suite by Brent Pedersen
  see: "https://github.com/brentp/hts-nim-tools"
  Static binary thanks to  "https://github.com/brentp/hts-nim"

 
   
  0.1.0   Initial release
]#


var
  gffIdentifier = "ID"
  gffSeparator  = ";"
  gffField      = "CDS"

type
  region_t = ref object
    chrom: string
    start: int
    stop: int
    name: string
    count: int

proc inc_count(r:region_t) = inc(r.count)
proc start(r: region_t): int {.inline.} = return r.start
proc stop(r: region_t): int {.inline.} = return r.stop
proc `$`(m:region_t): string = return "($#-$#:$#, $#)" % [$m.chrom, $m.start, $m.stop, $m.name]
proc overlapLen(reference, feature: region_t): int =
    if feature.start < reference.start:
        return feature.stop - reference.start
    elif feature.stop > reference.stop:
        return reference.stop - feature.start
    else:
        return feature.stop - feature.start
 
# Converts a GFF line to region object
proc gff_line_to_region(line: string): region_t =
  var
   cse = line.strip().split('\t')

  if len(cse) < 5:
    stderr.write_line("[warning] skipping GFF line (fields not found):", line.strip())
    return nil

  # Skip non CDS fields (or user provided)
  if cse[2] != gffField:
    return nil

  var
    s = parse_int(cse[3])  - 1
    e = parse_int(cse[4])
    reg = region_t(chrom: cse[0], start: s, stop: e, count:0)
  
  # In the future, 8th field could be requireed [TODO]
  if len(cse) == 9:
    for gffAnnotPart in cse[8].split(gffSeparator):
      if gffAnnotPart.startsWith(gffIdentifier):
        reg.name = gffAnnotPart.split("=")[1] 
        break
  return reg

# Convert a BED line to region object
proc bed_line_to_region(line: string): region_t =
  var
   cse = line.strip().split('\t', 5)

  if len(cse) == 9:
    stderr.writeLine("[warning] GFF format detected.")
    return gff_line_to_region(line)
    

  if len(cse) < 3:
    stderr.write_line("[warning] skipping bad bed line:", line.strip())
    return nil
  var
    s = parse_int(cse[1])
    e = parse_int(cse[2])
    reg = region_t(chrom: cse[0], start: s, stop: e, count:0)
  if len(cse) > 3:
   reg.name = cse[3]
  return reg

# Convert BED file to table
proc bed_to_table(bed: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("track "):
      continue
    if $kstr.s[0] == "#":
      continue
    var v = bed_line_to_region($kstr.s)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)
  
  for chrom, ivs in bed_regions.mpairs:     # since it is read into mem, can also well sort. (BP)
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions



proc gff_to_table(bed: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("##FASTA"):
      break
    if $kstr.s[0] == "#":
      continue

    var v = gff_line_to_region($kstr.s)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions
 
 
proc processCoverage(f: File, target: TableRef[string, seq[region_t]]) =
    var
        line: string
        chroms = initCountTable[string]()
        featCount = initCountTable[string]()
        lap: Lapper[region_t] 
        res = new_seq[region_t]() 
    #for chromosome in target.keys:

    while f.readLine(line):
      let interval = bed_line_to_region(line)
      chroms.inc(interval.chrom)
      if chroms[interval.chrom] == 1:
          lap = lapify(target[interval.chrom])
          res = @[]
      if lap.seek(interval.start, interval.stop, res):

        let counts = parseInt(interval.name) * overlapLen(res[0], interval)
        if counts > 0:
          featCount.inc(res[0].name, counts)
    
    for feature in featCount.keys:
      echo feature, "\t", featCount[feature]
        

      
#[
  proc each_seek[T: Interval](
    L: var Lapper[T]; 
    start: int; 
    stop: int; f
    n: proc (v: T)) {..}
call fn(x) for each interval x in L that overlaps start..stop this assumes that subsequent calls to this function will be in sorted order
]#

proc main(argv: var seq[string]): int =
  
  let doc = format("""
  covToCounts $version

  Usage: covtocounts [options] <Target> [<covtobed-output>]

Arguments:                                                                                                                                                 

  <Target>           the BED (or GFF) file containing regions in which to count reads
  <covtobed-output>  covtobed output, or STDIN if not provided

Options:

  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]
  -h, --help                   Show help
  """ % ["version", version])

  let args = docopt(doc, version=version, argv=argv)
  var prokkaGff : bool = args["--gff"]
   
  gffIdentifier = $args["--id"]
  gffField      = $args["--type"]
 
  if ($args["<Target>"]).contains(".gff"):
    prokkaGff = true

  var regions = if prokkaGff == true: gff_to_table($args["<Target>"])
                 else: bed_to_table($args["<Target>"])
   
  var
    f: File

  if $args["<covtobed-output>"] != "nil":
    try:
      f = open($args["<covtobed-output>"])
    except:
      stderr.writeLine("FATAL ERROR: Unable to read input file: ", $args["<covtobed-output>"]) 
  else:
    f = stdin

  processCoverage(f, regions)
  #print_alignments_count(bam, uint8(mapq), eflag, regions)
  return 0

when isMainModule:
  var args = commandLineParams()
  discard main(args)