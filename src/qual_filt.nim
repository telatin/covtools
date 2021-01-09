import nimbioseq, os, strutils, iterutils

type
  Task = iterator (r: var int)

proc runTasks(t: varargs[Task]) =
  var ticker = 0
  var r: int
  while true:
    var x = t[ticker mod t.len]
    x(r)
    echo "#", r
    if finished(x): break
    inc ticker
 

iterator cReadSeqs(filename: string): Record {.closure.} =
  for rec in readSeqs(filename):
    yield rec

iterator letters: char {.closure.} =
  for c in 'a' .. 'z':
    yield c

iterator numbers(s: int): int {.closure.}=
  var n = s
  while true:
    yield n
    inc n


proc seqSummary(R1, R2: string, total = false) =
  let progName = split(getAppFilename(), "/")[getAppFileName().count("/")]
  if R1 != "" and R2 != "":
    for (s,z) in zip(cReadSeqs(filename: R1), cReadSeqs(filename: R2)):
      echo s.id
  else:
    echo progName & ": need file name"

when isMainModule: import cligen;dispatch(seqSummary)