# Package
version       = "2.0.0"
author        = "Andrea Telatin"
description   = "covtobed"
license       = "MIT"

# Dependencies
requires "nim >= 0.17.2", "c2nim >= 0.9.10", "docopt", "lapper", "hts"

srcDir = "src"
bin = @["covtocounts"]

task named_build, "custom build":
  mkdir "bin"
  exec "nimble c --out:bin/covtobed2 src/covtobed2"
