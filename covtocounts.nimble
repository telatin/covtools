# Package
version       = "0.4.0"
author        = "Andrea Telatin"
description   = "covtocounts"
license       = "MIT"

# Dependencies
requires "nim >= 0.17.2", "c2nim >= 0.9.10", "docopt", "lapper", "hts"

srcDir = "src"
bin = @["covtocounts"]

task named_build, "custom build":
  mkdir "bin"
  exec "nimble c --out:bin/covtocounts src/covtocounts"
