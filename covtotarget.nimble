# Package
version       = "0.4.0"
author        = "Andrea Telatin"
description   = "covtotarget"
license       = "MIT"

# Dependencies
requires "nim >= 0.17.2", "c2nim >= 0.9.10", "docopt", "lapper"

srcDir = "src"
bin = @["covtotarget"]

task named_build, "custom build":
  mkdir "bin"
  exec "nimble c --out:bin/covtotarget src/covtotarget"
