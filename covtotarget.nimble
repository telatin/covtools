# Package
version       = "0.1.1"
author        = "Andrea Telatin"
description   = "covtotarget"
license       = "MIT"

# Dependencies
requires "nim >= 0.17.2", "docopt", "lapper"

srcDir = "src"
bin = @["covtotarget"]

task named_build, "custom build":
  mkdir "bin"
  exec "nimble c --out:bin/covtotarget src/covtotarget"
