.onLoad <- function (lib, pkg) {
  prostata.init()
}

.onUnload <- function (libpath) {
  prostata.exit()
  library.dynam.unload("prostata", libpath)
}
