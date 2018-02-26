.onLoad <- function (lib, pkg) {
    if(.Platform$OS.type == "unix")
        microsimulation.init(pkg)
}

.onUnload <- function (libpath) {
    pkg <- utils::packageName()
    if(.Platform$OS.type == "unix")
        microsimulation.exit(pkg)
    library.dynam.unload(pkg, libpath)
}
