##' @import rje
##' @importFrom utils head
##' @importFrom stats rgamma

.onLoad <- function(libname, pkgname) {
  vig_list = tools::vignetteEngine(name='knitr', package = 'knitr')
  vweave <- vig_list[['knitr::knitr']][c('weave')][[1]]
  vtangle <- vig_list[['knitr::knitr']][c('tangle')][[1]]
  tools::vignetteEngine(pkgname, weave = vweave, tangle = vtangle,
                        pattern = "[.]Rmd$", package = pkgname)
  #register_vignette_engines(pkgname)
}
