# Vignettes have been precompiled, due to long evaluation time

library(knitr)
opts_knit$set(base.dir = "vignettes")
knit("vignettes/low-dimensional-examples.Rmd.orig", "vignettes/low-dimensional-examples.Rmd")
knit("vignettes/demonstrating-the-hmer-package.Rmd.orig", "vignettes/demonstrating-the-hmer-package.Rmd")
knit("vignettes/emulationhandbook.Rmd.orig", "vignettes/emulationhandbook.Rmd")
knit("vignettes/stochasticandbimodalemulation.Rmd.orig", "vignettes/stochasticandbimodalemulation.Rmd")

library(devtools)
build_vignettes()
