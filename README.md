# corrRNASeq

This is the repo for the corrRNA R package, which provides functionality for analyzing correlated RNA-sequencing data.

To install the package, the devtools package needs to first be installed:

install.packages(pkgs = c('devtools'))

Then, the lmerSeq package can be installed from github:

devtools::install_github("stop-pre16/lmerSeq", build_vignettes = T)

Next, the glmmADMB package needs to be installed (see https://glmmadmb.r-forge.r-project.org/ for more information):

install.packages("R2admb")

install.packages("glmmADMB", 
    repos=c("http://glmmadmb.r-forge.r-project.org/repos",
            getOption("repos")),
    type="source")

Finally, the corrRNASeq package can be installed:

devtools::install_github("ewynn610/corrRNASeq", build_vignettes = T)

To read the detailed vignette with an example analysis, run:

vignette("corrRNASeq-vignette")
