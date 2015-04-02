
## ----vignetteSetup, echo=FALSE, message=FALSE, warning = FALSE-----------
## Track time spent on making the vignette
startTime <- Sys.time()

## Bib setup
library('knitcitations')

## Load knitcitations with a clean bibliography
cleanbib()
cite_options(hyperlink = 'to.doc', citation_format = 'text', style = 'html')
# Note links won't show for now due to the following issue
# https://github.com/cboettig/knitcitations/issues/63

## Write bibliography information
bibs <- c(knitcitations = citation('knitcitations'),
    knitrBootstrap = citation('knitrBootstrap'), 
    knitr = citation('knitr')[3],
    rmarkdown = citation('rmarkdown'),
    R = citation(),
    BiocParallel = citation('BiocParallel')
)

write.bibtex(bibs,
    file = 'BiocParallel-knitrBootstrap.bib')
bib <- read.bibtex('BiocParallel-knitrBootstrap.bib')

## Assign short names
names(bib) <- names(bibs)

library('BiocParallel')
runSlow <- TRUE


## ----'install', eval = FALSE---------------------------------------------
## ## Install BiocParallel
## source('http://bioconductor.org/biocLite.R')
## biocLite('BiocParallel')


## ----'install2', eval = FALSE--------------------------------------------
## ## If needed:
## # install.packages('devtools')
## devtools::install_github('jimhester/knitrBootstrap')


## ----'rhelp', eval = FALSE-----------------------------------------------
## help(package = 'BiocParallel')
## help(package = 'knitrBootstrap')


## ------------------------------------------------------------------------
plot(y = 10 / (1:10), 1:10, xlab = 'Number of cores', ylab = 'Time',
    main = 'Ideal scenario', type = 'o', col = 'blue',
    cex = 2, cex.axis = 2, cex.lab = 1.5, cex.main = 2, pch = 16)


## ------------------------------------------------------------------------
plot(y = 10 / (1:10), 1:10, xlab = 'Number of cores', ylab = 'Time',
    main = 'Reality', type = 'o', col = 'blue',
    cex = 2, cex.axis = 2, cex.lab = 1.5, cex.main = 2, pch = 16)
lines(y = 10 / (1:10) * c(1, 1.05^(2:10) ), 1:10, col = 'red',
    type = 'o', cex = 2)


## ------------------------------------------------------------------------
birthday <- function(n) {
    m <- 10000
    x <- numeric(m)
    for(i in 1:m) {
        b <- sample(1:365, n, replace = TRUE)
        x[i] <- ifelse(length(unique(b)) == n, 0, 1)
    }
    mean(x)
}


## ----, eval = runSlow----------------------------------------------------
system.time( lapply(1:100, birthday) )


## ----, eval = runSlow----------------------------------------------------
library('doMC')
registerDoMC(2)
system.time( x <- foreach(j = 1:100) %dopar% birthday(j) )


## ----, eval = runSlow----------------------------------------------------
library('BiocParallel')
system.time( y <- bplapply(1:100, birthday) )


## ------------------------------------------------------------------------
registered()


## ------------------------------------------------------------------------
## Test in serial mode
system.time( y.serial <- bplapply(1:10, birthday,
    BPPARAM = SerialParam()) )

## Try Snow
system.time( y.snow <- bplapply(1:10, birthday, 
    BPPARAM = SnowParam(workers = 2)) )


## ----, eval = FALSE------------------------------------------------------
## cluster.functions = makeClusterFunctionsSGE("~/simple.tmpl")
## mail.start = "none"
## mail.done = "none"
## mail.error = "none"
## staged.queries = TRUE
## fs.timeout = 10


## ----, eval = FALSE------------------------------------------------------
## library('BiocParallel')
## library('BatchJobs')
## 
## # define birthday() function
## 
## ## Register cluster
## funs <- makeClusterFunctionsSGE("~/simple.tmpl")
## param <- BatchJobsParam(workers = 10, resources = list(ncpus = 1),
##     cluster.functions = funs)
## register(param)
## 
## ## Run
## system.time( xx <- bplapply(1:100, birthday) )
## 
## ## Jobs spend a little bit of time in the queue
## #   user  system elapsed
## #  0.597   0.350  31.644


## ----, eval = FALSE------------------------------------------------------
## rmarkdown::render('myFile.Rmd')


## ----'citation'----------------------------------------------------------
## Citation info
citation('BiocParallel')



## ----'citation2'---------------------------------------------------------
citation('knitrBootstrap')


## ----createVignette, eval=FALSE------------------------------------------
## ## Create this page
## library('rmarkdown')
## render('index.Rmd')
## 
## ## Clean up
## file.remove('BiocParallel-knitrBootstrap.bib')
## 
## ## Extract the R code
## library('knitr')
## knit('index.Rmd', tangle = TRUE)


## ----reproducibility1, echo=FALSE----------------------------------------
## Date the vignette was generated
Sys.time()


## ----reproducibility2, echo=FALSE----------------------------------------
## Processing time in seconds
totalTime <- diff(c(startTime, Sys.time()))
round(totalTime, digits=3)


## ----reproducibility3, echo=FALSE----------------------------------------
## Session info
library('devtools')
options(width = 120)
session_info()$platform


## ----reproducibility4, echo=FALSE----------------------------------------
## Session info packages
session_info()$packages


## ----vignetteBiblio, results = 'asis', echo = FALSE, warning = FALSE-----
## Print bibliography
bibliography()


