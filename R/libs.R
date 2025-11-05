suppressPackageStartupMessages(library("cmdstanr"))
suppressPackageStartupMessages(library("config"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("git2r"))

# for devtools:
# https://forum.posit.co/t/failling-install-devtools-error-onload-failed-in-loadnamespace-for-pkgload/64787

# relies on ‘gt’ version ‘1.0.0.9000’, use 
# devtools::install_github("rstudio/gt")
# to install.
suppressPackageStartupMessages(suppressWarnings(library("gt")))
suppressPackageStartupMessages(library("jsonlite"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("logger"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("pracma"))
suppressPackageStartupMessages(library("marginaleffects"))
suppressPackageStartupMessages(library("ggthemes"))
suppressPackageStartupMessages(library("bayesplot"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("qs"))
suppressPackageStartupMessages(library("extraDistr"))
suppressPackageStartupMessages(library("mice"))
suppressPackageStartupMessages(library("pbapply"))
suppressPackageStartupMessages(library("poisson"))
suppressPackageStartupMessages(library("tictoc"))
suppressPackageStartupMessages(library("ggh4x"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(mvtnorm))

suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(nlme))

# install.packages(
#   c(
#     "marginaleffects",
#     "mice"
#   ),
#   Ncpus = 20
# )
