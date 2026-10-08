#!/usr/bin/env Rscript
# Inverse of tools/rda_to_arma.R: build the SAIGE .rda null model (list `modglmm`)
# that R SAIGE 1.5.2 step 2 (ReadModel in readInGLMM.R) reads, from the .arma +
# nullmodel.json directory our C++ step 1 writes -- so that R and C++ step 2
# run on the *identical* null model.
#
#   arma_to_rda.R <armaDir> <out.rda>
#
# Works for full-GRM and sparse-GRM models (not LOCO). For a sparse-GRM model, give R step 2
# the same sparse GRM files (--sparseGRMFile, --sparseGRMSampleIDFile) that step 1 used.
#
# Fields R step 2 uses (readInGLMM.R / SAIGE_Test_main.R): theta, traitType, y, X,
# fitted.values, residuals, sampleID, LOCO, offset, obj.noK{XV,XVX,XXVX_inv,
# XVX_inv,S_a,XVX_inv_XV,V}. Everything else ReadModel sets to NULL or defaults.
args <- commandArgs(trailingOnly = TRUE)
IND <- args[1]; OUT <- args[2]

rarma <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  hdr <- readLines(con, n = 2)
  stopifnot(hdr[1] == "ARMA_MAT_BIN_FN008")
  d <- as.integer(strsplit(hdr[2], " ")[[1]])
  v <- readBin(con, what = "double", n = d[1] * d[2], size = 8, endian = "little")
  matrix(v, nrow = d[1], ncol = d[2])
}
js <- paste(readLines(file.path(IND, "nullmodel.json"), warn = FALSE), collapse = "\n")
getnum <- function(key) {
  m <- regmatches(js, regexpr(paste0('"', key, '"\\s*:\\s*\\[[^]]*\\]'), js))
  as.numeric(strsplit(gsub('.*\\[|\\]', '', m), ",")[[1]])
}
getstr <- function(key) sub(paste0('.*"', key, '"\\s*:\\s*"([^"]*)".*'), "\\1", js)
getbool <- function(key) grepl(paste0('"', key, '"\\s*:\\s*true'), js)
ids <- regmatches(js, regexpr('"sampleIDs"\\s*:\\s*\\[[^]]*\\]', js))
ids <- gsub('"', '', strsplit(gsub('.*\\[|\\]', '', ids), ",\\s*")[[1]])

theta <- getnum("theta"); alpha <- getnum("alpha"); traitType <- getstr("traitType")
stopifnot(!getbool("loco"))      # LOCO models are not converted
sparse <- getbool("flagSparseGRM")  # sparse-GRM model: R step 2 also needs --sparseGRMFile / --sparseGRMSampleIDFile

mu  <- rarma(file.path(IND, "mu.arma"));  res <- rarma(file.path(IND, "res.arma"))
y   <- as.vector(rarma(file.path(IND, "y.arma")))
X   <- rarma(file.path(IND, "X.arma"));   off <- rarma(file.path(IND, "offset.arma"))
V   <- as.vector(rarma(file.path(IND, "V.arma")))
N <- length(y); p <- ncol(X)
stopifnot(length(ids) == N, nrow(mu) == N, nrow(res) == N, nrow(X) == N, length(V) == N)
colnames(X) <- c("minus1", paste0("x", seq_len(p - 1)))[seq_len(p)]
names(y) <- as.character(seq_len(N))

obj.noK <- list(
  XV         = rarma(file.path(IND, "XV.arma")),          # p x N
  XVX        = rarma(file.path(IND, "XVX.arma")),         # p x p
  XXVX_inv   = rarma(file.path(IND, "XXVX_inv.arma")),    # N x p
  XVX_inv    = rarma(file.path(IND, "XVX_inv.arma")),     # p x p
  S_a        = as.vector(rarma(file.path(IND, "S_a.arma"))),
  XVX_inv_XV = rarma(file.path(IND, "XVX_inv_XV.arma")),  # N x p
  V          = V)
stopifnot(dim(obj.noK$XV) == c(p, N), dim(obj.noK$XXVX_inv) == c(N, p),
          dim(obj.noK$XVX_inv_XV) == c(N, p), length(obj.noK$S_a) == p)

modglmm <- list(
  theta = theta,
  coefficients = matrix(alpha, ncol = 1),
  linear.predictors = NULL,
  fitted.values = mu,                 # N x 1
  Y = NULL,
  residuals = res,                    # N x 1
  cov = NULL,
  converged = TRUE,
  sampleID = ids,
  obj.noK = obj.noK,
  y = y,
  X = X,
  traitType = traitType,
  isCovariateOffset = TRUE,
  LOCO = FALSE,
  obj.glm.null = NULL,
  offset = off,                       # N x 1
  useSparseGRMtoFitNULL = sparse)
save(modglmm, file = OUT)
cat("arma_to_rda: wrote ", OUT, " (n=", N, " p=", p, " trait=", traitType,
    " tau=", paste(formatC(theta, digits = 15, format = "g"), collapse = ","), ")\n", sep = "")
