#!/usr/bin/env Rscript
# Convert a non-LOCO SAIGE .rda null model into the .arma + nullmodel.json
# directory that the C++ saige-step2 loads, so that the R and C++ step-2 runs
# use the *identical* null model and any p-value difference is attributable to
# step 2 alone.
#
#   rda_to_arma.R <null.rda> <outDir>
#
# Derived from test_data/Step_2_Feb_11/test/R/convert_rda_to_arma.R
# (hard-coded paths removed, LOCO branch removed).
#   rda_to_arma.R <null.rda> <outDir> [<sparseGRM.mtx> <sparseGRM.sampleIDs.txt>]
# With the two sparse-GRM files (the ones R's step 1 / step 2 take) the GRM K
# is written as sparseGRM_locationMat.arma / sparseGRM_valueVec.arma --
# reindexed to the model's samples, both triangles, as this project's step 1
# writes it (null_model_loader.cpp scales it by tau1 and adds the diagonal) --
# and nullmodel.json says flagSparseGRM true.
args <- commandArgs(trailingOnly = TRUE)
RDA <- args[1]; OUTD <- args[2]
MTX <- if (length(args) >= 4) args[3] else ""
IDS <- if (length(args) >= 4) args[4] else ""
# A fifth argument "stored" writes only the triangle the .mtx stores -- which is
# what R's own step 2 hands its C++ (setSparseSigma_new passes @i / @j / @x of
# the symmetric Tsparse matrix straight to sp_mat, so R's Sigma is lower-
# triangular plus the diagonal). For reproducing R, not for analysis.
TRI <- if (length(args) >= 5) args[5] == "stored" else FALSE
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
wumat <- function(path, m) {   # 2 x nnz, arma::umat (64-bit words)
  con <- file(path, "wb"); on.exit(close(con))
  writeLines(c("ARMA_MAT_BIN_IU008", paste(nrow(m), ncol(m))), con)
  v <- as.numeric(m)
  # 64-bit unsigned little-endian: low and high 32-bit halves
  lo <- v %% 4294967296; hi <- (v - lo) / 4294967296
  raw <- as.vector(rbind(lo, hi))
  writeBin(as.integer(ifelse(raw >= 2147483648, raw - 4294967296, raw)), con, size = 4, endian = "little")
}

wvec <- function(path, v) {
  v <- as.double(v); con <- file(path, "wb"); on.exit(close(con))
  writeLines(c("ARMA_MAT_BIN_FN008", paste(length(v), 1)), con)
  writeBin(v, con, size = 8, endian = "little")
}
wmat <- function(path, m) {
  m <- as.matrix(m); storage.mode(m) <- "double"
  con <- file(path, "wb"); on.exit(close(con))
  writeLines(c("ARMA_MAT_BIN_FN008", paste(nrow(m), ncol(m))), con)
  writeBin(as.double(m), con, size = 8, endian = "little")
}

env <- new.env(); load(RDA, envir = env)
modglmm <- NULL
for (nm in ls(env)) {
  o <- get(nm, envir = env)
  if (is.list(o) && any(c("theta", "fitted.values", "obj.noK") %in% names(o))) {
    modglmm <- o; break
  }
}
stopifnot(!is.null(modglmm))
if (isTRUE(modglmm$LOCO)) stop("rda_to_arma.R only handles non-LOCO models")

tau <- modglmm$theta
traitType <- modglmm$traitType
y  <- as.vector(modglmm$y)
X  <- as.matrix(modglmm$X)
mu <- as.vector(modglmm$fitted.values)
res <- as.vector(modglmm$residuals)
off <- as.vector(modglmm$offset)
if (length(off) != length(y)) off <- rep(0, length(y))
N <- length(y); p <- ncol(X)

V <- if (traitType == "binary") mu * (1 - mu) else rep(1.0 / tau[1], N)

o <- modglmm$obj.noK
if (is.null(o)) {
  XV <- t(X * V); XVX <- t(X) %*% t(XV); XVX_inv <- solve(XVX)
  XXVX_inv <- X %*% XVX_inv; XVX_inv_XV <- XXVX_inv * V
  S_a <- as.vector(colSums(X * res))
} else {
  XV <- o$XV; XVX <- o$XVX
  XVX_inv <- if (is.null(o$XVX_inv)) solve(o$XVX) else o$XVX_inv
  XXVX_inv <- o$XXVX_inv; XVX_inv_XV <- o$XVX_inv_XV; S_a <- o$S_a
}

wvec(file.path(OUTD, "mu.arma"),  mu)
wvec(file.path(OUTD, "res.arma"), res)
wvec(file.path(OUTD, "y.arma"),   y)
wvec(file.path(OUTD, "V.arma"),   V)
wvec(file.path(OUTD, "S_a.arma"), S_a)
wvec(file.path(OUTD, "offset.arma"), off)
wmat(file.path(OUTD, "X.arma"),          X)
wmat(file.path(OUTD, "XV.arma"),         XV)
wmat(file.path(OUTD, "XVX.arma"),        XVX)
wmat(file.path(OUTD, "XVX_inv.arma"),    XVX_inv)
wmat(file.path(OUTD, "XXVX_inv.arma"),   XXVX_inv)
wmat(file.path(OUTD, "XVX_inv_XV.arma"), XVX_inv_XV)

hasSparse <- FALSE
if (nzchar(MTX)) {
  suppressPackageStartupMessages(library(Matrix))
  K <- Matrix::readMM(MTX)
  ids <- readLines(IDS)
  stopifnot(length(ids) == nrow(K))
  pos <- match(as.character(modglmm$sampleID), ids)   # model sample k -> GRM row
  stopifnot(!any(is.na(pos)))
  Ks <- as(K[pos, pos], "TsparseMatrix")             # symmetric storage: one triangle
  i <- Ks@i; j <- Ks@j; x <- Ks@x
  off <- if (TRI) rep(FALSE, length(i)) else (i != j)
  loc <- rbind(c(i, j[off]), c(j, i[off]))            # both triangles, 0-based (TRI: as stored)
  val <- c(x, x[off])
  wumat(file.path(OUTD, "sparseGRM_locationMat.arma"), loc)
  wvec(file.path(OUTD, "sparseGRM_valueVec.arma"), val)
  hasSparse <- TRUE
  cat("rda_to_arma: sparse GRM ", nrow(K), " samples, ", length(val), " stored entries after reindexing\n", sep = "")
}

alpha <- modglmm$coefficients
if (is.null(alpha)) alpha <- rep(0.0, p)
FMT <- function(x) formatC(as.numeric(x), digits = 15, format = "g")
writeLines(c(
  "{",
  paste0('  "trait": "', traitType, '",'),
  paste0('  "traitType": "', traitType, '",'),
  paste0('  "n": ', N, ','),
  paste0('  "p": ', p, ','),
  paste0('  "theta": [', paste(FMT(tau), collapse = ", "), '],'),
  paste0('  "tau": [', paste(FMT(tau), collapse = ", "), '],'),
  paste0('  "alpha": [', paste(FMT(alpha), collapse = ", "), '],'),
  '  "loco": false,',
  '  "lowmem_loco": false,',
  '  "loco_chroms": [],',
  # R keeps these as step-2 options; they are recorded here with R's step-2
  # defaults and saige-step2 does not read them (it has its own config keys,
  # with the same defaults).
  '  "SPA_Cutoff": 2,',
  '  "impute_method": "best_guess",',
  paste0('  "flagSparseGRM": ', if (hasSparse) 'true' else 'false', ','),
  '  "isFastTest": false,',
  '  "isnoadjCov": true,',
  '  "pval_cutoff_for_fastTest": 0.05,',
  '  "isCondition": false,',
  '  "is_Firth_beta": false,',
  '  "pCutoffforFirth": 0.01,',
  paste0('  "sampleIDs": [',
         paste0('"', as.character(modglmm$sampleID), '"', collapse = ", "), ']'),
  "}"), file.path(OUTD, "nullmodel.json"))
cat("rda_to_arma: wrote ", OUTD, " (n=", N, " p=", p, " trait=", traitType,
    " tau=", paste(FMT(tau), collapse = ","), ")\n", sep = "")
