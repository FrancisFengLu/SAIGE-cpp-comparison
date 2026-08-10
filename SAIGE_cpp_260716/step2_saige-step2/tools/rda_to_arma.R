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
args <- commandArgs(trailingOnly = TRUE)
RDA <- args[1]; OUTD <- args[2]
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

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
  '  "SPA_Cutoff": 2,',
  '  "impute_method": "mean",',
  '  "flagSparseGRM": false,',
  '  "isFastTest": true,',
  '  "isnoadjCov": false,',
  '  "pval_cutoff_for_fastTest": 0.05,',
  '  "isCondition": false,',
  '  "is_Firth_beta": false,',
  '  "pCutoffforFirth": 0.01,',
  paste0('  "sampleIDs": [',
         paste0('"', as.character(modglmm$sampleID), '"', collapse = ", "), ']'),
  "}"), file.path(OUTD, "nullmodel.json"))
cat("rda_to_arma: wrote ", OUTD, " (n=", N, " p=", p, " trait=", traitType,
    " tau=", paste(FMT(tau), collapse = ","), ")\n", sep = "")
