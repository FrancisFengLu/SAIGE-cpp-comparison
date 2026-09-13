#!/usr/bin/env Rscript
# LOCO companion to rda_to_arma.R.
#
#   rda_to_arma_loco.R <null.rda> <outDir>
#
# Writes the same top-level .arma set as rda_to_arma.R (the genome-wide fit),
# plus one chr<j>/ directory per chromosome that step 1 actually produced a
# LOCO fit for, exactly as LOCO_FORMAT.md specifies:
#
#   chr<j>/{mu,res,V,offset,XV,XVX,XVX_inv,XVX_inv_XV,XXVX_inv,S_a}.arma
#
# and sets loco/loco_chroms in nullmodel.json.
#
# WHY IT EXISTS: the parity harness (tests/parity/parity.py) compares R and C++
# step 2 on the *same* step-1 fit.  Without this converter the LOCO path cannot
# be compared at all, because the only alternative -- fitting step 1 twice --
# folds AI-REML's own run-to-run spread into every p-value.
#
# WHAT R ACTUALLY SWAPS PER CHROMOSOME (readInGLMM.R:92-101):
#   fitted.values, residuals, obj.noK          -- always
#   offset                                     -- ONLY when is_Firth_beta AND
#                                                 LOCOResult[[j]]$offset exists
# The offset is read by exactly one thing, the Firth fit (C++
# saige_test.cpp:1170 -> fast_logistf_fit_simple), so under is_Firth_beta=FALSE
# it never reaches a number.  That asymmetry cannot be expressed in a file set:
# the C++ model directory holds ONE offset per chromosome, while R picks
# between two depending on a step-2 flag.  So the default here is the
# genome-wide offset, which is what R uses on every non-Firth run and on every
# run against a SAIGE <= 1.3.3 model (those .rda files carry no per-chromosome
# offset at all).  Pass --loco-offset to write LOCOResult[[j]]$offset instead,
# which is what R uses when is_Firth_beta=TRUE and the model has one.
args <- commandArgs(trailingOnly = TRUE)
USE_LOCO_OFFSET <- "--loco-offset" %in% args
args <- args[args != "--loco-offset"]
if (length(args) < 2)
  stop("usage: rda_to_arma_loco.R [--loco-offset] <null.rda> <outDir>")
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
if (!isTRUE(modglmm$LOCO))
  stop("this .rda was not fit with LOCO=TRUE; use rda_to_arma.R instead")

tau <- modglmm$theta
traitType <- modglmm$traitType
y  <- as.vector(modglmm$y)
X  <- as.matrix(modglmm$X)
N <- length(y); p <- ncol(X)
off_global <- as.vector(modglmm$offset)
if (length(off_global) != N) off_global <- rep(0, N)

# ---- one fit (genome-wide or one chromosome) -> the 10-file set -------------
writeFit <- function(dir, mu, res, objnoK, off) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  mu <- as.vector(mu); res <- as.vector(res)
  V <- if (traitType == "binary") mu * (1 - mu) else rep(1.0 / tau[1], N)
  if (!is.null(objnoK$V)) {
    d <- max(abs(as.vector(objnoK$V) - V))
    if (d > 1e-10)
      cat("  NOTE: obj.noK$V differs from the recomputed V by ", d,
          " in ", dir, "; using the recomputed V\n", sep = "")
  }
  if (is.null(objnoK)) {
    XV <- t(X * V); XVX <- t(X) %*% t(XV); XVX_inv <- solve(XVX)
    XXVX_inv <- X %*% XVX_inv; XVX_inv_XV <- XXVX_inv * V
    S_a <- as.vector(colSums(X * res))
  } else {
    XV <- objnoK$XV; XVX <- objnoK$XVX
    XVX_inv <- if (is.null(objnoK$XVX_inv)) solve(objnoK$XVX) else objnoK$XVX_inv
    XXVX_inv <- objnoK$XXVX_inv
    XVX_inv_XV <- objnoK$XVX_inv_XV
    S_a <- objnoK$S_a
  }
  wvec(file.path(dir, "mu.arma"),  mu)
  wvec(file.path(dir, "res.arma"), res)
  wvec(file.path(dir, "V.arma"),   V)
  wvec(file.path(dir, "offset.arma"), off)
  wvec(file.path(dir, "S_a.arma"), S_a)
  wmat(file.path(dir, "XV.arma"),         XV)
  wmat(file.path(dir, "XVX.arma"),        XVX)
  wmat(file.path(dir, "XVX_inv.arma"),    XVX_inv)
  wmat(file.path(dir, "XXVX_inv.arma"),   XXVX_inv)
  wmat(file.path(dir, "XVX_inv_XV.arma"), XVX_inv_XV)
}

# top level: the genome-wide fit + the chromosome-invariant X and y
writeFit(OUTD, modglmm$fitted.values, modglmm$residuals, modglmm$obj.noK,
         off_global)
wvec(file.path(OUTD, "y.arma"), y)
wmat(file.path(OUTD, "X.arma"), X)

# per chromosome
loco_chroms <- c()
for (j in seq_along(modglmm$LOCOResult)) {
  L <- modglmm$LOCOResult[[j]]
  if (is.null(L)) next
  if (!is.null(L$isLOCO) && !isTRUE(L$isLOCO)) next
  if (is.null(L$fitted.values) || is.null(L$obj.noK)) next
  off_j <- off_global
  if (USE_LOCO_OFFSET) {
    if (!is.null(L$offset)) {
      off_j <- as.vector(L$offset)
    } else if (length(loco_chroms) == 0) {
      cat("NOTE: --loco-offset asked for, but LOCOResult carries no",
          "per-chromosome offset (SAIGE <= 1.3.3); writing the genome-wide",
          "offset.\n")
    }
  }
  writeFit(file.path(OUTD, paste0("chr", j)),
           L$fitted.values, L$residuals, L$obj.noK, off_j)
  loco_chroms <- c(loco_chroms, j)
}
if (length(loco_chroms) == 0)
  stop("LOCO=TRUE but no usable LOCOResult entries were found")

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
  '  "loco": true,',
  '  "lowmem_loco": false,',
  paste0('  "loco_chroms": [', paste(loco_chroms, collapse = ", "), '],'),
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
cat("rda_to_arma_loco: wrote ", OUTD, " (n=", N, " p=", p, " trait=", traitType,
    " tau=", paste(FMT(tau), collapse = ","),
    " loco_chroms=", paste(loco_chroms, collapse = ","), ")\n", sep = "")
