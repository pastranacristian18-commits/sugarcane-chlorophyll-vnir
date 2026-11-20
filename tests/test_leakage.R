# ==============================================================================
# tests/test_leakage.R -- invariants that must hold for every validation scheme
# Run:  Rscript tests/test_leakage.R      (stops with an error if any check fails)
# ==============================================================================
source("config.R"); source("R/functions.R")
D <- load_soca(verbose = FALSE); meta <- D$meta; n <- nrow(meta)
chk <- function(ok, msg) { if (!isTRUE(ok)) stop("FAILED: ", msg); cat("ok  ", msg, "\n") }

for (s in c("random", "lopo", "lobo", "loco", "strict")) {
  sp <- make_splits(meta, s)
  te <- unlist(lapply(sp, `[[`, "test"))
  chk(!anyDuplicated(te) && length(te) == n, sprintf("[%s] every leaf is validated exactly once", s))
  chk(all(sapply(sp, function(z) length(intersect(z$train, z$test)) == 0)), sprintf("[%s] train and test never overlap", s))
}
sp <- make_splits(meta, "lopo");   chk(all(sapply(sp, function(z) length(intersect(meta$plot[z$train],  meta$plot[z$test])) == 0)), "[lopo] no plot is in both train and test")
sp <- make_splits(meta, "lobo");   chk(all(sapply(sp, function(z) length(intersect(meta$block[z$train], meta$block[z$test])) == 0)), "[lobo] no block is in both train and test")
sp <- make_splits(meta, "loco");   chk(all(sapply(sp, function(z) length(intersect(meta$campaign[z$train], meta$campaign[z$test])) == 0)), "[loco] no campaign is in both train and test")
sp <- make_splits(meta, "strict"); chk(all(sapply(sp, function(z) length(intersect(meta$plot[z$train], meta$plot[z$test])) == 0 &&
                                                                 length(intersect(meta$campaign[z$train], meta$campaign[z$test])) == 0)), "[strict] neither plot nor campaign shared")

ref <- colMeans(D$X); x <- D$X[1:5, ]
chk(max(abs(msc_apply(x, ref) - t(sapply(1:5, function(i) { b <- coef(lm(x[i, ] ~ ref)); (x[i, ] - b[1]) / b[2] })))) < 1e-9, "MSC implementation equals a per-spectrum linear fit")
o <- meta$Chla; p <- o + rnorm(n, 0, 0.3); m <- metrics(o, p)
chk(abs(m$R2 - (1 - (n / (n - 1)) / m$RPD^2)) < 1e-8, "R2 (determination) and RPD carry the same information: R2 = 1 - (n/(n-1))/RPD^2")
chk(metrics(o, rep(mean(o), n))$R2 < 1e-12, "a constant mean predictor has R2 <= 0 (r2 is undefined/0)")
Q <- epo_Q(D$X[meta$campaign != "Dia 1", ], meta$campaign[meta$campaign != "Dia 1"], 3)
chk(ncol(Q) == 3 && max(abs(crossprod(Q) - diag(3))) < 1e-8, "EPO projector is orthonormal and estimated from the given rows only")
cat("\nAll leakage / consistency checks passed.\n")
