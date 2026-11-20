# ==============================================================================
# R/functions.R -- data loading, preprocessing, models, validation, metrics
# ==============================================================================
suppressPackageStartupMessages({ library(pls) })

DDC_MAP <- c("Dia 1" = 60, "Dia 2" = 90, "Dia 3" = 120,
             "Dia 4" = 150, "Dia 5" = 180, "Dia 6" = 210)
PIGMENTS  <- c("Chla", "Chlb")
PIG_LABEL <- c(Chla = "Chl a", Chlb = "Chl b")

# ------------------------------------------------------------------------------
# Data
# ------------------------------------------------------------------------------
#' Load chemistry + spectra, remove duplicates, build metadata.
#' Rows keep the order of the raw file (needed to reproduce the manuscript split).
load_soca <- function(cfg = CFG, verbose = TRUE) {
  d <- read.table(cfg$data_txt, header = TRUE, sep = "\t", dec = ".",
                  check.names = FALSE)
  spec_cols <- grep("^[0-9]", colnames(d))
  key <- paste(d$Dia, d$Unidad)
  dup <- which(duplicated(key))
  if (length(dup)) {
    if (verbose) message("[data] dropping duplicated key(s): ",
                         paste(key[dup], collapse = ", "), " (row ", paste(dup, collapse = ","), ")")
    d <- d[-dup, ]; key <- key[-dup]
  }
  if (isTRUE(cfg$drop_p6s1_150ddc)) {
    drop <- which(d$Dia == "Dia 4" & d$Unidad == "P6S1")
    if (verbose && length(drop)) message("[data] dropping Day 4 / P6S1 (150 DDC) -- reason must be documented")
    if (length(drop)) { d <- d[-drop, ]; key <- key[-drop] }
  }
  meta <- data.frame(
    id       = key,
    campaign = factor(d$Dia, levels = names(DDC_MAP)),
    DDC      = unname(DDC_MAP[d$Dia]),
    plot     = as.integer(sub("^P(\\d+)S\\d+$", "\\1", d$Unidad)),   # experimental unit (1-20)
    leaf     = as.integer(sub("^P\\d+S(\\d+)$", "\\1", d$Unidad)),   # sub-sample within plot (1-4)
    block    = d$Bloque_Num,
    trt      = d$Tratamiento,
    Chla     = d$Chla,
    Chlb     = d$Chlb,
    stringsAsFactors = FALSE
  )
  X <- as.matrix(d[, spec_cols]); rownames(X) <- key
  if (verbose) message("[data] n = ", nrow(meta), " leaves | ", ncol(X), " bands | ",
                       length(unique(meta$plot)), " plots")
  list(meta = meta, X = X, wl = as.numeric(colnames(X)))
}

# ------------------------------------------------------------------------------
# Metrics
#   r2        squared Pearson correlation (what the original code called "R2")
#   R2        coefficient of determination, 1 - SSE/SST (can be negative)
#   RPD       SD(obs) / RMSE
#   RPIQ      IQR(obs) / RMSE
#   R2_within 1 - SSE / SS around the CAMPAIGN means (does the model beat "knowing the campaign"?)
# ------------------------------------------------------------------------------
metrics <- function(obs, pred, campaign = NULL) {
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  sse <- sum((obs - pred)^2); rmse <- sqrt(mean((obs - pred)^2))
  out <- data.frame(
    n = length(obs),
    r2 = suppressWarnings(cor(obs, pred))^2,
    R2 = 1 - sse / sum((obs - mean(obs))^2),
    RMSE = rmse, RPD = sd(obs) / rmse, RPIQ = IQR(obs) / rmse, Bias = mean(pred - obs))
  if (!is.null(campaign)) {
    cm <- ave(obs, campaign[ok]); out$R2_within <- 1 - sse / sum((obs - cm)^2)
  }
  out
}

# ------------------------------------------------------------------------------
# Spectral preprocessing that must be re-estimated inside each training fold
# ------------------------------------------------------------------------------
msc_apply <- function(X, ref) {                       # equals MATLAB polyfit(ref, x, 1)
  rc <- ref - mean(ref); rm <- rowMeans(X)
  b  <- as.vector((X %*% rc) - rm * sum(rc)) / sum(rc^2)
  a  <- rm - b * mean(ref)
  (X - a) / b
}
pp_msc <- function(Xtr, Xte, mtr, ytr) {
  ref <- colMeans(Xtr)
  list(train = msc_apply(Xtr, ref), test = msc_apply(Xte, ref))
}

#' EPO projector from group (campaign) means of the TRAINING data only.
epo_Q <- function(Xtr, groups, k) {
  g  <- factor(groups)
  mu <- do.call(rbind, lapply(levels(g), function(l) colMeans(Xtr[g == l, , drop = FALSE])))
  ng <- nrow(mu); D <- NULL
  for (i in 1:(ng - 1)) for (j in (i + 1):ng) D <- rbind(D, mu[i, ] - mu[j, ])
  k <- min(k, ng - 1)
  svd(D)$v[, seq_len(k), drop = FALSE]
}
epo_apply <- function(X, Q) X - (X %*% Q) %*% t(Q)
pp_epo <- function(k) function(Xtr, Xte, mtr, ytr) {
  Q <- epo_Q(Xtr, mtr$campaign, k)
  list(train = epo_apply(Xtr, Q), test = epo_apply(Xte, Q))
}

# ------------------------------------------------------------------------------
# PLSR with latent-variable selection by INNER cross-validation grouped by plot
# ------------------------------------------------------------------------------
fit_pls <- function(Xtr, ytr, groups = NULL, max_lv = CFG$max_lv,
                    inner_k = CFG$inner_k, seed = 1) {
  max_lv <- min(max_lv, nrow(Xtr) - 2)
  dtr <- data.frame(y = ytr); dtr$X <- I(Xtr)
  set.seed(seed)
  seg <- 10
  if (!is.null(groups)) {
    ug <- sort(unique(groups)); kk <- min(inner_k, length(ug))
    a  <- sample(rep(seq_len(kk), length.out = length(ug)))
    seg <- unname(split(seq_along(ytr), a[match(groups, ug)]))
  }
  m  <- pls::plsr(y ~ X, ncomp = max_lv, data = dtr, method = "oscorespls",
                  validation = "CV", segments = seg)
  cv <- pls::RMSEP(m, estimate = "CV")$val[1, 1, -1]
  list(model = m, lv = which.min(cv), cv_rmse = cv)
}
predict_pls <- function(fit, Xte) {
  as.vector(predict(fit$model, newdata = data.frame(X = I(Xte)), ncomp = fit$lv))
}

#' Generic cross-validation runner. `splits` = list(list(train=, test=), ...)
#' `preprocess(Xtr, Xte, meta_tr, ytr)` is re-estimated on each training set.
cv_pls <- function(X, y, meta, splits, preprocess = NULL, seed = 1, ...) {
  pred <- rep(NA_real_, length(y)); lv <- integer(length(splits)); kk <- rep(NA_real_, length(splits))
  for (i in seq_along(splits)) {
    tr <- splits[[i]]$train; te <- splits[[i]]$test
    Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
    if (!is.null(preprocess)) {
      p <- preprocess(Xtr, Xte, meta[tr, , drop = FALSE], y[tr]); Xtr <- p$train; Xte <- p$test
      if (!is.null(p$k)) kk[i] <- p$k
    }
    fit <- fit_pls(Xtr, y[tr], groups = meta$plot[tr], seed = seed + i, ...)
    pred[te] <- predict_pls(fit, Xte); lv[i] <- fit$lv
  }
  list(pred = pred, lv = lv, k = kk)
}

#' PLSR on CAMPAIGN-CENTRED spectra and responses ("within-campaign" skill).
#'   center_test = "train": test campaign also present in training -> its mean spectrum/response come from training rows
#'   center_test = "self" : test campaign unseen (LOCO) -> its mean SPECTRUM is taken from the (unlabelled) test spectra;
#'                          the response deviation is then evaluated against the campaign's own mean (oracle offset).
#' Returns predicted and observed DEVIATIONS from the campaign mean.
cv_pls_centered <- function(X, y, meta, splits, center_test = c("train", "self"), seed = 1, ...) {
  center_test <- match.arg(center_test); n <- length(y)
  dev_pred <- dev_obs <- rep(NA_real_, n)
  cmean <- function(M, g) { u <- levels(droplevels(g)); r <- t(sapply(u, function(l) colMeans(M[g == l, , drop = FALSE]))); rownames(r) <- u; r }
  for (i in seq_along(splits)) {
    tr <- splits[[i]]$train; te <- splits[[i]]$test
    gtr <- droplevels(meta$campaign[tr]); gte <- as.character(meta$campaign[te])
    mX <- cmean(X[tr, , drop = FALSE], gtr); my <- tapply(y[tr], gtr, mean)
    Xtr <- X[tr, , drop = FALSE] - mX[as.character(gtr), , drop = FALSE]; ytr <- y[tr] - my[as.character(gtr)]
    if (center_test == "train") { Xte <- X[te, , drop = FALSE] - mX[gte, , drop = FALSE]; ote <- y[te] - my[gte]
    } else { Xte <- X[te, , drop = FALSE] - matrix(colMeans(X[te, , drop = FALSE]), length(te), ncol(X), byrow = TRUE); ote <- y[te] - mean(y[te]) }
    fit <- fit_pls(Xtr, as.numeric(ytr), groups = meta$plot[tr], seed = seed + i, ...)
    dev_pred[te] <- predict_pls(fit, Xte); dev_obs[te] <- ote
  }
  list(dev_pred = dev_pred, dev_obs = dev_obs)
}
dev_metrics <- function(dev_obs, dev_pred) {
  ok <- is.finite(dev_obs) & is.finite(dev_pred); o <- dev_obs[ok]; p <- dev_pred[ok]
  data.frame(n = length(o), r2_dev = suppressWarnings(cor(o, p))^2,
             R2_dev = 1 - sum((o - p)^2) / sum(o^2), RMSE_dev = sqrt(mean((o - p)^2)), RPD_dev = sqrt(mean(o^2)) / sqrt(mean((o - p)^2)))
}

#' Nested EPO: k (0..kmax) chosen by an inner leave-one-campaign-out inside the training set only.
pp_epo_nested <- function(kmax = 3, seed = 1) function(Xtr, Xte, mtr, ytr) {
  cs  <- sort(unique(as.character(mtr$campaign)))
  isp <- lapply(cs, function(c) list(train = which(mtr$campaign != c), test = which(mtr$campaign == c)))
  sc  <- sapply(0:kmax, function(k) {
    r <- cv_pls(Xtr, ytr, mtr, isp, preprocess = if (k == 0) NULL else pp_epo(k), seed = seed)
    sqrt(mean((ytr - r$pred)^2)) })
  kbest <- (0:kmax)[which.min(sc)]
  if (kbest == 0) return(list(train = Xtr, test = Xte, k = 0))
  Q <- epo_Q(Xtr, mtr$campaign, kbest)
  list(train = epo_apply(Xtr, Q), test = epo_apply(Xte, Q), k = kbest)
}

# ------------------------------------------------------------------------------
# Validation schemes (each returns a list of train/test index sets)
#   random : random 5-fold stratified by campaign, ignores plots  (naive "80/20")
#   lopo   : leave-one-plot-out   -> new experimental units, same campaigns
#   lobo   : leave-one-block-out  -> new blocks
#   loco   : leave-one-campaign-out -> new time, same plots
#   strict : new campaign AND new plots (train excludes both)
# ------------------------------------------------------------------------------
make_splits <- function(meta, scheme, seed = CFG$seed, k = CFG$k_random) {
  n <- nrow(meta); idx <- seq_len(n)
  from_fold <- function(f) lapply(sort(unique(f)), function(g)
    list(train = idx[f != g], test = idx[f == g], id = g))
  plot_group <- function() {
    set.seed(seed + 1)
    g <- sample(rep(seq_len(k), length.out = length(unique(meta$plot))))
    g[match(meta$plot, sort(unique(meta$plot)))]
  }
  switch(scheme,
    random = { set.seed(seed); f <- integer(n)
               for (cm in levels(droplevels(meta$campaign))) {
                 i <- which(meta$campaign == cm)
                 f[i] <- sample(rep(seq_len(k), length.out = length(i))) }
               from_fold(f) },
    lopo   = from_fold(meta$plot),
    lobo   = from_fold(meta$block),
    loco   = from_fold(as.integer(meta$campaign)),
    strict = { pg <- plot_group(); ca <- as.integer(meta$campaign); out <- list()
               for (c in sort(unique(ca))) for (g in seq_len(k))
                 out[[length(out) + 1]] <- list(train = idx[ca != c & pg != g],
                                                test  = idx[ca == c & pg == g], id = paste(c, g))
               out },
    stop("unknown scheme: ", scheme))
}

#' Reference predictors evaluated on the same splits (no spectra involved).
#'   mean          : mean of the training responses
#'   campaign_mean : mean of the training responses of the SAME campaign ("knowing the campaign")
cv_baseline <- function(y, meta, splits, type = c("mean", "campaign_mean")) {
  type <- match.arg(type); pred <- rep(NA_real_, length(y))
  for (s in splits) {
    tr <- s$train
    for (i in s$test) {
      if (type == "mean") { pred[i] <- mean(y[tr]) } else {
        same <- tr[meta$campaign[tr] == meta$campaign[i]]
        pred[i] <- if (length(same)) mean(y[same]) else mean(y[tr]) }
    }
  }
  pred
}

# ------------------------------------------------------------------------------
# Cluster (plot-level) bootstrap on stored out-of-fold predictions -- no refits.
# Resamples the 20 plots with replacement, keeping all their leaves and campaigns.
# ------------------------------------------------------------------------------
boot_cluster <- function(y, preds, cluster, B = CFG$B_boot, seed = CFG$seed) {
  set.seed(seed); cl <- sort(unique(cluster)); ix <- split(seq_along(y), cluster)
  nm <- names(preds)
  out <- array(NA_real_, c(B, length(preds), 3), list(NULL, nm, c("R2", "RMSE", "RPD")))
  for (b in seq_len(B)) {
    s <- unlist(ix[as.character(sample(cl, length(cl), replace = TRUE))], use.names = FALSE)
    for (j in seq_along(preds)) {
      p <- preds[[j]][s]; o <- y[s]; ok <- is.finite(p)
      rm <- sqrt(mean((o[ok] - p[ok])^2))
      out[b, j, ] <- c(1 - sum((o[ok] - p[ok])^2) / sum((o[ok] - mean(o[ok]))^2), rm, sd(o[ok]) / rm)
    }
  }
  out
}
ci <- function(x, p = c(.025, .975)) quantile(x, p, na.rm = TRUE, names = FALSE)

# ------------------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------------------
VIP <- function(object) {
  SS  <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wn  <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wn, "*")
  t(sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS)))
}
write_table <- function(df, name, cfg = CFG) {
  dir.create(cfg$out_tables, recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(cfg$out_tables, name), row.names = FALSE)
  message("[table] ", name)
}
