robust_mahalanobis_test <- function(mat,
                                    mcd.alpha = 0.75,
                                    seed = 1,
                                    method = c("GM14", "HR05")) {

  method <- match.arg(method)

  mat <- as.data.frame(mat)
  numeric_cols <- vapply(mat, is.numeric, logical(1))
  if (!all(numeric_cols)) {
    warning(
      "Dropping non-numeric columns: ",
      paste(names(mat)[!numeric_cols], collapse = ", ")
    )
    mat <- mat[, numeric_cols, drop = FALSE]
  }

  if (ncol(mat) < 2L) {
    stop("Need at least 2 numeric variables for Mahalanobis distance.")
  }

  mat_num <- as.matrix(mat)
  if (anyNA(mat_num) || any(!is.finite(mat_num))) {
    stop("mat contains NA/Inf; clean data before Mahalanobis test.")
  }

  storage.mode(mat_num) <- "double"
  n <- nrow(mat_num)
  p <- ncol(mat_num)

  if (n <= p) {
    stop("Need more observations than variables for robust Mahalanobis estimation.")
  }

  set.seed(seed)
  mcd <- robustbase::covMcd(mat_num, alpha = mcd.alpha)
  D2 <- mahalanobis(mat_num, center = mcd$raw.center, cov = mcd$raw.cov)

  if (!is.numeric(D2) || anyNA(D2) || any(!is.finite(D2))) {
    stop("Computed D2 contains non-finite values.")
  }

  get_pvalue <- function(d2) {
    if (!is.numeric(d2) || length(d2) != 1L || !is.finite(d2)) {
      stop(
        "d2 is not a finite numeric scalar: class=", paste(class(d2), collapse = ","),
        " length=", length(d2),
        " value=", paste(d2, collapse = ",")
      )
    }

    f <- function(a) {
      if (!is.numeric(a) || length(a) != 1L || !is.finite(a)) {
        stop("a is not a finite numeric scalar.")
      }

      res <- CerioliOutlierDetection::hr05CutoffMvnormal(
        n.obs = n,
        p.dim = p,
        mcd.alpha = mcd.alpha,
        signif.alpha = a,
        method = method
      )

      cutoff <- res$cutoff.pred

      if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff)) {
        stop(
          "cutoff is not a finite numeric scalar: class=", paste(class(cutoff), collapse = ","),
          " length=", length(cutoff),
          " value=", paste(cutoff, collapse = ",")
        )
      }

      cutoff - d2
    }

    upper_val <- f(1 - 1e-8)
    lower_val <- f(1e-8)

    if (upper_val >= 0) return(1)
    if (lower_val <= 0) return(1e-8)

    uniroot(f, interval = c(1e-8, 1 - 1e-8))$root
  }

  pvalue <- vapply(D2, get_pvalue, numeric(1))
  qvalue <- p.adjust(pvalue, method = "BH")

  list(
    D2 = D2,
    pvalue = pvalue,
    qvalue = qvalue,
    mcd = mcd,
    method = method,
    n = n,
    p = p,
    mcd.alpha = mcd.alpha
  )
}
