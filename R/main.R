
ss <- function(x) sum(x^2)
significance_codes <- function(pvals) {
  symnum(pvals, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
         symbols = c("***", "**", "*", ".", " "))
}

make_twi_string <- function(var, vars) {
  strs <- c()
  for (v in vars) {
    if (var != v) {
      strs <- c(strs, paste(var, v, sep=":"))
    }
  }
  paste(strs, collapse=" + ")
}

#' ANOVA with FO, PQ, and TWI terms.
#'
#' @param formula Complete formula for the model.
#' @param formula Complete formula for the model.
#' @param data Data frame.
#' @param vars Character vector of variable names for ANOVA.
#' @return A data frame with ANOVA statistics.
#'
#' The ANOVA is not fully sequential. The FO and PQ terms are
#' added for the first variable, followed by the TWI terms.
#' The FO and PQ terms for the second variable are added only
#' to the FO and PQ model for the first variable since there
#' is overlap in the TWI terms. Thus the FO and PQ terms for
#' the final variable are added to all other FO and PQ terms
#' but no TWI terms.
#'
#' @export
qanova <- function(formula, data, vars) {
  rhs_vars <- vars
  nvars <- length(rhs_vars)
  tbl <- data.frame(variable=c(rep(rhs_vars, each=2), "residuals"),
                    model=c(rep(c("FO+PQ", "+TWI"), times=length(rhs_vars)), NA),
                    df=0, sum_sq=0, mean_sq=0, Fvalue=0, PrF=0,
                    stringsAsFactors = FALSE)
  ri <- nrow(tbl)

  full_model <- lm(formula, data)
  tbl$df[ri] <- full_model$df.residual
  tbl$sum_sq[ri] <- ss(full_model$residuals)
  tbl$mean_sq[ri] <- tbl$sum_sq[ri] / tbl$df[ri]
  tbl$Fvalue[ri] <- NA
  tbl$PrF[ri] <- NA

  lhs_var <- formula.tools::lhs.vars(formula)
  sst <- var(data[ ,lhs_var]) * (nrow(data)-1)
  base_formula <- paste(lhs_var, "~", "1")
  sse_base <- 0
  df_base <- 1

  for (i in 1:nvars) {
    var <- rhs_vars[i]

    # FO + PQ
    idx <- 2*(i-1) + 1
    base_formula <- paste(base_formula, "+", var, "+", paste("I(", var, "^2)", sep=""))
    model <- lm(base_formula, data)
    tbl$df[idx] <- length(model$coefficients) - df_base

    sse_new <- sst - ss(model$residuals)
    tbl$sum_sq[idx] <- sse_new - sse_base
    tbl$mean_sq[idx] <- tbl$sum_sq[idx] / tbl$df[idx]
    tbl$Fvalue[idx] <- tbl$mean_sq[idx] / tbl$mean_sq[ri]
    tbl$PrF[idx] <- 1 - pf(tbl$Fvalue[idx], tbl$df[idx], tbl$df[ri])
    sse_base <- sse_new
    df_base <- length(model$coefficients)

    # TWI
    idx <- 2*(i-1) + 2
    new_formula <- paste(base_formula, "+", make_twi_string(var, rhs_vars))
    model <- lm(new_formula, data)
    tbl$df[idx] <- length(model$coefficients) - df_base
    tbl$sum_sq[idx] <- sst - ss(model$residuals) - sse_base
    tbl$mean_sq[idx] <- tbl$sum_sq[idx] / tbl$df[idx]
    tbl$Fvalue[idx] <- tbl$mean_sq[idx] / tbl$mean_sq[ri]
    tbl$PrF[idx] <- 1 - pf(tbl$Fvalue[idx], tbl$df[idx], tbl$df[ri])
  }

  tbl$sig <- significance_codes(tbl$PrF)
  tbl
}

