

#' @title metamed
#' @description Estimation of a summary mediation proportion and summary total
#' effect
#'
#'
#' @param te total effect estimates from each individual study
#' @param te_lb lower bound of total effect confidence intervals
#' @param te_ub upper bound of total effect confidence intervals
#' @param de direct effect estimates from each individual study
#' @param de_lb lower bound of direct effect confidence intervals
#' @param de_ub upper bound of total effect confidence intervals
#' @param corr correlation between direct
#' @param author label of each individual study
#' @param serv serving size of each individual study
#' @param stdserv serving size of pooled total effect
#' @param rr relative risk function, either "exponential" or identity"
#' @param prec precision of pooled estimates
#' @param weight.prec precision of weights on individual studies
#'
#' @return list of data frame of meta-analysis of the total effects, pooled
#' mediation proportion and total effect results, and heterogeneity statistics
#' @import stats
#'
#'
#' @examples
#' S <- 30
#' ST <- 10
#' SB <- 10
#' SD <- 10
#' set.seed(1)
#' beta.T <- rnorm(ST + SB, mean = 0.3, sd = 0.05)
#' set.seed(1)
#' beta.D <- rnorm(SD + SB, mean = 0.2, sd = 0.05)
#' set.seed(1)
#' var.beta.T <- rnorm(ST + SB, mean = 0.05, sd = 0.01)
#' set.seed(1)
#' var.beta.D <- rnorm(SD + SB, mean = 0.05, sd = 0.01)
#' df <- data.frame(beta.T = c(beta.T, rep(NA, SD)),
#'                  beta.D = c(rep(NA, ST), beta.D),
#'                  var.beta.T = c(var.beta.T, rep(NA, SD)),
#'                  var.beta.D = c(rep(NA, ST), var.beta.D))
#' res <- metamed(te = df$beta.T,
#'                te_lb = df$beta.T - qnorm(0.975) * sqrt(df$var.beta.T),
#'                de = df$beta.D,
#'                de_lb = df$beta.D - qnorm(0.975) * sqrt(df$var.beta.D),
#'                author = paste("Study", seq(1:S)),
#'                serv = 1,
#'                stdserv = 1,
#'                rr = "identity")
metamed <- function(te,
                    te_lb,
                    te_ub = NULL,
                    de,
                    de_lb,
                    de_ub = NULL,
                    corr = 0.9855538,
                    author,
                    serv = 1,
                    stdserv = 1,
                    rr = "exp",
                    prec = 3,
                    weight.prec = 2){

  # argument checks
  if (!rr %in% c("exp", "identity"))
    stop("Relative risk function form must be specified as either \"exp\" or
         \"identity\"")

  # get indices from type of study
  sb_ind <- !is.na(te) & !is.na(de)
  sd_ind <- is.na(te) & !is.na(de)
  st_ind <- !is.na(te) & is.na(de)

  if (rr == "exp"){
    te <- log(te)
    te_lb <- log(te_lb)
    de <- log(de)
    de_lb <- log(de_lb)
    if (!is.null(te_ub)) te_ub <- log(te_ub)
    if (!is.null(de_ub)) de_ub <- log(de_ub)
  }

  # standardize serving
  te <- te * stdserv / serv
  te_lb <- te_lb * stdserv / serv
  de <- de * stdserv / serv
  de_lb <- de_lb * stdserv / serv
  if (!is.null(te_ub)) te_ub <- te_ub * stdserv / serv
  if (!is.null(de_ub)) de_ub <- de_ub * stdserv / serv

  # get variances
  if (is.null(te_ub)){
    te_ub <- te + abs(te-te_lb)
    te_var <- ((te-te_lb) / qnorm(0.975))^2
  } else{
    te_var <- (((te_ub - te_lb)/2) / qnorm(0.975))^2
  }
  if (is.null(de_ub)){
    de_ub <- de + abs(de-de_lb)
    de_var <- ((de-de_lb) / qnorm(0.975))^2
  } else{
    de_var <- (((de_ub - de_lb)/2) / qnorm(0.975))^2
  }

  # pooled mediation proportion
  de.bar <- sum(de[sb_ind] / de_var[sb_ind]) / sum(1 / de_var[sb_ind])
  te.bar <- sum(te[sb_ind] / te_var[sb_ind]) / sum(1 / te_var[sb_ind])
  debar_var <- 1 / sum(1 / de_var[sb_ind])
  tebar_var <- 1 / sum(1 / te_var[sb_ind])
  cov <- corr * sqrt(debar_var) * sqrt(tebar_var)
  pbar <- 1 - de.bar / te.bar
  pbar_var <- (de.bar^2) * tebar_var / (te.bar^4) +
    debar_var / (te.bar^2) - 2 * cov * de.bar / (te.bar^3)
  pbar_ci <- c(pbar - qnorm(0.975) * sqrt(pbar_var),
               pbar + qnorm(0.975) * sqrt(pbar_var))

  # heterogeneity measures for MP
  p <- 1 - de/te
  cov <- corr * sqrt(de_var) * sqrt(te_var)
  p_var <- (de^2)*te_var/(te^4) + de_var/(te^2) - 2*de*cov/(te^3)
  p_w.fix <- 1 / p_var
  Q_p <- sum(p_w.fix * p^2, na.rm = TRUE) - sum(p_w.fix * p, na.rm = TRUE)^2 /
    sum(p_w.fix, na.rm = T)
  df_p <- sum(sb_ind)
  C_p <- sum(p_w.fix, na.rm = TRUE) - sum(p_w.fix^2, na.rm = TRUE) /
    sum(p_w.fix, na.rm = TRUE)
  tausq_p <- max(c((Q_p - df_p ) / C_p, 0))
  Isq_p <- (Q_p - df_p) / Q_p
  #print(paste("Isq stat:", Isq_p))

  # MP random effects
  p_w.rand <- 1 / (p_var + tausq_p)
  pbar.rand <- sum(p * p_w.rand, na.rm = T) / sum(p_w.rand, na.rm = T)
  pbar_var.rand <- 1 / sum(p_w.rand, na.rm = T)
  pbar_ci.rand <- c(pbar.rand - qnorm(0.975) * sqrt(pbar_var.rand),
                    pbar.rand + qnorm(0.975) * sqrt(pbar_var.rand))

  CVb_p <- sqrt(tausq_p) / pbar.rand

  # back-calculated total effect
  te[sd_ind] <- de[sd_ind] / (1 - pbar)
  te_var[sd_ind] <- pbar_var * de[sd_ind]^2 / (1 - pbar)^4 +
    de_var[sd_ind] / (1 - pbar)^2
  te_lb[sd_ind] <- te[sd_ind] - qnorm(0.975) * sqrt(te_var[sd_ind])
  te_ub[sd_ind] <- te[sd_ind] + qnorm(0.975) * sqrt(te_var[sd_ind])
  te_w.fix <- 1 / te_var
  #tebcbar <- sum(te[sd_ind] / te_var[sd_ind]) / sum(1 / te_var[sd_ind])

  # TE fixed effects
  te.fix <- sum(te * te_w.fix) / sum(te_w.fix)
  te_var.fix <- 1 / sum(te_w.fix)
  te_ci.fix <- c(te.fix - qnorm(0.975) * sqrt(te_var.fix),
                 te.fix + qnorm(0.975) * sqrt(te_var.fix))

  # heterogeneity stats for TE
  Q_te <- sum(te_w.fix * te^2) - sum(te_w.fix * te)^2 / sum(te_w.fix)
  df_te <- length(te)
  C_te <- sum(te_w.fix) - sum(te_w.fix^2)/sum(te_w.fix)
  tausq_te <- max(c((Q_te - df_te ) / C_te, 0))
  #print(paste("p-val:", 1 - pchisq(Q_te, df_te-1)))
  Isq_te <- (Q_te - df_te) / Q_te
  #print(paste("Isq stat:", Isq_te))

  # TE random effects
  te_w.rand <- 1 / (te_var + tausq_te)
  te.rand <- sum(te * te_w.rand) / sum(te_w.rand)
  te_var.rand <- 1 / sum(te_w.rand)
  te_lb.rand <- te.rand - qnorm(0.975) * sqrt(te_var.rand)
  te_ub.rand <- te.rand + qnorm(0.975) * sqrt(te_var.rand)
  te_ci.rand <- c(te_lb.rand, te_ub.rand)

  CVb_te <- sqrt(tausq_te) / te.rand
  #print(paste("CVb stat:", CVb_te))

  # transform back
  if (rr == "exp"){
    te <- exp(te)
    te_lb <- exp(te_lb)
    te_ub <- exp(te_ub)
    te.fix <- exp(te.fix)
    te.rand <- exp(te.rand)
    te_ci.fix <- exp(te_ci.fix)

    te_ci.rand <- exp(te_ci.rand)
  }

  # rounding
  te <- sprintf(paste0("%.", prec, "f"), te)
  te_ci <- paste0("(", sprintf(paste0("%.", prec, "f"), te_lb), ", ",
                  sprintf(paste0("%.", prec, "f"), te_ub), ")")
  te_w.fix <- sprintf(paste0("%.", weight.prec, "f"), te_w.fix * 100 / sum(te_w.fix))
  te_w.rand <- sprintf(paste0("%.", weight.prec, "f"), te_w.rand * 100 / sum(te_w.rand))

  tedf <- data.frame(author, te, te_ci, te_w.fix, te_w.rand)
  names(tedf) <- c("Study", "TE", "95% CI", "Weight (fixed)", "Weight (random)")

  return(list("TE" = tedf,
         "pbar" = pbar,
         "pbar.ci" = pbar_ci,
         "pbar.tau2" = tausq_p,
         "te.fix" = te.fix,
         "te.ci.fix" = te_ci.fix,
         "te.rand" = te.rand,
         "te.ci.rand" = te_ci.rand,
         "te.tau2" = tausq_te)
  )
}
