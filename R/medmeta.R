

#' @title metamed
#' @description Estimates a pooled mediation proportion and total effects and
#' outputs results in data frame
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
#' @param pmodel model used for back-calculated mediation proportion, either
#' "fixed" or "random" effects model
#' @param prec precision of pooled estimates
#' @param weight.prec precision of weights on individual studies
#'
#' @return data frames of meta-analysis tables for mediation proportion and
#' total effect, pooled mediation proportion and total effect results, and
#' heterogeneity statistics
#' @export
#' @import stats
#'
#' @examples
medmeta <- function(te,
                    te_lb,
                    te_ub = NULL,
                    de,
                    de_lb,
                    de_ub = NULL,
                    corr = 0.9855538,
                    author,
                    serv,
                    stdserv,
                    rr = "exp",
                    pmodel = "fixed",
                    prec = 3,
                    weight.prec = 2){
  #TODO check args

  # get indices from type of study
  sb_ind <- !is.na(te) & !is.na(de)
  sd_ind <- is.na(te) & !is.na(de)
  st_ind <- is.na(te) & !is.na(de)

  if (rr == "exp"){
    te <- log(te)
    te_lb <- log(te_lb)
    de <- log(de)
    de_lb <- log(de_lb)
    if (!is.null(te_ub)) te_ub <- log(te_ub)
    if (!is.null(de_ub)) de_ub <- log(de_ub)
  }

  # standardize
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

  # individual MP's and CI's (fixed effects)
  p <- 1 - de/te
  cov <- corr * sqrt(de_var) * sqrt(te_var)
  p_var <- (de^2)*te_var/(te^4) + de_var/(te^2) - 2*de*cov/(te^3)
  p_w.fix <- 1 / p_var
  p_lb <- p - qnorm(0.975) * sqrt(p_var)
  p_ub <- p + qnorm(0.975) * sqrt(p_var)

  # MP fixed effects
  pbar.fix <- sum(p * p_w.fix, na.rm = T) / sum(p_w.fix, na.rm = T)
  pbar_var.fix <- 1 / sum(p_w.fix, na.rm = T)
  pbar_lb.fix <- pbar.fix - qnorm(0.975) * sqrt(pbar_var.fix)
  pbar_ub.fix <- pbar.fix + qnorm(0.975) * sqrt(pbar_var.fix)
  pbar_ci.fix <- c(pbar_lb.fix, pbar_ub.fix)

  # heterogeneity stats for MP
  Q_p <- sum(p_w.fix * p^2, na.rm = T) - sum(p_w.fix*p, na.rm = T)^2 /
    sum(p_w.fix, na.rm = T)
  #print(paste("Q_p:", Q_p))
  df_p <- sum(sb_ind)
  #print(paste("p-val:", 1 - pchisq(Q_p, df_p-1)))
  C_p <- sum(p_w.fix, na.rm = T) - sum(p_w.fix^2, na.rm = T) / sum(p_w.fix, na.rm = T)
  tausq_p <- (Q_p - df_p ) / C_p
  Isq_p <- (Q_p - df_p) / Q_p
  #print(paste("Isq stat:", Isq_p))

  # MP random effects
  p_w.rand <- 1 / (p_var + tausq_p)
  pbar.rand <- sum(p * p_w.rand, na.rm = T) / sum(p_w.rand, na.rm = T)
  pbar_var.rand <- 1 / sum(p_w.rand, na.rm = T)
  pbar_lb.rand <- pbar.rand - qnorm(0.975) * sqrt(pbar_var.rand)
  pbar_ub.rand <- pbar.rand + qnorm(0.975) * sqrt(pbar_var.rand)
  pbar_ci.rand <- c(pbar_lb.rand, pbar_ub.rand)

  CVb_p <- sqrt(tausq_p) / pbar.rand
  #print(paste("CVb stat:", CVb_p))

  # rounding
  p <- sprintf(paste0("%.", prec, "f"), p)
  p_ci <- paste0("(", sprintf(paste0("%.", prec, "f"), p_lb[sb_ind]), ", ",
                 sprintf(paste0("%.", prec, "f"), p_ub[sb_ind]), ")")
  p_w.fix <- sprintf(paste0("%.", weight.prec, "f"), p_w.fix[sb_ind] * 100 /
                       sum(p_w.fix[sb_ind]))
  p_w.rand <- sprintf(paste0("%.", weight.prec, "f"), p_w.rand[sb_ind] * 100 /
                        sum(p_w.rand[sb_ind]))

  pdf <- data.frame(author[sb_ind], p[sb_ind], p_ci, p_w.fix,  p_w.rand)
  names(pdf) <- c("Study", "MP", "95% CI", "Weight (fixed)", "Weight (random)")

  # back-calculated total effect
  if (pmodel == "fixed"){
    tebc <- de[sd_ind] / (1 - pbar.fix)
    tebc_var <- pbar_var.fix * (de[sd_ind]^2) / (1-pbar.fix)^4 +
      de_var[sd_ind] / (1-pbar.fix)^2
  } else{
    tebc <- de[sd_ind] / (1 - pbar.rand)
    tebc_var <- pbar_var.rand * (de[sd_ind]^2) / (1-pbar.rand)^4 +
      de_var[sd_ind] / (1-pbar.rand)^2
  }
  te[sd_ind] <- tebc
  te_var[sd_ind] <- tebc_var
  te_lb[sd_ind] <- te[sd_ind] - qnorm(0.975) * sqrt(te_var[sd_ind])
  te_ub[sd_ind] <- te[sd_ind] + qnorm(0.975) * sqrt(te_var[sd_ind])
  te_w.fix <- 1 / te_var

  # TE fixed effects
  te.fix <- sum(te * te_w.fix) / sum(te_w.fix)
  te_var.fix <- 1 / sum(te_w.fix)
  te_lb.fix <- te.fix - qnorm(0.975) * sqrt(te_var.fix)
  te_ub.fix <- te.fix + qnorm(0.975) * sqrt(te_var.fix)
  te_ci.fix <- c(te_lb.fix, te_ub.fix)

  # heterogeneity stats for MP
  Q_te <- sum(te_w.fix * te^2) - sum(te_w.fix * te)^2 / sum(te_w.fix)
  #Q_te <- sum(te_w.fix*c(te[!sd_ind], tebc)^2) - sum(te_w.fix*c(te[!sd_ind], tebc))^2 / sum(te_w.fix)
  df_te <- length(te)
  C_te <- sum(te_w.fix) - sum(te_w.fix^2)/sum(te_w.fix)
  tausq_te <- (Q_te - df_te ) / C_te
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

  ### PLOTTING
  #require(gplots)
  # par(mar = c(5, 6, 3, 2), mgp = c(2.3, 0.6, 0))
  # plotCI(c(pbar.rand, pbar.fix, rev(p[both_ind])), 1:nrow(pdf),
  #        ui = c(pbar_ub.rand, pbar_ub.fix, rev(p_ub[both_ind])),
  #        li = c(pbar_lb.rand, pbar_lb.fix, rev(p_lb[both_ind])),
  #        err = "x", xlab = "Mediation Proportion", ylab = "", axes = F, pch = 20)
  # axis(side = 1, cex.axis = 0.75, tck = -0.015)
  # axis(side = 2, at = 1:nrow(pdf), label = rev(pdf$Study), las = 2,
  #      cex.axis = 0.75, tck = -0.015)


  # max(ceiling(log2(te_all_ub)))
  # plotCI(c(te.rand, te.fix, rev(log(te_all[sind]))), 1:nrow(tedf),
  #        ui = c(te_ub.rand, te_ub.fix, rev(log(te_all_ub[sind]))),
  #        li = c(te_lb.rand, te_lb.fix, rev(log(te_all_lb[sind]))),
  #        err = "x", xlab = "Total Effect (RR)", ylab = "", axes = F, pch = 20,
  #        xlim = c(log(0.5), log(2^max(ceiling(log2(te_all_ub))))))
  # aty <- 2^(-1:max(ceiling(log2(te_all_ub))))
  # axis(side = 1, at = log(aty), labels = aty, cex.axis = 0.75, tck = -0.015)
  # axis(side = 2, at = 1:nrow(tedf), label = rev(tedf$Study), las = 2,
  #      cex.axis = 0.75, tck = -0.015)
  # abline(v = 0, lty = 2)

  return(list("TE" = tedf,
              "MP" = pdf,
              "pbar.fix" = pbar.fix,
              "pbar.ci.fix" = pbar_ci.fix,
              "pbar.rand" = pbar.rand,
              "pbar.ci.rand" = pbar_ci.rand,
              "te.fix" = te.fix,
              "te.ci.fix" = te_ci.fix,
              "te.rand" = te.rand,
              "te.ci.rand" = te_ci.rand
  ))
}

