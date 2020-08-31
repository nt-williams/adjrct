
#' Plot Restricted Mean Survival Time Estimates
#'
#' @param x An object of class "rmst" produced by a call to \code{rmst()}
#' @param ... Extra arguments to be passed to plot, not used.
#'
#' @return Object of class \code{ggplot} containing a patchwork step function plot of the
#'  treatment and control arm RMST point estimates and their differences across a sequence
#'  of time horizons.
#'
#' @method plot rmst
#'
#' @importFrom ggplot2 ggplot aes_ geom_step geom_line labs scale_color_manual geom_ribbon scale_fill_manual theme element_blank
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survrct)
#'   veteran <- survival::veteran
#'   veteran$trt <- veteran$trt - 1
#'   veteran$celltype <- factor(veteran$celltype)
#'   surv <- survrct(Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior,
#'                   target = "trt", data = veteran, coarsen = 7, estimator = "tmle")
#'   estm <- rmst(surv, 10:40)
#'   plot(estm)
#'   }
#' }
plot.rmst <- function(x, ...) {
  estm <- all_estimates(x)
  clrs <- c("Treatment" = "blue", "Control" = "red")
  flls <- c("95% Uniform" = "lightblue", "95% Pointwise" = "blue")
  top <- ggplot(estm, aes_(x = ~horizon)) +
    geom_step(aes_(y = ~treatment, color = "Treatment")) +
    geom_step(aes_(y = ~control, color = "Control")) +
    labs(x = NULL,
                  y = "Restricted Mean Survival Time",
                  color = NULL) +
    scale_color_manual(values = clrs) +
    theme(axis.ticks.x.bottom = element_blank(),
          axis.text.x.bottom = element_blank())

  bottom <- ggplot(estm, aes_(x = ~horizon, y = ~theta)) +
    geom_ribbon(aes_(ymin = ~unif.conf.low, ymax = ~unif.conf.high, fill = "95% Uniform"), alpha = 0.4) +
    geom_ribbon(aes_(ymin = ~theta.conf.low, ymax = ~theta.conf.high, fill = "95% Pointwise"), alpha = 0.35) +
    geom_line() +
    scale_fill_manual(values = flls) +
    labs(x = "Time horizon",
                  y = "Difference",
                  fill = NULL)

  wrap_plots(top, bottom, ncol = 1, heights = c(3, 1))
}

#' Plot Survival Probability Estimates
#'
#' @param x An object of class "survprob" produced by a call to \code{survprob()}
#' @param ... Extra arguments to be passed to plot, not used.
#'
#' @return Object of class \code{ggplot} containing a patchwork step function plot of the
#'  treatment and control arm survpval probability point estimates and their differences
#'  across a sequence of time horizons.
#'
#' @method plot survprob
#'
#' @importFrom ggplot2 ggplot aes_ geom_step geom_line labs scale_color_manual geom_ribbon scale_fill_manual theme element_blank
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survrct)
#'   veteran <- survival::veteran
#'   veteran$trt <- veteran$trt - 1
#'   veteran$celltype <- factor(veteran$celltype)
#'   surv <- survrct(Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior,
#'                   target = "trt", data = veteran, coarsen = 7, estimator = "tmle")
#'   estm <- survprob(surv, 10:40)
#'   plot(estm)
#'   }
#' }
plot.survprob <- function(x, ...) {
  estm <- all_estimates(x)
  clrs <- c("Treatment" = "blue", "Control" = "red")
  flls <- c("95% Uniform" = "lightblue", "95% Point-wise" = "blue")
  top <- ggplot(estm, aes_(x = ~horizon)) +
    geom_step(aes_(y = ~treatment, color = "Treatment")) +
    geom_step(aes_(y = ~control, color = "Control")) +
    labs(x = NULL,
                  y = "Survival probability",
                  color = NULL) +
    scale_color_manual(values = clrs) +
    theme(axis.ticks.x.bottom = element_blank(),
                   axis.text.x.bottom = element_blank())

  bottom <- ggplot(estm, aes_(x = ~horizon, y = ~theta)) +
    geom_ribbon(aes_(ymin = ~unif.conf.low, ymax = ~unif.conf.high, fill = "95% Uniform"), alpha = 0.4) +
    geom_ribbon(aes_(ymin = ~theta.conf.low, ymax = ~theta.conf.high, fill = "95% Point-wise"), alpha = 0.35) +
    geom_line() +
    scale_fill_manual(values = flls) +
    labs(x = "Time",
                  y = "Difference",
                  fill = NULL)

  wrap_plots(top, bottom, ncol = 1, heights = c(3, 1))
}
