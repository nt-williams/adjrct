
#' Plot Restricted Mean Survival Time Estimates
#'
#' @param x
#' @param ...
#'
#' @return
#'
#' @importFrom ggplot2 ggplot aes geom_step labs scale_color_manual geom_ribbon scale_fill_manual theme element_blank
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @examples
plot.rmst <- function(x, ...) {
  estm <- all_estimates(x)
  clrs <- c("Treatment" = "blue", "Control" = "red")
  flls <- c("95% Uniform" = "lightblue", "95% Pointwise" = "blue")
  top <- ggplot2::ggplot(estm, aes(x = horizon)) +
    ggplot2::geom_step(aes(y = treatment, color = "Treatment")) +
    ggplot2::geom_step(aes(y = control, color = "Control")) +
    ggplot2::labs(x = NULL,
                  y = "Restricted Mean Survival Time",
                  color = NULL) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::theme(axis.ticks.x.bottom = element_blank(),
                   axis.text.x.bottom = element_blank())

  bottom <- ggplot2::ggplot(estm, aes(x = horizon, y = theta)) +
    ggplot2::geom_ribbon(aes(ymin = unif.conf.low, ymax = unif.conf.high, fill = "95% Uniform"), alpha = 0.4) +
    ggplot2::geom_ribbon(aes(ymin = theta.conf.low, ymax = theta.conf.high, fill = "95% Pointwise"), alpha = 0.35) +
    ggplot2::geom_line() +
    ggplot2::scale_fill_manual(values = flls) +
    ggplot2::labs(x = "Time horizon",
                  y = "Difference",
                  fill = NULL)

  patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(3, 1))
}

#' Plot Survival Probability Estimates
#'
#' @param x
#' @param ...
#'
#' @return
#'
#' @importFrom ggplot2 ggplot aes geom_step labs scale_color_manual geom_ribbon scale_fill_manual theme element_blank
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @examples
plot.survprob <- function(x, ...) {
  estm <- all_estimates(x)
  clrs <- c("Treatment" = "blue", "Control" = "red")
  flls <- c("95% Uniform" = "lightblue", "95% Point-wise" = "blue")
  top <- ggplot2::ggplot(estm, aes(x = horizon)) +
    ggplot2::geom_step(aes(y = treatment, color = "Treatment")) +
    ggplot2::geom_step(aes(y = control, color = "Control")) +
    ggplot2::labs(x = NULL,
                  y = "Survival probability",
                  color = NULL) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::theme(axis.ticks.x.bottom = element_blank(),
                   axis.text.x.bottom = element_blank())

  bottom <- ggplot2::ggplot(estm, aes(x = horizon, y = theta)) +
    ggplot2::geom_ribbon(aes(ymin = unif.conf.low, ymax = unif.conf.high, fill = "95% Uniform"), alpha = 0.4) +
    ggplot2::geom_ribbon(aes(ymin = theta.conf.low, ymax = theta.conf.high, fill = "95% Point-wise"), alpha = 0.35) +
    ggplot2::geom_line() +
    ggplot2::scale_fill_manual(values = flls) +
    ggplot2::labs(x = "Time",
                  y = "Difference",
                  fill = NULL)

  patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(3, 1))
}
