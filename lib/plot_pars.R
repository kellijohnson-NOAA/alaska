#' Plot parameter values using histograms.
#'
#' @param data
#' @param pars The parameter values you want to plot.
#' @param subpopulations The subpopulations you want to plot.
#'
plot_pars <- function(data,
  pars = c("alpha", "sigma[epsilon]", "sigma[omega]", "scale", "rho"),
  subpopulations = 1, type = c("historgram", "boxplot")) {

  type <- match.arg(type, choices = c("historgram", "boxplot"),
    several.ok = FALSE)

  # Subset the data and reshape it.
  scen1 <- reshape_results(data[data$Year == 1 &
    data$Site == 1 &
    data$subpopulations %in% subpopulations, ])
  scen1 <- droplevels(scen1[scen1$par %in% pars, ])

  scen1$variable[which(scen1$variable == "Simulated_example")] <- "Poisson"
  scen1$variable[which(scen1$variable == "zeroinflatedlnorm")] <- "Poissonln"

  # Find if there are any really bad values
  bad <- which(scen1$em > scen1$om * 300)
  if (length(bad) > 0) {
    scen1 <- scen1[-bad, ]
    bad <- scen1[bad, ]
    bad <- aggregate(re ~ variable + model, data = bad, FUN = length)
    bad$par <- pars[1]
    bad$text <- paste0("(removed ", NROW(bad$re), ")")
    bad$em <- min(scen1$em)
    bad$y <- 0.05
  }

  if (type == "boxplot"){
    message("Need to finish coding this.")
    g <- ggplot(data = scen1, aes(y = em)) +
        facet_grid(par ~ variable,
          scales = "free_y",
          labeller = label_parsed
          ) +
        xlab("estimate") + ylab("density") +
        geom_boxplot(aes(x = model))
      means <- calc_means(data = scen1, digits = 2, g = g, yscalar = 1.1)
    g <- g + geom_hline(data = means, aes(yintercept = mean), linetype = "dashed",
        col = "red")
    return(g)
  }

  g <- ggplot(data = scen1, aes(x = em)) +
    facet_grid(model + variable ~ par, labeller = label_parsed,
      scales = "free_x") +
    xlab("estimate") + ylab("density") +
    geom_histogram(aes(y = (..count..)/sum(..count..)),
      colour = "black", fill = "white", bins = 30) +
    scale_x_continuous(labels = fmt) +
    theme(plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "white", fill = NA, size = 1),
      legend.key = element_rect(colour = "white"),
      legend.title = element_text(size = 0, face = "bold"),
      legend.text = element_text(size = 7, face = "bold"),
      legend.position = "none",
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 6))

  means <- calc_means(data = scen1, digits = 2, g = g, yscalar = 1.1)
  g <- g +
    geom_text(data = means, mapping = aes(x = x, y = y, label = MARE),
      hjust = 1.0, vjust = 1.90, size = 3) +
    geom_text(data = means, mapping = aes(x = x, y = y, label = MRE),
      hjust = 1.0, vjust = 0.60, size = 3) +
    geom_vline(data = means, aes(xintercept = mean), linetype = "dashed",
      col = "red")
  if (NROW(bad) > 0) {
    g <- g +
      geom_text(data = bad, mapping = aes(x = em, y = y, label = text),
        hjust = 0.05, size = 3)
  }
  return(g)
}

