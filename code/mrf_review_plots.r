# Simulation Study results
library(tidyverse)
library(patchwork)
source("mcmc_functions.r")


# Response Functions Plot
# Monte Carlo Data
filename_i <- "mrf_sim_out/ising_n64_a0_0.csv"
ss_func_i_df <- read.csv(filename_i)
ss_func_i_df$alpha <- 0.0

p1 <- ggplot() +
  geom_segment(aes(y = 0.5, x = 0.0, xend = 1.2), linewidth = .5) +
  geom_vline(xintercept = 0.88, linetype = 2) +
  labs(x = expression(psi), title = "Proportion Black Units", y = "Mean") +
  theme_minimal()
p2 <- ggplot() +
  geom_line(data = ss_func_i_df, aes(x = beta, y = sqrt(vT1) / N)) +
  geom_vline(xintercept = 0.88, linetype = 2) +
  labs(x = expression(psi), y = "Standard Deviation") +
  theme_minimal()
p3 <- ggplot() +
  geom_line(data = ss_func_i_df, aes(x = beta, y = eT2i / n_edges)) +
  geom_vline(xintercept = 0.88, linetype = 2) +
  labs(x = expression(psi), title = "Proportion Matches", y = NULL) +
  theme_minimal()
p4 <- ggplot(data = ss_func_i_df) +
  geom_smooth(aes(x = beta, y = sqrt(vT2i) / n_edges),
    color = "black", linewidth = 0.5, se = FALSE, span = 0.2
  ) +
  geom_point(aes(x = beta, y = sqrt(vT2i) / n_edges), size = 0.25) +
  geom_vline(xintercept = 0.88, linetype = 2) +
  labs(x = expression(psi), y = NULL) +
  theme_minimal()

p_prf_alpha_0 <- (p1 | p3) / (p2 | p4) +
  plot_annotation(
    title = "Response Functions",
    # expression(Ising ~ Model ~ Response ~ Functions * "," * ~ alpha == 0),
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  ) +
  plot_layout(axis_titles = "collect", axes = "collect")

# Plot dimensions for saving
h <- 3
w <- 4
h_rf <- h * 1.5
w_rf <- w * 1.5
h_rf_long <- h_rf * 1.5

ggsave("figures_mrf_review/p_response_functions_ising_alpha_0.png",
  p_prf_alpha_0,
  height = h_rf, width = w_rf
)

# Plot of covariates used in simulation study
n_r <- 64
x <- rep(1:n_r, times = n_r)
y <- rep(1:n_r, each = n_r)
xy <- (x + y - n_r - 1) / (n_r - 1)
p_covariates <- ggplot() +
  geom_tile(aes(x = x, y = y, fill = xy)) +
  scale_fill_gradient2(low = "#882255", mid = "white", high = "#117733", midpoint = 0) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    axis.title = element_blank(), axis.ticks = element_blank(),
    axis.text = element_blank(), panel.grid = element_blank(),
    panel.background = element_blank(), panel.border = element_rect(color = "gray", fill = NA)
  ) +
  labs(fill = "x")
ggsave("figures_mrf_review/p_covariates.png", p_covariates, height = 4, width = 5)


# Prior predictive response functions for the MRF with covariates
# Function to generate dataframe for plotting
generate_pprf_df <- function(tag) {
  autol <- read.csv(paste0("mrf_sim_out/grad_autol_n64", tag))
  ising <- read.csv(paste0("mrf_sim_out/grad_ising_n64", tag))
  cente <- read.csv(paste0("mrf_sim_out/grad_cente_n64", tag))
  autol$Formulation <- "Autologistic"
  ising$Formulation <- "Ising"
  cente$Formulation <- "Centered Autol."
  rbind(autol, ising, cente)
}

# Function to generate prior predictive response function plot
generate_pprf_plot <- function(df_all, title = NULL, subtitle = NULL) {
  color_scheme <- scale_color_manual(values = c("Autologistic" = "#88CCEE", "Ising" = "#117733", "Centered Autol." = "#882255"))
  eT1 <- ggplot(data = df_all) +
    geom_smooth(aes(x = psi, y = eT1 / N, color = Formulation),
      span = 0.5, se = FALSE, size = 0.75
    ) +
    # geom_point(aes(x = psi, y = eT1 / N, color = Formulation), size = 0.25) +
    color_scheme +
    labs(x = expression(psi), title = "Proportion Black Units", y = "Mean") +
    guides(color = "none")
  vT1 <- ggplot(data = df_all) +
    geom_smooth(aes(x = psi, y = sqrt(vT1) / N, color = Formulation),
      span = 0.3, se = FALSE, size = 0.75
    ) +
    geom_point(aes(x = psi, y = sqrt(vT1) / N, color = Formulation), size = 0.25) +
    color_scheme +
    labs(x = expression(psi), y = "Standard Deviation")
  eT2i <- ggplot(data = df_all) +
    geom_line(aes(x = psi, y = eT2i / n_edges, color = Formulation), size = 0.75) +
    # geom_point(aes(x = psi, y = eT2i / n_edges, color = Formulation), size = 0.25) +
    color_scheme +
    labs(x = expression(psi), title = "Proportion Matches", y = NULL) +
    guides(color = "none")
  vT2i <- ggplot(data = df_all) +
    geom_smooth(aes(x = psi, y = sqrt(vT2i) / n_edges, color = Formulation),
      span = 0.25, se = FALSE, size = 0.75
    ) +
    geom_point(aes(x = psi, y = sqrt(vT2i) / n_edges, color = Formulation), size = 0.25) +
    # geom_point(aes(x=psi, y=sqrt(vT2i)/n_edges), size=0.75) +
    color_scheme +
    labs(x = expression(psi), y = NULL)
  MR <- ggplot(data = df_all) +
    geom_smooth(aes(x = psi, y = eMR, color = Formulation),
      span = 0.3, se = FALSE, size = 0.75
    ) +
    # geom_point(aes(x = psi, y = eMR, color = Formulation), size = 0.25) +
    # geom_point(aes(x=psi, y=eMR), size=0.75) +
    color_scheme +
    labs(x = expression(psi), title = "Misclassification Rate", y = NULL) +
    guides(color = "none")
  DC <- ggplot(data = df_all) +
    geom_smooth(aes(x = psi, y = eDC / N, color = Formulation),
      span = 0.3, se = FALSE, size = 0.75
    ) +
    # geom_point(aes(x = psi, y = eDC / N, color = Formulation), size = 0.25) +
    color_scheme +
    labs(x = expression(psi), title = "Dominant Color", y = NULL) +
    guides(color = "none")

  (eT1 | eT2i) / (vT1 | vT2i) / (DC | MR) + plot_layout(guides = "collect") +
    plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    ) &
    theme_minimal() &
    theme(legend.position = "bottom")
}

# Intercept, norm
title <- "Prior Predictive Response Functions"
df_int_norm <- generate_pprf_df("_adj0_0_b0_0_sd1_0int_norm.csv")
p_int_norm <- generate_pprf_plot(
  df_int_norm, title, NULL
  #"Direct Data Models with Covariate External Field"
  #expression(Centered ~ Covariates * "," ~ bold(beta) ~ "~" ~ N(0, 1))
)
# p_int_norm
ggsave("figures_mrf_review/p_prior_predictive_response_functions_int_norm.png", p_int_norm, height = h_rf_long, width = w_rf)
