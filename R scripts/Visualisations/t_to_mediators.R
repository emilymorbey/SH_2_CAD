library(ggplot2)


library(ggplot2)
library(ggpubr)

# Data frame for grouped Blood Pressure measurements
data_bp <- data.frame(
  Outcome = rep(c("Diastolic Blood Pressure", "Systolic Blood Pressure"), each = 2),
  Method = rep(c("IVW", "MR-Egger"), 2),
  Estimate = c(0.049642372, 0.034789588,
               0.03340528, 0.034108999),
  StdError = c(0.015143273, 0.027970136,
               0.015166549, 0.028082301),
  LCI = c(0.019962102, -0.020030872,
          0.003679391, -0.020931299),
  UCI = c(0.079322641, 0.089610048,
          0.063131169, 0.089149298),
  Pvalue = c(0.001044788, 0.213568801,
             0.02762556, 0.224515475)
)

data_ivw <- subset(data_bp, Method == "IVW")

# Create the plot
ivw_plot <- ggplot(data_ivw, aes(x = Outcome, y = Estimate, ymin = LCI, ymax = UCI, color = Method)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "MR from Testosterone to Blood\n Pressure (Diastolic and Systolic)",
       x = "Outcome",
       y = "Estimate") +
  theme_classic() +
  theme(
    text = element_text(size = 15, family = "Arial"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  ) +
  scale_color_manual(values = c("darkblue")) +
  coord_flip()

# Print the plot
print(ivw_plot)

# Data frame for grouped lipid measures
data_chol_trig <- data.frame(
  Outcome = rep(c("HDL Cholesterol", "LDL Cholesterol", "Triglycerides"), each = 2),
  Method = rep(c("IVW", "MR-Egger"), 3),
  Estimate = c(-0.038087951, -0.006039741,
               -0.00985302, 0.00333768,
               0.019912407, 0.045863243),
  StdError = c(0.016176409, 0.035780225,
               0.010824837, 0.023978193,
               0.013846076, 0.030619251),
  LCI = c(-0.069793131, -0.076167694,
          -0.031069311, -0.043658715,
          -0.007225403, -0.014149387),
  UCI = c(-0.006382771, 0.064088212,
          0.011363271, 0.050334074,
          0.047050216, 0.105875872),
  Pvalue = c(0.018545807, 0.86595311,
             0.362704719, 0.889294906,
             0.1503982, 0.134170554)
)

data_ivw_2 <- subset(data_chol_trig, Method == "IVW")




# Create the plot
ivw_plot_2 <- ggplot(data_ivw_2, aes(x = Outcome, y = Estimate, ymin = LCI, ymax = UCI, color = Method)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "MR from Testosterone to Lipids",
       x = "Outcome",
       y = "Estimate") +
  theme_classic() +
  theme(
    text = element_text(size = 15, family = "Arial"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  ) +
  scale_color_manual(values = c("darkblue")) +
  coord_flip()

print(ivw_plot_2)


plot <- ggarrange(ivw_plot, ivw_plot_2, nrow = 1, ncol = 2)


ggsave("t_to_mediators.png", plot = plot , width = 14, height = 6, dpi = 400)
