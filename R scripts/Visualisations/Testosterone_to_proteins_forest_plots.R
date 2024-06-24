# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create a data frame with the provided information
data <- data.frame(
  Protein = c("Leptin", "Sclerostin", "Fatty acid-binding protein, adipocyte", "Interleukin-1 receptor-like 1", 
              "Contactin-1", "Carbonic anhydrase 1", "Prokineticin-1", "C-type natriuretic peptide", 
              "T-cell differentiation antigen CD6", "Cytokine receptor-like factor 1", "Glycodelin", 
              "WNT1-inducible-signaling pathway protein 2", "Stromelysin-1", "Flavin reductase", 
              "Carbonic anhydrase 2", "Stanniocalcin-2", "Corticoliberin", 
              "ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 1", "T-cell surface glycoprotein CD5"),
  Estimate = c(-0.193745432, 0.198693639, -0.191443167, 0.171792927, -0.173797397, 
               0.146394379, 0.267706035, 0.160373517, -0.141951956, 0.201803014, 
               0.141792523, 0.219632791, 0.231993435, 0.14577985, 0.142940945, 
               0.205042566, -0.24111977, 0.152005171, -0.162033324),
  SE = c(0.039547085, 0.03904318, 0.04290187, 0.033148395, 0.034205222, 
         0.03391447, 0.036514221, 0.035402537, 0.029876269, 0.045309963, 
         0.025087521, 0.035461876, 0.032985449, 0.031346464, 0.031683879, 
         0.042767868, 0.045099107, 0.033590086, 0.034655065)
)

# Calculate OR and 95% CI
data <- data %>%
  mutate(OR = exp(Estimate),
         CI_lower = exp(Estimate - 1.96 * SE),
         CI_upper = exp(Estimate + 1.96 * SE))

# Create forest plot with enhanced aesthetics
ggplot(data, aes(x = reorder(Protein, OR), y = OR)) +
  geom_point(size = 3, shape = 21, fill = "blue", color = "black") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.3, color = "black") +
  coord_flip() +
  labs(title = "MR estimates for testosterone\nassociations with proteins in men",
       x = "Proteins",
       y = "Odds Ratio (95% Confidence Interval)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_log10(breaks = c(0.5, 1, 2, 3), labels = c("0.5", "1", "2", "3"))