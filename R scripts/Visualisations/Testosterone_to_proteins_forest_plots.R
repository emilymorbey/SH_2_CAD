library(ggplot2)
library(dplyr)



# Create a data frame with the provided data
data <- data.frame(
  Protein = c("Pancreatic secretory granule membrane major glycoprotein GP2", "Leptin", 
              "Sclerostin", "Fatty acid-binding protein, adipocyte", "Interleukin-1 receptor-like 1", 
              "Ectonucleotide pyrophosphatase/phosphodiesterase family member 2", "Contactin-1", 
              "Ecto-ADP-ribosyltransferase 3", "Prokineticin-1", "C-type natriuretic peptide", 
              "T-cell differentiation antigen CD6", "Cytokine receptor-like factor 1", "Glycodelin", 
              "WNT1-inducible-signaling pathway protein 2", "Stromelysin-1", "Flavin reductase", 
              "Stanniocalcin-2", "Corticoliberin", "ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 1", 
              "T-cell surface glycoprotein CD5"),
  Estimate = c(0.113730591, -0.16027249, 0.158691837, -0.148691414, 0.138023751, -0.101105742, 
               -0.129101931, 0.144332248, 0.210832951, 0.122876322, -0.110021639, 0.175474231, 
               0.102591413, 0.180581409, 0.188497155, 0.108626313, 0.175720436, -0.190962729, 
               0.123478842, -0.127832838),
  SE = c(0.025248618, 0.030323828, 0.030454623, 0.033313111, 0.025391417, 0.023232773, 
         0.027896456, 0.032257285, 0.028341783, 0.0275248, 0.024042844, 0.03454515, 
         0.019854443, 0.027625419, 0.02617611, 0.024597298, 0.032547019, 0.034908003, 
         0.026181014, 0.026807648),
  OR = c(1.12, 0.85, 1.17, 0.86, 1.15, 0.9, 0.9, 1.16, 1.23, 1.13, 0.9, 1.19, 
         1.11, 1.2, 1.21, 1.11, 1.19, 0.83, 1.13, 0.88),
  LCI = c(1.07, 0.8, 1.1, 0.81, 1.09, 0.86, 0.83, 1.08, 1.17, 1.07, 0.85, 1.11, 
          1.07, 1.13, 1.15, 1.06, 1.12, 0.77, 1.07, 0.83),
  UCI = c(1.18, 0.9, 1.24, 0.92, 1.21, 0.95, 0.93, 1.23, 1.31, 1.19, 0.94, 1.28, 
          1.15, 1.26, 1.27, 1.17, 1.27, 0.88, 1.19, 0.93),
  P_Value = c(1.54e-05, 5.64e-07, 7.84e-07, 1.83e-05, 2.91e-07, 2.85e-05, 
              9.42e-06, 1.74e-05, 1.62e-11, 1.82e-05, 1.16e-05, 1.39e-06, 
              9.49e-07, 1.60e-09, 5.68e-11, 2.21e-05, 3.44e-07, 2.46e-07, 
              6.54e-06, 5.23e-06)
)

# Create the forest plot
ggplot(data, aes(x = reorder(Protein, OR), y = OR)) +
  geom_point(shape = 15, size = 3) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
  scale_y_continuous(trans = 'log10', 
                     breaks = c(0.75, 0.85, 1, 1.15, 1.25, 1.35),
                     labels = c(0.75, 0.85, 1, 1.15, 1.25, 1.35)) +
  coord_flip() +
  labs(
    title = "MR of Testosterone Effects on Proteins",
    subtitle = "Odds Ratios (OR) with 95% Confidence Intervals",
    x = "Protein",
    y = "Odds Ratio (log scale)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),
    plot.caption = element_text(size = 10)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")





### proteins to CAD 

library(ggplot2)

# Example data (from your original request)
data <- data.frame(
  protein = c("Ecto-ADP-ribosyltransferase 3", 
              "WNT1-inducible-signaling pathway protein 2", 
              "ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 1", 
              "T-cell differentiation antigen CD6", 
              "Contactin-1", 
              "Corticoliberin", 
              "Fatty acid-binding protein, adipocyte", 
              "Interleukin-1 receptor-like 1", 
              "Stromelysin-1", 
              "Glycodelin", 
              "Prokineticin-1", 
              "Sclerostin"),
  beta = c(-0.008596, -0.017662, 0.010184, 0.000726, -0.004038, -0.0022, 
           -0.008891, -0.043223, 0.013935, 0.001677, 0.00223, -0.00038),
  se = c(0.009274, 0.008419, 0.0096, 0.008351, 0.01187, 0.022382, 
         0.013264, 0.027441, 0.008013, 0.010773, 0.010805, 0.008294)
)

# Calculate odds ratios and confidence intervals
data$odds_ratio <- exp(data$beta)
data$lower_ci <- exp(data$beta - 1.96 * data$se)
data$upper_ci <- exp(data$beta + 1.96 * data$se)

# Create a new order variable for plotting
data <- data[order(data$odds_ratio), ]

# Create the forest plot in ggplot2
ggplot(data, aes(x = reorder(protein, odds_ratio), y = odds_ratio)) +
  geom_point(shape = 15, size = 3, color = "#377EB8") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, color = "#377EB8") +
  scale_y_continuous(trans = 'log10',
                     breaks = c(0.75, 0.85, 1, 1.15, 1.25, 1.35),
                     labels = c(0.75, 0.85, 1, 1.15, 1.25, 1.35)) +
  coord_flip() +
  labs(
    title = "Effect of pQTLs on CAD",
    subtitle = "95% Confidence Intervals",
    x = "Protein",
    y = "Odds Ratio (log scale)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),
    plot.caption = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_blank()
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.5)