library(forestplot)
library(dplyr)









data <- data.frame(
  Hormone = c("SHBG", "SHBG", "Testosterone", "Testosterone"),
  Sex = c("Male", "Female", "Male", "Female"),
  SNPs = c(306, 322, 96, 220),
  OR = c(0.237, 0.334, -0.165, -0.01),
  se = c(0.061, 0.073, 0.040, 0.039),
  p = c(0.0001, 6.16e-6, 7.92e-5, 0.806)
  
)





# Create a data frame for the forest plot
forest_data <- data.frame(
  label = paste(data$Hormone, data$Sex),
  estimate = data$OR,
  lo = data$OR - 1.96 * data$se,
  hi = data$OR + 1.96 * data$se,
  pvalue = data$p
)



forestplot(
  mean = forest_data$estimate,
  lower = forest_data$lo,
  upper = forest_data$hi,
  labeltext = forest_data$label,
  title = "Forest Plot",
  xlab = "Odds Ratio",
  col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
  pvalues = forest_data$pvalue,
  txt_gp = fpTxtGp(label = gpar(fontsize = 10)),
  xticks = c(-0.5, 0, 0.5),
  xticklabels = c("0.5", "1", "1.5")
)






# Create a data frame for the forest plot
forest_data <- data.frame(
  label = paste(data$Hormone, data$Sex),
  estimate = data$OR,
  lo = data$OR - 1.96 * data$se,
  hi = data$OR + 1.96 * data$se,
  pvalue = data$p
)

# Define colors for males and females
male_color <- "lightblue"
female_color <- "pink"

# Create the forest plot with different colors for males and females
forestplot(
  mean = forest_data$estimate,
  lower = forest_data$lo,
  upper = forest_data$hi,
  labeltext = forest_data$label,
  title = "Associations between sex hormones and CAD",
  xlab = "Odds Ratio",
  col = fpColors(
    box = "darkblue",
    line = "black",
  ),
  pvalues = forest_data$pvalue,
  txt_gp = fpTxtGp(label = gpar(fontsize = 10)),
  xticks = c(-0.5, 0, 0.5),
  xticklabels = c("0.5", "1", "1.5"),
)

