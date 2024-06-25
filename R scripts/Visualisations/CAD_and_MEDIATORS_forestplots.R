library(readxl)
library(ggplot2)
library(ggplot2)


############### CAD 

cad <- read_excel("Forest plots/cad.xlsx")

fp <- ggplot(data = cad, aes(x = Sex, y = OR, ymin = LOWERCI, ymax = UPPERCI, colour = Hormone)) +
  geom_pointrange(size = 1.2, position = position_dodge(width = 0.5)) + # Adjust width for dodge
  geom_errorbar(aes(ymin = LOWERCI, ymax = UPPERCI), width = 0.05, position = position_dodge(width = 0.5)) + # Add error bars
  geom_hline(yintercept = 1, lty = 2, lwd = 1) +
  coord_flip() +
  labs(y = "Odds ratio (95% Confidence Interval)", x = "Sex", title = "Associations of sex hormone levels with CAD") +
  theme_minimal() + # Minimal theme
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13, angle = 0, hjust = 0.5), # Adjust x-axis text angle and position
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + # Adjust plot title position and style
  scale_colour_manual(values = c("darkred", "darkturquoise")) +
  guides(colour = guide_legend(override.aes = list(size = 1))) # Adjust legend size

print(fp)






cad <- head(cad, 2)

fp <- ggplot(data = cad, aes(x = Sex, y = OR, ymin = LOWERCI, ymax = UPPERCI, colour = Hormone)) +
  geom_pointrange(size = 2, position = position_dodge(width = 0.5)) + # Adjust width for dodge
  geom_errorbar(aes(ymin = LOWERCI, ymax = UPPERCI), width = 0.05, position = position_dodge(width = 0.5)) + # Add error bars
  geom_hline(yintercept = 1, lty = 2, lwd = 1) +
  coord_flip() +
  labs(y = "Odds ratio (95% Confidence Interval)", x = "Sex", title = "Associations of sex hormone levels with CAD\nin males") +
  theme_minimal() + # Minimal theme
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 17, angle = 0, hjust = 0.5), # Adjust x-axis text angle and position
        axis.text.y = element_text(size = 17),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 19, hjust = 0.5, face = "bold")) + # Adjust plot title position and style
  scale_colour_manual(values = c("darkred", "darkturquoise")) +
  guides(colour = guide_legend(override.aes = list(size = 1))) # Adjust legend size

print(fp)





cad <- tail(cad, 2)

fp <- ggplot(data = cad, aes(x = Sex, y = OR, ymin = LOWERCI, ymax = UPPERCI, colour = Hormone)) +
  geom_pointrange(size = 2, position = position_dodge(width = 0.5)) + # Adjust width for dodge
  geom_errorbar(aes(ymin = LOWERCI, ymax = UPPERCI), width = 0.05, position = position_dodge(width = 0.5)) + # Add error bars
  geom_hline(yintercept = 1, lty = 2, lwd = 1) +
  coord_flip() +
  labs(y = "Odds ratio (95% Confidence Interval)", x = "Sex", title = "Associations of sex hormone levels with CAD\nin females") +
  theme_minimal() + # Minimal theme
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 17, angle = 0, hjust = 0.5), # Adjust x-axis text angle and position
        axis.text.y = element_text(size = 17),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size = 19, hjust = 0.5, face = "bold")) + # Adjust plot title position and style
  scale_colour_manual(values = c("darkred", "darkturquoise")) +
  guides(colour = guide_legend(override.aes = list(size = 1))) # Adjust legend size

print(fp)


############## MALE SHBG MEDIATORS ###################


  cad <- read_excel("Forest plots/maleshbg.xlsx")
  cad
  
  fp <- ggplot(data=cad, aes(x=Mediator, y=beta, ymin=LOWERCI, ymax=UPPERCI, colour = Mediator))+
    geom_pointrange(size=0.7, show.legend = FALSE)+
    geom_hline(yintercept = 0, lty=2, lwd=1)+
    coord_flip()+
    labs(x="Mediator", y="Unit change in mediator\nwith 1 unit increase in SHBG", title = "Male SHBG")+
    theme_bw()+
    theme(axis.title = element_text(size = 16),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          title = element_text(size=18))+
    theme(axis.title.y = element_text(vjust = -3))
  
  print(fp)
fp + scale_colour_manual(values = c("darkorange", "darkturquoise", "royalblue", "forestgreen", "red"))
  
  
  
  
  
  ############## FEMALE SHBG MEDIATORS ###################
  
  
  cad <- read_excel("Forest plots/femaleshbg.xlsx")
  cad
  
  fp <- ggplot(data=cad, aes(x=Mediator, y=beta, ymin=LOWERCI, ymax=UPPERCI, colour = Mediator))+
    geom_pointrange(size=0.7, show.legend = FALSE, position = position_dodge(width = 0.5))+
    geom_hline(yintercept = 0, lty=2, lwd=1)+
    coord_flip()+
    labs(x="Mediator", y="Unit change in mediator\nwith 1 unit increase in SHBG", title = "Female SHBG")+
    theme_bw()+
    theme(axis.title = element_text(size = 16),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          title = element_text(size=18))+
    theme(axis.title.y = element_text(vjust = -3))
  
  print(fp)
  
fp + scale_colour_manual(values = c("darkorange", "darkturquoise", "royalblue", "forestgreen", "red"))
  
  
  
  
  
  ###### Male testosterone mediators #####################
  
  
  
cad <- read_excel("Forest plots/malet.xlsx")
cad
  
cad$Mediator <- as.factor(cad$Mediator)
  
fp <- ggplot(data=cad, aes(x=Mediator, y=beta, ymin=LOWERCI, ymax=UPPERCI, colour = Mediator))+
    geom_pointrange(size=2, show.legend = FALSE)+
    geom_errorbar(aes(ymin = LOWERCI, ymax = UPPERCI), width = 0.1, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, lty=2, lwd=1)+
    coord_flip()+
    labs(x="Mediator", y="Unit change in mediator with\n1 unit increase in testosterone", title = "Male Testosterone")+
    theme_minimal()+
    theme(axis.title = element_text(size = 20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          title = element_text(size=22))
  
print(fp)
fp + scale_colour_manual(values = c("red", "darkturquoise", "forestgreen")) +
  coord_cartesian(clip = "on")
  


cad <- read_excel("Forest plots/maletbp.xlsx")
cad

cad$Mediator <- as.factor(cad$Mediator)

fp <- ggplot(data=cad, aes(x=Mediator, y=beta, ymin=LOWERCI, ymax=UPPERCI, colour = Mediator))+
  geom_pointrange(size=2, show.legend = FALSE)+
  geom_errorbar(aes(ymin = LOWERCI, ymax = UPPERCI), width = 0.1, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, lty=2, lwd=1)+
  coord_flip()+
  labs(x="Mediator", y="Unit change in mediator with\n1 unit increase in testosterone")+
  theme_minimal()+
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        title = element_text(size=22))

print(fp)
fp + scale_colour_manual(values = c("orange", "royalblue")) +
  coord_cartesian(clip = "on")
  
  
  
  
  
  
  
  ### female testosterone mediators ########
  
  
  
  
  cad <- read_excel("Forest plots/femalet.xlsx")
  cad
  
  fp <- ggplot(data=cad, aes(x=Mediator, y=beta, ymin=LOWERCI, ymax=UPPERCI, colour = Mediator))+
    geom_pointrange(size=0.7, show.legend = FALSE)+
    geom_hline(yintercept = 0, lty=2, lwd=1)+
    coord_flip()+
    labs(x="Mediator", y="Unit change in mediator\nwith 1 unit increase in testosterone", title = "Female Testosterone")+
    theme_bw()+
    theme(axis.title = element_text(size = 16),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          title = element_text(size=18))+
    theme(axis.title.y = element_text(vjust = -3))
  
  print(fp)
fp + scale_colour_manual(values = c("darkorange", "darkturquoise", "royalblue", "forestgreen", "red"))
  
  
