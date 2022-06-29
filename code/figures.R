# MappingDirect
# Figures
# 2022 Katrin Petermann

## grey, green, blue, yellow, pink
cols <- c("#677078", "#009870", "#6ca5da", "#f8cf3b", "#ea8698")

#######################################################
### Baseline Values
#######################################################

p_sex <- ggplot(data = MDSTNdata, aes(x = gender))+
  geom_bar(alpha = 0.5, fill = cols[2], color = cols[2])+
  labs(x = "Gender", y = "Number of Patients")+
  theme_pubr(legend = "none")

p_age <- ggplot(data = MDSTNdata, aes(x = age.at.surgery))+
  geom_histogram(breaks = seq(30, 90, 5),
                 fill = cols[2], color = cols[2], 
                 alpha = 0.5)+
  scale_x_continuous(breaks = seq(30, 80, 10), limits = c(30, 80))+
  labs(x = "Age at Surgery", y = "Number of Patients")+
  theme_pubr(legend = "none")

p_duration <- ggplot(data = MDSTNdata, aes(x = duration))+
  geom_histogram(breaks = seq(0, 30, 2),
                 fill = cols[2], color = cols[2], 
                 alpha = 0.5)+
  scale_x_continuous(breaks = seq(0, 30, 10), limits = c(0, 25))+
  scale_y_continuous(breaks = seq(0, 14, 2), limits =c(0, 10))+
  labs(x = "Disease Duration [years]", y = "Number of Patients")+
  theme_pubr(legend = "none")

p_mapping <- ggplot(data = MDSTNdata, aes(x = time.after.surgery))+
  geom_histogram(breaks = seq(14, 36, 2),
                 fill = cols[2], color = cols[2],
                 alpha = 0.5)+
  scale_x_continuous(breaks = seq(14, 36, 2))+
  labs(x = "Time between Surgery and Mapping [weeks]", y = "Number of Patients")+
  theme_pubr(legend = "none")

p_mu3_pre_off <- ggplot(data = MDSTNdata, aes(x = MU3.pre.off))+
  geom_histogram(breaks = seq(0, 90, 10),
                 fill = cols[2], color = cols[2],
                 alpha = 0.5)+
  scale_x_continuous(breaks = seq(0, 90, 10))+
  labs(x = "MDS-UPDRS-III: Preoperative, Off-Medication", y = "Number of Patients")+
  theme_pubr(legend = "none")

p_mu3_pre_on <- ggplot(data = MDSTNdata, aes(x = MU3.pre.on))+
  geom_histogram(breaks = seq(0, 90, 10),
                 fill = cols[2], color = cols[2],
                 alpha = 0.5)+
  scale_x_continuous(breaks = seq(0, 90, 10))+
  labs(x = "MDS-UPDRS-III: Preoperative, On-Medication", y = "Number of Patients")+
  theme_pubr(legend = "none")

p_mu3_improvement_pre <- ggplot(data = MDSTNdata, aes(x = mu3_improvement_pre))+
  geom_histogram(breaks = seq(30, 100, 10),
                 fill = cols[2], color = cols[2],
                 alpha = 0.5)+
  scale_x_continuous(breaks = seq(30, 100, 10))+
  labs(x = "MDS-UPDRS-III: Preoperative, % Improvement", y = "Number of Patients")+
  theme_pubr(legend = "none")

p_ledd_pre <- ggplot(data = MDSTNdata, aes(x = ledd.pre))+
  geom_histogram(breaks = seq(0, 2200, 200),
                 fill = cols[2], color = cols[2],
                 alpha = 0.5)+
  scale_x_continuous(breaks = seq(0, 2200, 400))+
  labs(x = "LEDD: Preoperative", y = "Number of Patients")+
  theme_pubr(legend = "none")


pdf(file = "../figures/baseline_data.pdf", width = 10, height = 8)
print(ggarrange(p_sex, p_mu3_pre_off, p_age, 
                p_mu3_pre_on, p_duration, p_mu3_improvement_pre, 
                p_mapping, p_ledd_pre, 
                nrow = 4, ncol = 2))
dev.off()


# ************************************
# Figure 1: Difference in TW
# ************************************

lmeTWdiff <- lmer(diffTW ~ direction + (1 | hemisphere), data = longData)
(summaryLmeBLTW <- summary(lmeTWdiff))
(lmeConfint <- confint(lmeTWdiff))

p_values <- summaryLmeBLTW$coefficients[2:4, "Pr(>|t|)"]
p_stars <- symnum(p_values, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                  symbols = c("****", "***", "**", "*", "ns"))

df1 <- data.frame(a = c(1, 1, 2, 2), b = c(7.1, 7.2, 7.2, 7.1))
df2 <- data.frame(a = c(1, 1, 3, 3), b = c(7.7, 7.9, 7.9, 7.7))

plot_1 <- ggplot(data = longData[longData$direction != "ring level", ], 
                       aes(x = direction, y = diffTW, fill = direction))+
  geom_boxjitter(jitter.shape = 21, jitter.color = NA,
                 outlier.colour = NULL, outlier.shape = 19,
                 outlier.size = 1, width = 0.3,
                 jitter.params = list(width = 0.06, height = 0.05),
                 jitter.alpha = 0.4, jitter.size = 1,
                 errorbar.draw = T,
                 errorbar.length = 0.2)+
  geom_flat_violin(position = position_nudge(x = 0.17, y = 0), alpha = 0.8,
                   width = 0.4) +
  geom_hline(yintercept = 0, lty = 2, col = "lightgrey")+
  scale_fill_manual(values = c("#ffa600", "#bc5090", "#003f5c"))+
  scale_color_manual(values = c("#ffa600", "#bc5090", "#003f5c"))+
  annotate("text", x = 1:3, y = 4.5, label = p_stars) +
  annotate("rect", xmin = 0.82, xmax = 1.35, ymin = 0.25, ymax = 4, 
           alpha = .2)+
  theme_pubr()+
  theme(legend.position = "")+
  labs(x = "", y = expression(Delta~"Therapeutic Window [mA]"))+
  scale_y_continuous(breaks = seq(-4, 4.5, by = 1), limits = c(-4, 4.5))

pdf(file = "../figures/fig_1.pdf", width = 8, height = 4)
print(plot_1)
dev.off()

jpeg(file = "../figures/fig_1.jpg", width = 8, height = 4, units = "in", res = 300)
print(plot_1)
dev.off()


# ************************************
# Figure 2: Distribution of ET, SET, TW
# ************************************

measure.labs =  c("Therapeutic Window", "Side Effect Threshold", "Effect Threshold")
names(measure.labs) = c("TW", "SET", "ET")

plot_2 <- ggplot(data = dd_long, aes(x = TW_larger, y = amplitude, fill = TW_larger))+
  geom_boxjitter(jitter.shape = 21, jitter.color = NA,
                 outlier.colour = NULL, outlier.shape = 19,
                 outlier.size = 1,
                 jitter.params = list(width = 0.15, height = 0.05),
                 jitter.alpha = 0.4, jitter.size = 1,
                 errorbar.draw = T,
                 errorbar.length = 0.2)+
  facet_grid(cols = vars(measure), labeller = labeller(measure = measure.labs))+
  scale_fill_manual(values = coldef[c(1,3)])+
  scale_color_manual(values = coldef[c(1,3)])+
  labs(x = "Benefit from directional Testing", y = "Stimulation Amplitude (mA)")+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 12))

pdf(file = "../figures/fig_2.pdf", width = 8, height = 3)
print(plot_2)
dev.off()

jpeg(file = "../figures/fig_2.jpg", width = 8, height = 3, units = "in", res = 300)
print(plot_2)
dev.off()


# ************************************
# Figure 3: Prediction Model
# ************************************

## ROC plot
plot_3a <- ggroc(list(roc_ET, roc_SET, roc_TW, roc_TW_ET), size = 1)+ 
  theme_pubr(legend = c(0.75, 0.3))+
  theme(legend.background = element_rect(fill = NA))+
  grids()+
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + 
  labs(color = "")+
  scale_color_manual(values = coldef, 
                     labels = c("Therapeutic Window", 
                                "Side Effect Threshold", 
                                "Effect Threshold", 
                                "Effect Threshold + \nTherapeutic Window"))+
  coord_equal()

for(i in 1:4) {
  plot_3a <- plot_3a +
    geom_ribbon(
      data = dd_roc[[i]],
      aes(x = x, ymin = lower, ymax = upper),
      fill = coldef[i],
      alpha = 0.1,
      inherit.aes = F)
}

plot_3a <- plot_3a + labs(x = "Specificity", y = "Sensitivity")

## Predictive Performance
plot_3b <- ggplot(data = results_CV, aes(x = result, y = estimate, col = measure))+
  geom_point(position = position_dodge(width = 0.5), 
             size = 3,
             shape = 19)+
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.5),
                size = 1, width = 0.5)+
  lims(y = c(0, 1))+
  labs(x = "", y = "", col = "")+
  scale_color_manual(values = coldef, 
                     labels = c("Effect Threshold", 
                                "Side Effect Threshold", 
                                "Therapeutic Window", 
                                "Effect Threshold + \nTherapeutic Window"))+
  theme_pubr(legend = "right")+
  grids()

pdf(file = "../figures/fig_3.pdf", width = 11, height = 3.8)
print(ggarrange(plot_3a, plot_3b, nrow = 1, ncol = 2, 
                widths = c(1, 1.4), labels = c("A", "B")))
dev.off()

jpeg(file = "../figures/fig_3.jpg", width = 11, height = 3.8, units = "in", res = 300)
print(ggarrange(plot_3a, plot_3b, nrow = 1, ncol = 2, 
                widths = c(1, 1.4), labels = c("A", "B")))
dev.off()
