#Install or load needed libraries

if(!require(readr)){
  install.packages("readr")
  library(readr)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(rstatix)){
  install.packages("rstatix")
  library(rstatix)
}
if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}
if(!require(lmerTest)){
  install.packages("lmerTest")
  library(lmerTest)
}
if(!require(emmeans)){
  install.packages("emmeans")
  library(emmeans)
}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(bpnreg)){
  install.packages("bpnreg")
  library(bpnreg)
}

#Load spatiotemporal data from Github
orienting_pcpa_spatiotemporal <- read_csv("https://raw.githubusercontent.com/lanec-unifesspa/5-HT-CAS/master/orienting/orienting_pcpa_spatiotemporal.csv",
col_types = cols(Drug = col_factor(levels = c("VEH",
"pCPA")), Treatment = col_factor(levels = c("CTRL",
"SA")), Type = col_factor(levels = c("0",
"1"))))
View(orienting_pcpa_spatiotemporal)

#Conduct repeated-measures ANOVA on 'Time spent near the stimulus'
#Summary statistics
orienting_pcpa_spatiotemporal %>% group_by(Drug, Treatment, Trial, Type) %>% get_summary_stats(TimeNearStimulus, type = "mean_sd")

#Fit a mixed model to the data, run ANOVA and post-hoc tests
fit_lm_Time_pcpa <- lmer(TimeNearStimulus ~ Drug*Treatment*Trial*Type + (1|Subject), data = orienting_pcpa_spatiotemporal)
anova(fit_lm_Time_pcpa)
emmeans(fit_lm_Time_pcpa, list(pairwise ~ Drug), adjust = "tukey")
emmeans(fit_lm_Time_pcpa, list(pairwise ~ Treatment), adjust = "tukey")
emmeans(fit_lm_Time_pcpa, "Trial", adjust = "tukey")
emmeans(fit_lm_Time_pcpa, list(pairwise ~ Type), adjust = "tukey")

#Plot
grouped_pcpa <- group_by(orienting_pcpa_spatiotemporal, Drug, Treatment, Trial, Type) #Prepare data for plotting by grouping by independent variables
Time_grouped_pcpa <- summarise(grouped_pcpa, mean = mean(TimeNearStimulus, na.rm = TRUE), sd = sd(TimeNearStimulus, na.rm = TRUE)) #Extract means and standard deviations for Time Spent Near Stimulus
ggplot(Time_grouped_pcpa, aes(Trial))  + geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Drug, alpha = 0.2)) + geom_line(aes(y = mean, color = Drug)) + facet_wrap(~Treatment) + scale_fill_discrete(name = "Drug") + scale_color_discrete(guide = FALSE) + scale_alpha_continuous(guide = "none") + labs(x = "Trial", y = "Time near stimulus (s)") + theme_bw() #Produce line and ribbon plots with means per trial as lines and SDs as ribbons

#Conduct repeated-measures ANOVA on 'Erratic Swimming'
#Summary statistics
orienting_pcpa_spatiotemporal %>% group_by(Drug, Treatment, Trial, Type) %>% get_summary_stats(AbsoluteTurnAngle, type = "mean_sd")

#Fit a mixed model to the data, run ANOVA and post-hoc tests
fit_lm_Erratic_pcpa <- lmer(AbsoluteTurnAngle ~ Drug*Treatment*Trial*Type + (1|Subject), data = orienting_pcpa_spatiotemporal)
anova(fit_lm_Erratic_pcpa)
emmeans(fit_lm_Erratic_pcpa, list(pairwise ~ Drug), adjust = "tukey")
emmeans(fit_lm_Erratic_pcpa, list(pairwise ~ Treatment), adjust = "tukey")
emmeans(fit_lm_Erratic_pcpa, "Trial", adjust = "tukey")
emmeans(fit_lm_Erratic_pcpa, list(pairwise ~ Type), adjust = "tukey")

#Plot
Erratic_grouped_pcpa <- summarise(grouped_pcpa, mean = mean(AbsoluteTurnAngle, na.rm = TRUE), sd = sd(AbsoluteTurnAngle, na.rm = TRUE)) #Extract means and standard deviations for Erratic swimming
ggplot(Erratic_grouped_pcpa, aes(Trial))  + geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Drug, alpha = 0.2)) + geom_line(aes(y = mean, color = Drug)) + facet_wrap(~Treatment) + scale_fill_discrete(name = "Drug") + scale_color_discrete(guide = FALSE) + scale_alpha_continuous(guide = "none") + labs(x = "Trial", y = "Erratic swimming (Abs. turn angle)") + theme_bw() #Produce line and ribbon plots with means per trial as lines and SDs as ribbons

#Conduct repeated-measures ANOVA on 'Swimming speed'
#Summary statistics
orienting_pcpa_spatiotemporal %>% group_by(Drug, Treatment, Trial, Type) %>% get_summary_stats(Speed, type = "mean_sd")

#Fit a mixed model to the data, run ANOVA and post-hoc tests
fit_lm_Speed_pcpa <- lmer(Speed ~ Drug*Treatment*Trial*Type + (1|Subject), data = orienting_pcpa_spatiotemporal)
anova(fit_lm_Speed_pcpa)
emmeans(fit_lm_Speed_pcpa, list(pairwise ~ Drug), adjust = "tukey")
emmeans(fit_lm_Speed_pcpa, list(pairwise ~ Treatment), adjust = "tukey")
emmeans(fit_lm_Speed_pcpa, "Trial", adjust = "tukey")
emmeans(fit_lm_Speed_pcpa, list(pairwise ~ Type), adjust = "tukey")

#Plot
Speed_grouped_pcpa <- summarise(grouped_pcpa, mean = mean(Speed, na.rm = TRUE), sd = sd(Speed, na.rm = TRUE)) #Extract means and standard deviations for Speed
ggplot(Speed_grouped_pcpa, aes(Trial))  + geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = Drug, alpha = 0.2)) + geom_line(aes(y = mean, color = Drug)) + facet_wrap(~Treatment) + scale_fill_discrete(name = "Drug") + scale_color_discrete(guide = FALSE) + scale_alpha_continuous(guide = "none") + labs(x = "Trial", y = "Swimming speed (cm/s)") + theme_bw() #Produce line and ribbon plots with means per trial as lines and SDs as ribbons


#Load body orientation data from Github
orienting_pcpa_rproj <- read_csv("https://raw.githubusercontent.com/lanec-unifesspa/5-HT-CAS/master/orienting/orienting_pcpa_rproj.csv",
                                          col_types = cols(Drug = col_factor(levels = c("VEH",
                                                                                        "pCPA")), Treatment = col_factor(levels = c("CTRL",
                                                                                                                                    "SA")), Type = col_factor(levels = c("OFF",
                                                                                                                                                                         "ON"))))
View(orienting_pcpa_rproj)

#Conduct ANOVAs on 'Directional focus (Rproj)'
#Summary statistics
grouped_pcpa_rproj <- group_by(orienting_pcpa_rproj, Drug, Treatment, Trial, Type) #Prepare data for summary statistics
Rproj_grouped_pcpa <- summarise(grouped_pcpa_rproj, mean = mean(Rproj, na.rm = TRUE), sd = sd(Rproj, na.rm = TRUE))

#Fit a linear model with fixed effects
fit_lm_Rproj_pcpa <- aov(Rproj ~ Drug*Treatment*Type, data = orienting_pcpa_rproj)
anova(fit_lm_Rproj_pcpa)
emmeans(fit_lm_Rproj_pcpa, list(pairwise ~ Drug), adjust = "tukey")
emmeans(fit_lm_Rproj_pcpa, list(pairwise ~ Treatment), adjust = "tukey")
emmeans(fit_lm_Rproj_pcpa, list(pairwise ~ Type), adjust = "tukey")

#Conduct circular ANOVA
#Fit a circular model on the data
fit_cir_angle_pcpa <- bpnr(pred.I = Angle ~ 1 + Drug*Treatment*Type, data = orienting_pcpa_rproj, its = 1000, burn = 100, n.lag = 3, seed = 101)
fit_cir_angle_pcpa

#Plot data
ggplot(orienting_pcpa_rproj, aes(x = Angle, y = Rproj)) + coord_polar() + geom_point(aes(shape = Drug, color = Type)) + scale_color_manual(values = c(ON = "red", OFF = "blue"), name = "Trial type") + scale_shape(name = "Drug") + facet_wrap(.~Treatment) + theme_bw() + labs(y = "Directional focus (Rproj)", x ="")