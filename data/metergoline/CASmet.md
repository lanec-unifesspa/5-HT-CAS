---
title: "R Notebook for effects of metergoline on zebrafish alarm reactions (LaNeC)"
author:
- Caio Maximino^[Universidade Federal do Sul e Sudeste do Pará]
- Monica Gomes Lima^[Universidade do Estado do Pará]
- Sueslene Prado Rocha^[Universidade do Estado do Pará]
output:
  github_document 
subtitle: From project "Role of serotonin on behavioral responses to alarm substance in zebrafish: A putative model for panic disorder" (DOI: 10.17605/OSF.IO/BK85D)
tags:
- fear
- zebrafish
- metergoline
abstract: |
  Behavioral data during and after exposure to conspecific alarm substance in zebrafish treated with metergoline (1 mg/kg)
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook for the  library

When you execute code within the notebook, the results appear beneath the code. 

* Load needed libraries:
```{r}
if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
}
if(!require(coin)){
    install.packages("coin")
    library(coin)
}
if(!require(RCurl)){
    install.packages("RCurl")
    library(RCurl)
}

if(!require(plyr)){
    install.packages("plyr")
    library(plyr)
}

if(!require(rcompanion)){
    install.packages("rcompanion")
    library(rcompanion)
}

if(!require(WRS2)){
    install.packages("WRS2")
    library(WRS2)
}

if(!require(psych)){
    install.packages("psych")
    library(psych)
}

if(!require(multcompView)){
    install.packages("multcompView")
    library(multcompView)
}

if(!require(ggpubr)){
    install.packages("ggpubr")
    library(ggpubr)
}
```

* Load data
```{r}
x1 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/5-HT-CAS/master/data/metergoline/MET_during_exposure.csv")
MET_during_exposure <- read.csv(text = x1)
# Ensure factor variables are interpreted as such
MET_during_exposure$Treatment <- factor(MET_during_exposure$Treatment, levels = c("CTRL", "CAS"))
MET_during_exposure$Drug <- as.factor(MET_during_exposure$Drug)
View(MET_during_exposure)

x2 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/5-HT-CAS/master/data/metergoline/MET_after_exposure.csv")
MET_after_exposure <- read.csv(text = x2)
# Ensure factor variables are interpreted as such
MET_after_exposure$Treatment <- factor(MET_after_exposure$Treatment, levels = c("CTRL", "CAS"))
MET_after_exposure$Drug <- as.factor(MET_after_exposure$Drug)
View(MET_after_exposure)
```

* Detect and remove outliers, based on median absolute deviations
```{r}
gA <- MET_during_exposure[c(1:15), ]
median(gA$Time_On_Bottom) + 3*mad(gA$Time_On_Bottom)
median(gA$Time_On_Bottom) - 3*mad(gA$Time_On_Bottom)
gB <- MET_during_exposure[c(16:30), ]
median(gB$Time_On_Bottom) + 3*mad(gB$Time_On_Bottom)
median(gB$Time_On_Bottom) - 3*mad(gB$Time_On_Bottom)
gC <- MET_during_exposure[c(31:45), ]
median(gC$Time_On_Bottom) + 3*mad(gC$Time_On_Bottom)
median(gC$Time_On_Bottom) - 3*mad(gC$Time_On_Bottom)
gD <- MET_during_exposure[c(46:60), ]
median(gD$Time_On_Bottom) + 3*mad(gD$Time_On_Bottom)
median(gD$Time_On_Bottom) - 3*mad(gD$Time_On_Bottom)

#Subjects 4 from group CTRL+VEH, and 5 from group CTRL+MET detected as outlier; remove
MET_during_exposure <- MET_during_exposure[-c(4, 20), ]
MET_after_exposure <- MET_after_exposure[-c(4, 20), ]
View(MET_during_exposure)
View(MET_after_exposure)
```

* Run bootstrapped ANOVA for main and interaction effects on 2-way ANOVA, in data for observation DURING CAS exposure (based on https://rcompanion.org/rcompanion/d_08a.html)

-Create seed for computational reproducibility
```{r}
set.seed(42)
```

-Create interaction variable for pairwise permutation tests (post-hoc tests); the created variable will be used across all tests
```{r}
MET_during_exposure$MET_during_exposure.int <- interaction(MET_during_exposure$Treatment, MET_during_exposure$Drug)
MET_during_exposure$MET_during_exposure.int <- factor(MET_during_exposure$MET_during_exposure.int, levels = c("CTRL.VEH", "CAS.VEH", "CTRL.MET", "CAS.MET"))
```

A1) Time on top
```{r}
anova_TT <- pbad2way(Time_On_Top ~ Treatment + Drug + Treatment:Drug, data = MET_during_exposure, est = "mom", nboot = 5000)
```
-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
postTT <- mcp2a(Time_On_Top ~ Treatment + Drug + Treatment:Drug, data = MET_during_exposure, est = "mom", nboot = 5000)
PTTT <- pairwisePermutationTest(Time_On_Top ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTT
cldList(comparison = PTTT$Comparison, p.value = PTTT$p.adjust, threshold = 0.05)
```

A2) Time on bottom
```{r}
anova_TB <- pbad2way(Time_On_Bottom ~ Treatment + Drug + Treatment:Drug, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTB <- pairwisePermutationTest(Time_On_Bottom ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTB
cldList(comparison = PTTB$Comparison, p.value = PTTB$p.adjust, threshold = 0.05)
```

A3) Absolute turn angle
```{r}
anova_ES <- pbad2way(AbsoluteTurnAngle ~ Treatment + Drug + Treatment:Drug, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTEn <- pairwisePermutationTest(AbsoluteTurnAngle ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTEn
cldList(comparison = PTEn$Comparison, p.value = PTEn$p.adjust, threshold = 0.05)
```

A4) Freezing
```{r}
anova_Fr <- pbad2way(Freezing ~ Treatment + Drug + Treatment:Drug, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFR <- pairwisePermutationTest(Freezing ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFR
cldList(comparison = PTFR$Comparison, p.value = PTFR$p.adjust, threshold = 0.05)
```

A5) Speed
```{r}
anova_Sp <- pbad2way(Speed ~ Treatment + Drug + Treatment:Drug, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTSp <- pairwisePermutationTest(Speed ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTSp
cldList(comparison = PTFrN$Comparison, p.value = PTFrN$p.adjust, threshold = 0.05)
```


* Produce figures on ggplot2
```{r}
metTTplot <- ggplot(MET_during_exposure, aes(x = Drug, y = Time_On_Top, fill = Treatment)) +  geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 360) + scale_fill_brewer(palette = "Set1") + ylab("Time on top third (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 50, label = "a") + annotate("text", x = 1.2, y = 30, label = "b") + annotate("text", x = 1.8, y = 55, label = "a") + annotate("text", x = 2.2, y = 45, label = "b") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed") + theme_bw()

metTBplot <- ggplot(MET_during_exposure, aes(x = Drug, y = Time_On_Bottom, fill = Treatment)) +  geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 380) + scale_fill_brewer(palette = "Set1") +  ylab("Time on bottom third (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 360, label = "a") + annotate("text", x = 1.2, y = 380, label = "b") + annotate("text", x = 1.8, y = 350, label = "a") + annotate("text", x = 2.2, y = 380, label = "ab") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed") + theme_bw()

metESplot <- ggplot(MET_during_exposure, aes(x = Drug, y = AbsoluteTurnAngle, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + ylim(0, 60) + ylab("Absolute Turn Angle (º)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 35, label = "a") + annotate("text", x = 1.2, y = 50, label = "b") + annotate("text", x = 1.8, y = 32, label = "a") + annotate("text", x = 2.2, y = 48, label = "b") + theme_bw()

metFrplot <- ggplot(MET_during_exposure, aes(x = Drug, y = Freezing, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw() + ylab("Freezing (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 35, label = "a") + annotate("text", x = 1.2, y = 105, label = "b") + annotate("text", x = 1.8, y = 40, label = "a") + annotate("text", x = 2.2, y = 100, label = "b")

metSpplot <- ggplot(MET_during_exposure, aes(x = Drug, y = Speed, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw() + ylab("Speed (cm/s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 11, label = "a") + annotate("text", x = 1.2, y = 12, label = "a") + annotate("text", x = 1.8, y = 11, label = "a") + annotate("text", x = 2.2, y = 12, label = "a")

filler <- ggplot() + theme_void() #create blank space to insert image

met_during_plot <- ggarrange(filler, metTTplot, metTBplot, metESplot, metFrplot, metSpplot, ncol = 3, nrow = 2, common.legend = TRUE, labels = c("", "A", "B", "C", "D", "E"), legend = "bottom")

annotate_figure(met_during_plot, top = text_grob("Metergoline treatment - During exposure", color = "black", size = 14))
```

* Run bootstrapped ANOVA for main and interaction effects on 2-way ANOVA, in data for observation AFTER CAS exposure (based on https://rcompanion.org/rcompanion/d_08a.html)

-Create interaction variable for pairwise permutation tests (post-hoc tests); the created variable will be used across all tests
```{r}
MET_after_exposure$MET_after_exposure.int <- interaction(MET_after_exposure$Treatment, MET_after_exposure$Drug)
MET_after_exposure$MET_after_exposure.int <- factor(MET_after_exposure$MET_after_exposure.int, levels = c("CTRL.VEH", "CAS.VEH", "CTRL.MET", "CAS.MET"))
```

B1) Time on top
```{r}
anova_TT2 <- pbad2way(Time_On_Top ~ Treatment + Drug + Treatment:Drug, data = MET_after_exposure, est = "mom", nboot = 5000)
```
-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTT2 <- pairwisePermutationTest(Time_On_Top ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTT2
cldList(comparison = PTTT2$Comparison, p.value = PTTT2$p.adjust, threshold = 0.05)
```

B2) Time on bottom
```{r}
anova_TB2 <- pbad2way(Time_On_Bottom ~ Treatment + Drug + Treatment:Drug, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTB2 <- pairwisePermutationTest(Time_On_Bottom ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTB2
cldList(comparison = PTTB2$Comparison, p.value = PTTB2$p.adjust, threshold = 0.05)
```

B3) Absolute Turn Angle
```{r}
anova_ES2 <- pbad2way(AbsoluteTurnAngle ~ Treatment + Drug + Treatment:Drug, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTES2 <- pairwisePermutationTest(AbsoluteTurnAngle ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTES2
cldList(comparison = PTES2$Comparison, p.value = PTES2$p.adjust, threshold = 0.05)
```

B4) Freezing
```{r}
anova_Fr2 <- pbad2way(Freezing ~ Treatment + Drug + Treatment:Drug, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFr2 <- pairwisePermutationTest(Freezing ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFr2
cldList(comparison = PTFr2$Comparison, p.value = PTFr2$p.adjust, threshold = 0.05)
```

B5) Swimming speed
```{r}
anova_Sp2 <- pbad2way(Speed ~ Treatment + Drug + Treatment:Drug, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTSp2 <- pairwisePermutationTest(Speed ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTSp2
cldList(comparison = PTSp2$Comparison, p.value = PTSp2$p.adjust, threshold = 0.05)
```

* Produce figures on ggplot2
```{r}
metTTplot2 <- ggplot(MET_after_exposure, aes(x = Drug, y = Time_On_Top, fill = Treatment)) +  geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 360) + scale_fill_brewer(palette = "Accent") + theme_bw() + ylab("Time on top third (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 100, label = "a") + annotate("text", x = 1.2, y = 30, label = "a") + annotate("text", x = 1.8, y = 120, label = "a") + annotate("text", x = 2.2, y = 100, label = "a") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

metTBplot2 <- ggplot(MET_after_exposure, aes(x = Drug, y = Time_On_Bottom, fill = Treatment)) +  geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 380) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Time on bottom third (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 360, label = "a") + annotate("text", x = 1.2, y = 380, label = "b") + annotate("text", x = 1.8, y = 380, label = "a") + annotate("text", x = 2.2, y = 310, label = "a") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

metESplot2 <- ggplot(MET_after_exposure, aes(x = Drug, y = AbsoluteTurnAngle, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw() + ylab("Absolute Turn Angle (º)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 35, label = "a") + annotate("text", x = 1.2, y = 40, label = "a") + annotate("text", x = 1.8, y = 40, label = "a") + annotate("text", x = 2.2, y = 35, label = "a") + ylim(0, 60)

metFrplot2 <- ggplot(MET_after_exposure, aes(x = Drug, y = Freezing, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw() + ylab("Freezing (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 35, label = "a") + annotate("text", x = 1.2, y = 90, label = "b") + annotate("text", x = 1.8, y = 40, label = "ac") + annotate("text", x = 2.2, y = 60, label = "bc")

metSpplot2 <- ggplot(MET_after_exposure, aes(x = Drug, y = Speed, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw() + ylab("Speed (cm/s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 13, label = "a") + annotate("text", x = 1.2, y = 13, label = "a") + annotate("text", x = 1.8, y = 13, label = "a") + annotate("text", x = 2.2, y = 13, label = "a")

met_after_plot <- ggarrange(filler, metTTplot2, metTBplot2, metESplot2, metFrplot2, metSpplot2, ncol = 3, nrow = 2, common.legend = TRUE, labels = c("","A", "B", "C", "D", "E"), legend = "bottom")

annotate_figure(met_after_plot, top = text_grob("Metergoline treatment - After exposure", color = "black", size = 14))
```