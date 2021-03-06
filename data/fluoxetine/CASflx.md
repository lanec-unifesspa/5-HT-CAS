---
title: "R Notebook for effects of fluoxetine on zebrafish alarm reactions (LaNeC)"
author:
- Caio Maximino^[Universidade Federal do Sul e Sudeste do Pará]
- Monica Gomes Lima^[Universidade do Estado do Pará]
- Rhayra Xavier do Carmo Silva^[Universidade Federal do Pará]
output:
  github_document 
subtitle: From project "Role of serotonin on behavioral responses to alarm substance in zebrafish: A putative model for panic disorder" (DOI: 10.17605/OSF.IO/BK85D)
tags:
- fear
- zebrafish
- fluoxetine
abstract: |
  Behavioral data during and after exposure to conspecific alarm substance in zebrafish treated with fluoxetine (2.5 and 25 mg/kg)
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
x1 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/5-HT-CAS/master/data/fluoxetine/FLX_during_exposure.csv")
FLX_during_exposure <- read.csv(text = x1)
# Ensure factor variables are interpreted as such
FLX_during_exposure$Treatment <- factor(FLX_during_exposure$Treatment, levels = c("CTRL", "CAS"))
FLX_during_exposure$Dose <- as.factor(FLX_during_exposure$Dose)
View(FLX_during_exposure)

x2 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/5-HT-CAS/master/data/fluoxetine/FLX_after_exposure.csv")
FLX_after_exposure <- read.csv(text = x2)
# Ensure factor variables are interpreted as such
FLX_after_exposure$Treatment <- factor(FLX_after_exposure$Treatment, levels = c("CTRL", "CAS"))
FLX_after_exposure$Dose <- as.factor(FLX_after_exposure$Dose)
FLX_after_exposure$Treatment <- relevel(FLX_after_exposure$Treatment, "CTRL")
View(FLX_after_exposure)
```

* Detect and remove outliers, based on median absolute deviations
```{r}
gA <- FLX_during_exposure[c(1:15), ]
median(gA$Time_On_Bottom) + 3*mad(gA$Time_On_Bottom)
median(gA$Time_On_Bottom) - 3*mad(gA$Time_On_Bottom)
gB <- FLX_during_exposure[c(16:30), ]
median(gB$Time_On_Bottom) + 3*mad(gB$Time_On_Bottom)
median(gB$Time_On_Bottom) - 3*mad(gB$Time_On_Bottom)
gC <- FLX_during_exposure[c(31:45), ]
median(gC$Time_On_Bottom) + 3*mad(gC$Time_On_Bottom)
median(gC$Time_On_Bottom) - 3*mad(gC$Time_On_Bottom)
gD <- FLX_during_exposure[c(46:60), ]
median(gD$Time_On_Bottom) + 3*mad(gD$Time_On_Bottom)
median(gD$Time_On_Bottom) - 3*mad(gD$Time_On_Bottom)
gE <- FLX_during_exposure[c(61:75), ]
median(gE$Time_On_Bottom) + 3*mad(gE$Time_On_Bottom)
median(gE$Time_On_Bottom) - 3*mad(gE$Time_On_Bottom)
#1 outlier found
gF <- FLX_during_exposure[c(76:90), ]
median(gF$Time_On_Bottom) + 3*mad(gF$Time_On_Bottom)
median(gF$Time_On_Bottom) - 3*mad(gF$Time_On_Bottom)
FLX_during_exposure.nooutlier <- FLX_during_exposure[-c(75), ]
View(FLX_during_exposure.nooutlier)

FLX_after_exposure.nooutlier <- FLX_after_exposure[-c(75), ]
View(FLX_after_exposure.nooutlier)
```

* Run bootstrapped ANOVA for main and interaction effects on 2-way ANOVA, in data for observation DURING CAS exposure (based on https://rcompanion.org/rcompanion/d_08a.html)

-Create seed for computational reproducibility
```{r}
set.seed(42)
```

-Create interaction variable for pairwise permutation tests (post-hoc tests); the created variable will be used across all tests
```{r}
FLX_during_exposure.nooutlier$FLX_during_exposure.nooutlier.int <- interaction(FLX_during_exposure.nooutlier$Treatment, FLX_during_exposure.nooutlier$Dose)
FLX_during_exposure.nooutlier$FLX_during_exposure.nooutlier.int <- factor(FLX_during_exposure.nooutlier$FLX_during_exposure.nooutlier.int, levels = c("CTRL.0", "CAS.0", "CTRL.2.5", "CAS.2.5", "CTRL.25", "CAS.25"))
```

A1) Time on top
```{r}
anova_TT <- pbad2way(Time_On_Top ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure.nooutlier, est = "mom", nboot = 5000)
```
-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTT <- pairwisePermutationTest(Time_On_Top ~ FLX_during_exposure.nooutlier.int, data = FLX_during_exposure.nooutlier, est = "mom", nboot = 5000, method = "fdr")
PTTT
cldList(comparison = PTTT$Comparison, p.value = PTTT$p.adjust, threshold = 0.05)
```

A2) Time on bottom
```{r}
anova_TB <- pbad2way(Time_On_Bottom ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure.nooutlier, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTB <- pairwisePermutationTest(Time_On_Bottom ~ FLX_during_exposure.nooutlier.int, data = FLX_during_exposure.nooutlier, est = "mom", nboot = 5000, method = "fdr")
PTTB
cldList(comparison = PTTB$Comparison, p.value = PTTB$p.adjust, threshold = 0.05)
```

A3) Erratic swimming
```{r}
anova_ES <- pbad2way(AbsoluteTurnAngle ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure.nooutlier, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTES <- pairwisePermutationTest(AbsoluteTurnAngle ~ FLX_during_exposure.nooutlier.int, data = FLX_during_exposure.nooutlier, est = "mom", nboot = 5000, method = "fdr")
PTES
cldList(comparison = PTES$Comparison, p.value = PTES$p.adjust, threshold = 0.05)
```

A4) Freezing
```{r}
anova_Fr <- pbad2way(Freezing ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure.nooutlier, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFr <- pairwisePermutationTest(Freezing ~ FLX_during_exposure.nooutlier.int, data = FLX_during_exposure.nooutlier, est = "mom", nboot = 5000, method = "fdr")
PTFr
cldList(comparison = PTFr$Comparison, p.value = PTFr$p.adjust, threshold = 0.05)
```

A5) Locomotion
```{r}
anova_Sp <- pbad2way(Speed ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure.nooutlier, est = "mom", nboot = 5000)
```

* Produce figures on ggplot2
```{r}
blank <- grid.rect(gp=gpar(col="white"))
#Create a blank figure to produce space for a scheme. 

flxTTplot <- ggplot(FLX_during_exposure.nooutlier, aes(x = Dose, y = Time_On_Top, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 360) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Time on top third (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 200, label = "a") + annotate("text", x = 1.2, y = 100, label = "b") + annotate("text", x = 1.8, y = 300, label = "c") + annotate("text", x = 2.2, y = 200, label = "de") + annotate("text", x = 2.8, y = 270, label = "ad") + annotate("text", x = 3.2, y = 280, label = "e") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

flxTBplot <- ggplot(FLX_during_exposure.nooutlier, aes(x = Dose, y = Time_On_Bottom, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 400) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Time on bottom third (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 380, label = "a") + annotate("text", x = 1.2, y = 390, label = "b") + annotate("text", x = 1.8, y = 170, label = "c") + annotate("text", x = 2.2, y = 300, label = "d") + annotate("text", x = 2.8, y = 380, label = "d") + annotate("text", x = 3.2, y = 370, label = "ad") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

flxESplot <- ggplot(FLX_during_exposure.nooutlier, aes(x = Dose, y = AbsoluteTurnAngle, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Absolute Turn Angle (º)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 30, label = "a") + annotate("text", x = 1.2, y = 50, label = "b") + annotate("text", x = 1.8, y = 30, label = "a") + annotate("text", x = 2.2, y = 40, label = "c") + annotate("text", x = 2.8, y = 30, label = "ac") + annotate("text", x = 3.2, y = 30, label = "a")

flxFrplot <- ggplot(FLX_during_exposure.nooutlier, aes(x = Dose, y = Freezing, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Freezing (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 35, label = "a") + annotate("text", x = 1.2, y = 125, label = "b") + annotate("text", x = 1.8, y = 35, label = "a") + annotate("text", x = 2.2, y = 85, label = "c") + annotate("text", x = 2.8, y = 35, label = "a") + annotate("text", x = 3.2, y = 55, label = "c")

flxESplot <- ggplot(FLX_during_exposure.nooutlier, aes(x = Dose, y = AbsoluteTurnAngle, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Absolute Turn Angle (º)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 30, label = "a") + annotate("text", x = 1.2, y = 50, label = "b") + annotate("text", x = 1.8, y = 30, label = "a") + annotate("text", x = 2.2, y = 40, label = "c") + annotate("text", x = 2.8, y = 30, label = "ac") + annotate("text", x = 3.2, y = 30, label = "a")

flxSpplot <- ggplot(FLX_during_exposure.nooutlier, aes(x = Dose, y = Speed, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Speed (cm/s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 11, label = "a") + annotate("text", x = 1.2, y = 11, label = "a") + annotate("text", x = 1.8, y = 11, label = "a") + annotate("text", x = 2.2, y = 11, label = "a") + annotate("text", x = 2.8, y = 11, label = "a") + annotate("text", x = 3.2, y = 11, label = "a")

ggarrange(flxTTplot, flxTBplot, flxESplot, flxFrplot, flxSpplot, ncol = 3, nrow = 3, common.legend = TRUE, labels = c("A", "B", "C", "D", "E"), legend = "bottom")
```

* Run bootstrapped ANOVA for main and interaction effects on 2-way ANOVA, in data for observation AFTER CAS exposure (based on https://rcompanion.org/rcompanion/d_08a.html)

-Create interaction variable for pairwise permutation tests (post-hoc tests); the created variable will be used across all tests
```{r}
FLX_after_exposure.nooutlier$FLX_after_exposure.nooutlier.int <- interaction(FLX_after_exposure.nooutlier$Treatment, FLX_after_exposure.nooutlier$Dose)
FLX_after_exposure.nooutlier$FLX_after_exposure.nooutlier.int <- factor(FLX_after_exposure.nooutlier$FLX_after_exposure.nooutlier.int, levels = c("CTRL.0", "CAS.0", "CTRL.2.5", "CAS.2.5", "CTRL.25", "CAS.25"))
```

B1) Time on top
```{r}
anova_TT_after <- pbad2way(Time_On_Top ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure.nooutlier, est = "mom", nboot = 5000)
```
-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTTafter <- pairwisePermutationTest(Time_On_Top ~ FLX_after_exposure.nooutlier.int, data = FLX_after_exposure.nooutlier, est = "mom", nboot = 5000, method = "fdr")
PTTTafter
cldList(comparison = PTTTafter$Comparison, p.value = PTTTafter$p.adjust, threshold = 0.05)
```

B2) Time on bottom
```{r}
anova_TBafter <- pbad2way(Time_On_Bottom ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure.nooutlier, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTBafter <- pairwisePermutationTest(Time_On_Bottom ~ FLX_after_exposure.nooutlier.int, data = FLX_after_exposure.nooutlier, est = "mom", nboot = 5000, method = "fdr")
PTTBafter
cldList(comparison = PTTBafter$Comparison, p.value = PTTBafter$p.adjust, threshold = 0.05)
```

B3) Erratic swimming
```{r}
anova_ESafter <- pbad2way(AbsoluteTurnAngle ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure.nooutlier, est = "mom", nboot = 5000)
```

B4) Freezing
```{r}
anova_Frafter <- pbad2way(Freezing ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure.nooutlier, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFrafter <- pairwisePermutationTest(Freezing ~ FLX_after_exposure.nooutlier.int, data = FLX_after_exposure.nooutlier, est = "mom", nboot = 5000, method = "fdr")
PTFrafter
cldList(comparison = PTFrafter$Comparison, p.value = PTFrafter$p.adjust, threshold = 0.05)
```

B5) Locomotion
```{r}
anova_Spafter <- pbad2way(Speed ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure.nooutlier, est = "mom", nboot = 5000)
```

* Produce figures on ggplot2
```{r}
flxTTplot2 <- ggplot(FLX_after_exposure.nooutlier, aes(x = Dose, y = Time_On_Top, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 400) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Time on top third (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 120, label = "a") + annotate("text", x = 1.2, y = 50, label = "b") + annotate("text", x = 1.8, y = 160, label = "ac") + annotate("text", x = 2.2, y = 370, label = "abc") + annotate("text", x = 2.8, y = 130, label = "c") + annotate("text", x = 3.2, y = 130, label = "ac") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

flxTBplot2 <- ggplot(FLX_after_exposure.nooutlier, aes(x = Dose, y = Time_On_Bottom, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 400) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Time on bottom third (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 320, label = "a") + annotate("text", x = 1.2, y = 380, label = "b") + annotate("text", x = 1.8, y = 370, label = "a") + annotate("text", x = 2.2, y = 380, label = "ab") + annotate("text", x = 2.8, y = 330, label = "a") + annotate("text", x = 3.2, y = 350, label = "a") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

flxESplot2 <- ggplot(FLX_after_exposure.nooutlier, aes(x = Dose, y = AbsoluteTurnAngle, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Absolute Turn Angle (º)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 65, label = "a") + annotate("text", x = 1.2, y = 65, label = "a") + annotate("text", x = 1.8, y = 60, label = "a") + annotate("text", x = 2.2, y = 55, label = "a") + annotate("text", x = 2.8, y = 62, label = "a") + annotate("text", x = 3.2, y = 58, label = "a")

flxFrplot2 <- ggplot(FLX_after_exposure.nooutlier, aes(x = Dose, y = Freezing, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Freezing (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 26, label = "a") + annotate("text", x = 1.2, y = 60, label = "b") + annotate("text", x = 1.8, y = 12, label = "c") + annotate("text", x = 2.2, y = 58, label = "b") + annotate("text", x = 2.8, y = 12, label = "c") + annotate("text", x = 3.2, y = 56, label = "b")

flxSpplot2 <- ggplot(FLX_after_exposure.nooutlier, aes(x = Dose, y = Speed, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Speed (cm/s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 15, label = "a") + annotate("text", x = 1.2, y = 15, label = "a") + annotate("text", x = 1.8, y = 15, label = "a") + annotate("text", x = 2.2, y = 15, label = "a") + annotate("text", x = 2.8, y = 15, label = "a") + annotate("text", x = 3.2, y = 15, label = "a")

ggarrange(blank, flxTTplot2, flxTBplot2, flxESplot2, flxFrplot2, flxSpplot, ncol = 3, nrow = 2, common.legend = TRUE, labels = c("", "A", "B", "C", "E"), legend = "bottom")
```