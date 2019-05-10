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
MET_during_exposure$Treatment <- as.factor(MET_during_exposure$Treatment)
MET_during_exposure$Dose <- as.factor(MET_during_exposure$Dose)
relevel(MET_during_exposure$Treatment, "CTRL")
View(MET_during_exposure)

x2 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/5-HT-CAS/master/data/metergoline/MET_after_exposure.csv")
MET_after_exposure <- read.csv(text = x2)
# Ensure factor variables are interpreted as such
MET_after_exposure$Treatment <- as.factor(MET_after_exposure$Treatment)
MET_after_exposure$Dose <- as.factor(MET_after_exposure$Dose)
MET_after_exposure$Treatment <- relevel(MET_after_exposure$Treatment, "CTRL")
View(MET_after_exposure)
```

* Detect and remove outliers, based on median absolute deviations
```{r}
gA <- MET_during_exposure[c(1:15), ]
median(gA$TB) + 3*mad(gA$TB)
median(gA$TB) - 3*mad(gA$TB)
gB <- MET_during_exposure[c(16:30), ]
median(gB$TB) + 3*mad(gB$TB)
median(gB$TB) - 3*mad(gB$TB)
gC <- MET_during_exposure[c(31:45), ]
median(gC$TB) + 3*mad(gC$TB)
median(gC$TB) - 3*mad(gC$TB)
gD <- MET_during_exposure[c(46:60), ]
median(gD$TB) + 3*mad(gD$TB)
median(gD$TB) - 3*mad(gD$TB)

#Subject 4 from group CTRL+VEH detected as outlier; remove
MET_during_exposure <- MET_during_exposure[-c(4), ]
MET_after_exposure <- MET_after_exposure[-c(4), ]
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
MET_during_exposure$MET_during_exposure.int <- interaction(MET_during_exposure$Treatment, MET_during_exposure$Dose)
MET_during_exposure$MET_during_exposure.int <- factor(MET_during_exposure$MET_during_exposure.int, levels = c("CTRL.0", "CAS.0", "CTRL.1", "CAS.1"))
```

A1) Time on top
```{r}
anova_TT <- pbad2way(TT ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
```
-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
postTT <- mcp2a(TT ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
PTTT <- pairwisePermutationTest(TT ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTT
cldList(comparison = PTTT$Comparison, p.value = PTTT$p.adjust, threshold = 0.05)
```

A2) Time on bottom
```{r}
anova_TB <- pbad2way(TB ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTB <- pairwisePermutationTest(TB ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTB
cldList(comparison = PTTB$Comparison, p.value = PTTB$p.adjust, threshold = 0.05)
```

A3) Transitions to top
```{r}
anova_Entries <- pbad2way(Entries ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTEn <- pairwisePermutationTest(Entries ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTEn
cldList(comparison = PTEn$Comparison, p.value = PTEn$p.adjust, threshold = 0.05)
```

A4) Erratic swimming
```{r}
anova_ES <- pbad2way(NE ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTES <- pairwisePermutationTest(NE ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTES
cldList(comparison = PTES$Comparison, p.value = PTES$p.adjust, threshold = 0.05)
```

A5) Freezing frequency
```{r}
anova_FrN <- pbad2way(Freeze.N ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFrN <- pairwisePermutationTest(Freeze.N ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFrN
cldList(comparison = PTFrN$Comparison, p.value = PTFrN$p.adjust, threshold = 0.05)
```

A6) Freezing duration
```{r}
anova_FrS <- pbad2way(Freeze.Dur ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFrS <- pairwisePermutationTest(Freeze.Dur ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFrS
cldList(comparison = PTFrS$Comparison, p.value = PTFrS$p.adjust, threshold = 0.05)
```

A7) Thrashing
```{r}
anova_Thr <- pbad2way(Thrashing ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTThr <- pairwisePermutationTest(Thrashing ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTThr
cldList(comparison = PTThr$Comparison, p.value = PTThr$p.adjust, threshold = 0.05)
```

A8) Homebase time
```{r}
anova_HB <- pbad2way(Homebase ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTHB <- pairwisePermutationTest(Homebase ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTHB
cldList(comparison = PTHB$Comparison, p.value = PTHB$p.adjust, threshold = 0.05)
```

A9) Total locomotion
```{r}
anova_Sq <- pbad2way(Quad ~ Treatment + Dose + Treatment:Dose, data = MET_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTSq <- pairwisePermutationTest(Quad ~ MET_during_exposure.int, data = MET_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTSq
cldList(comparison = PTSq$Comparison, p.value = PTSq$p.adjust, threshold = 0.05)
```

* Produce figures on ggplot2
```{r}
metTTplot <- ggplot(MET_during_exposure, aes(x = Dose, y = TT, fill = Treatment)) +  geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 360) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Time on top third (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 50, label = "a") + annotate("text", x = 1.2, y = 30, label = "b") + annotate("text", x = 1.8, y = 55, label = "a") + annotate("text", x = 2.2, y = 45, label = "b") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

metTBplot <- ggplot(MET_during_exposure, aes(x = Dose, y = TB, fill = Treatment)) +  geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 380) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Time on bottom third (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 360, label = "a") + annotate("text", x = 1.2, y = 380, label = "b") + annotate("text", x = 1.8, y = 350, label = "a") + annotate("text", x = 2.2, y = 380, label = "ab") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

metEnplot <- ggplot(MET_during_exposure, aes(x = Dose, y = Entries, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Transitions to top (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 13, label = "a") + annotate("text", x = 1.2, y = 7, label = "b") + annotate("text", x = 1.8, y = 15, label = "a") + annotate("text", x = 2.2, y = 8, label = "b")

metESplot <- ggplot(MET_during_exposure, aes(x = Dose, y = NE, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Erratic swimming (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 5, label = "a") + annotate("text", x = 1.2, y = 15, label = "b") + annotate("text", x = 1.8, y = 7, label = "c") + annotate("text", x = 2.2, y = 18, label = "ac")

metFrNplot <- ggplot(MET_during_exposure, aes(x = Dose, y = Freeze.N, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Freezing frequency (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 9, label = "a") + annotate("text", x = 1.2, y = 12, label = "a") + annotate("text", x = 1.8, y = 25, label = "a") + annotate("text", x = 2.2, y = 5, label = "a")

metFrSplot <- ggplot(MET_during_exposure, aes(x = Dose, y = Freeze.Dur, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Freezing duration (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 55, label = "a") + annotate("text", x = 1.2, y = 350, label = "b") + annotate("text", x = 1.8, y = 350, label = "b") + annotate("text", x = 2.2, y = 360, label = "b")

metThrplot <- ggplot(MET_during_exposure, aes(x = Dose, y = Thrashing, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Thrashing (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 25, label = "a") + annotate("text", x = 1.2, y = 40, label = "a") + annotate("text", x = 1.8, y = 14, label = "a") + annotate("text", x = 2.2, y = 31, label = "a")

metHBplot <- ggplot(MET_during_exposure, aes(x = Dose, y = Homebase, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 100) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Time spent in the homebase (%)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 45, label = "a") + annotate("text", x = 1.2, y = 45, label = "a") + annotate("text", x = 1.8, y = 45, label = "a") + annotate("text", x = 2.2, y = 55, label = "a") + geom_hline(yintercept = 11.11111, color = 'coral', linetype = "dashed")

metSqplot <- ggplot(MET_during_exposure, aes(x = Dose, y = Quad, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Squares crossed (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 360, label = "a") + annotate("text", x = 1.2, y = 300, label = "a") + annotate("text", x = 1.8, y = 270, label = "a") + annotate("text", x = 2.2, y = 400, label = "a")

met_during_plot <- ggarrange(metTTplot, metTBplot, metEnplot, metESplot, metFrSplot, metFrNplot, metThrplot, metHBplot, metSqplot, ncol = 3, nrow = 3, common.legend = TRUE, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), legend = "bottom")

annotate_figure(met_during_plot, top = text_grob("Metergoline treatment - During exposure", color = "black", size = 14))
```

* Run bootstrapped ANOVA for main and interaction effects on 2-way ANOVA, in data for observation AFTER CAS exposure (based on https://rcompanion.org/rcompanion/d_08a.html)

-Create interaction variable for pairwise permutation tests (post-hoc tests); the created variable will be used across all tests
```{r}
MET_after_exposure$MET_after_exposure.int <- interaction(MET_after_exposure$Treatment, MET_after_exposure$Dose)
MET_after_exposure$MET_after_exposure.int <- factor(MET_after_exposure$MET_after_exposure.int, levels = c("CTRL.0", "CAS.0", "CTRL.1", "CAS.1"))
```

B1) Time on top
```{r}
anova_TT2 <- pbad2way(TT ~ Treatment + Dose + Treatment:Dose, data = MET_after_exposure, est = "mom", nboot = 5000)
```
-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTT2 <- pairwisePermutationTest(TT ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTT2
cldList(comparison = PTTT2$Comparison, p.value = PTTT2$p.adjust, threshold = 0.05)
```

B2) Time on bottom
```{r}
anova_TB2 <- pbad2way(TB ~ Treatment + Dose + Treatment:Dose, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTB2 <- pairwisePermutationTest(TB ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTB2
cldList(comparison = PTTB2$Comparison, p.value = PTTB2$p.adjust, threshold = 0.05)
```

B3) Transitions to top
```{r}
anova_Entries2 <- pbad2way(Entries ~ Treatment + Dose + Treatment:Dose, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTEn2 <- pairwisePermutationTest(Entries ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTEn2
cldList(comparison = PTEn2$Comparison, p.value = PTEn2$p.adjust, threshold = 0.05)
```

B4) Erratic swimming
```{r}
anova_ES2 <- pbad2way(NE ~ Treatment + Dose + Treatment:Dose, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTES2 <- pairwisePermutationTest(NE ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTES2
cldList(comparison = PTES2$Comparison, p.value = PTES2$p.adjust, threshold = 0.05)
```

B5) Freezing frequency
```{r}
anova_FrN2 <- pbad2way(Freeze.N ~ Treatment + Dose + Treatment:Dose, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFrN2 <- pairwisePermutationTest(Freeze.N ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFrN2
cldList(comparison = PTFrN2$Comparison, p.value = PTFrN2$p.adjust, threshold = 0.05)
```

B6) Freezing duration
```{r}
anova_FrS2 <- pbad2way(Freeze.Dur ~ Treatment + Dose + Treatment:Dose, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFrS2 <- pairwisePermutationTest(Freeze.Dur ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFrS2
cldList(comparison = PTFrS2$Comparison, p.value = PTFrS2$p.adjust, threshold = 0.05)
```

B7) Thrashing
```{r}
anova_Thr2 <- pbad2way(Thrashing ~ Treatment + Dose + Treatment:Dose, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTThr2 <- pairwisePermutationTest(Thrashing ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTThr2
cldList(comparison = PTThr2$Comparison, p.value = PTThr2$p.adjust, threshold = 0.05)
```

B8) Homebase time
```{r}
anova_HB2 <- pbad2way(Homebase ~ Treatment + Dose + Treatment:Dose, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTHB2 <- pairwisePermutationTest(Homebase ~ MET_after_exposure.int, data = MET_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTHB2
cldList(comparison = PTHB2$Comparison, p.value = PTHB2$p.adjust, threshold = 0.05)
```

B9) Total locomotion
```{r}
anova_Sq2 <- pbad2way(Quad ~ Treatment + Dose + Treatment:Dose, data = MET_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTSq2 <- pairwisePermutationTest(Squares ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTSq2
cldList(comparison = PTSq2$Comparison, p.value = PTSq2$p.adjust, threshold = 0.05)
```

* Produce figures on ggplot2
```{r}
metTTplot2 <- ggplot(MET_after_exposure, aes(x = Dose, y = TT, fill = Treatment)) +  geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 360) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Time on top third (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 100, label = "a") + annotate("text", x = 1.2, y = 30, label = "a") + annotate("text", x = 1.8, y = 120, label = "a") + annotate("text", x = 2.2, y = 100, label = "a") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

metTBplot2 <- ggplot(MET_after_exposure, aes(x = Dose, y = TB, fill = Treatment)) +  geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 380) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Time on bottom third (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 360, label = "a") + annotate("text", x = 1.2, y = 380, label = "b") + annotate("text", x = 1.8, y = 380, label = "ab") + annotate("text", x = 2.2, y = 310, label = "a") + geom_hline(yintercept = 120, color = 'coral', linetype = "dashed")

metEnplot2 <- ggplot(MET_after_exposure, aes(x = Dose, y = Entries, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Transitions to top (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 85, label = "ab") + annotate("text", x = 1.2, y = 30, label = "a") + annotate("text", x = 1.8, y = 45, label = "ab") + annotate("text", x = 2.2, y = 70, label = "b")

metESplot2 <- ggplot(MET_after_exposure, aes(x = Dose, y = NE, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Erratic swimming (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 7, label = "a") + annotate("text", x = 1.2, y = 11, label = "a") + annotate("text", x = 1.8, y = 12, label = "a") + annotate("text", x = 2.2, y = 8, label = "a") + ylim(0, 13)

metFrNplot2 <- ggplot(MET_after_exposure, aes(x = Dose, y = Freeze.N, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Freezing frequency (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 5, label = "a") + annotate("text", x = 1.2, y = 4, label = "a") + annotate("text", x = 1.8, y = 6, label = "a") + annotate("text", x = 2.2, y = 3, label = "a")

metFrSplot2 <- ggplot(MET_after_exposure, aes(x = Dose, y = Freeze.Dur, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Freezing duration (s)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 110, label = "a") + annotate("text", x = 1.2, y = 350, label = "b") + annotate("text", x = 1.8, y = 360, label = "ab") + annotate("text", x = 2.2, y = 360, label = "ab")

metThrplot2 <- ggplot(MET_after_exposure, aes(x = Dose, y = Thrashing, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Thrashing (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 7, label = "a") + annotate("text", x = 1.2, y = 20, label = "a") + annotate("text", x = 1.8, y = 14, label = "a") + annotate("text", x = 2.2, y = 7, label = "a")

metHBplot2 <- ggplot(MET_after_exposure, aes(x = Dose, y = Homebase, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 100) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Time spent in the homebase (%)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 53, label = "a") + annotate("text", x = 1.2, y = 75, label = "b") + annotate("text", x = 1.8, y = 55, label = "a") + annotate("text", x = 2.2, y = 55, label = "a") + geom_hline(yintercept = 11.11111, color = 'coral', linetype = "dashed")

metSqplot2 <- ggplot(MET_after_exposure, aes(x = Dose, y = Quad, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Squares crossed (N)") + xlab("Metergoline dose (mg/kg)") + annotate("text", x = 0.8, y = 530, label = "a") + annotate("text", x = 1.2, y = 480, label = "a") + annotate("text", x = 1.8, y = 350, label = "a") + annotate("text", x = 2.2, y = 450, label = "a")

met_after_plot <- ggarrange(metTTplot2, metTBplot2, metEnplot2, metESplot2, metFrSplot2, metFrNplot2, metThrplot2, metHBplot2, metSqplot2, ncol = 3, nrow = 3, common.legend = TRUE, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), legend = "bottom")

annotate_figure(met_after_plot, top = text_grob("Metergoline treatment - After exposure", color = "black", size = 14))
```