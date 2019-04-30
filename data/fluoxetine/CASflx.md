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
FLX_during_exposure$Treatment <- as.factor(FLX_during_exposure$Treatment)
FLX_during_exposure$Dose <- as.factor(FLX_during_exposure$Dose)
View(FLX_during_exposure)

x2 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/5-HT-CAS/master/data/fluoxetine/FLX_after_exposure.csv")
FLX_after_exposure <- read.csv(text = x2)
# Ensure factor variables are interpreted as such
FLX_after_exposure$Treatment <- as.factor(FLX_after_exposure$Treatment)
FLX_after_exposure$Dose <- as.factor(FLX_after_exposure$Dose)
FLX_after_exposure$Treatment <- relevel(FLX_after_exposure$Treatment, "CTRL")
View(FLX_after_exposure)
```

* Detect and remove outliers, based on median absolute deviations
```{r}
gA <- FLX_during_exposure[c(1:15), ]
median(gA$TB) + 3*mad(gA$TB)
median(gA$TB) - 3*mad(gA$TB)
gB <- durante_exposicao[c(16:30), ]
median(gB$TB) + 3*mad(gB$TB)
median(gB$TB) - 3*mad(gB$TB)
gC <- durante_exposicao[c(31:45), ]
median(gC$TB) + 3*mad(gC$TB)
median(gC$TB) - 3*mad(gC$TB)
gD <- durante_exposicao[c(46:60), ]
median(gD$TB) + 3*mad(gD$TB)
median(gD$TB) - 3*mad(gD$TB)
gE <- durante_exposicao[c(61:75), ]
median(gE$TB) + 3*mad(gE$TB)
median(gE$TB) - 3*mad(gE$TB)
gF <- durante_exposicao[c(76:90), ]
median(gF$TB) + 3*mad(gF$TB)
median(gF$TB) - 3*mad(gF$TB)

```

* Run bootstrapped ANOVA for main and interaction effects on 2-way ANOVA, in data for observation DURING CAS exposure (based on https://rcompanion.org/rcompanion/d_08a.html)

-Create seed for computational reproducibility
```{r}
set.seed(42)
```

-Create interaction variable for pairwise permutation tests (post-hoc tests); the created variable will be used across all tests
```{r}
FLX_during_exposure$FLX_during_exposure.int <- interaction(FLX_during_exposure$Treatment, FLX_during_exposure$Dose)
FLX_during_exposure$FLX_during_exposure.int <- factor(FLX_during_exposure$FLX_during_exposure.int, levels = c("CTRL.0", "CAS.0", "CTRL.2.5", "CAS.2.5", "CTRL.25", "CAS.25"))
```

A1) Time on top
```{r}
anova_TT <- pbad2way(TT ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
```
-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
postTT <- mcp2a(TT ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
PTTT <- pairwisePermutationTest(TT ~ FLX_during_exposure.int, data = FLX_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTT
cldList(comparison = PTTT$Comparison, p.value = PTTT$p.adjust, threshold = 0.05)
```

A2) Time on bottom
```{r}
anova_TB <- pbad2way(TB ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTB <- pairwisePermutationTest(TB ~ FLX_during_exposure.int, data = FLX_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTB
cldList(comparison = PTTB$Comparison, p.value = PTTB$p.adjust, threshold = 0.05)
```

A3) Transitions to top
```{r}
anova_Entries <- pbad2way(Entries ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTEn <- pairwisePermutationTest(Entries ~ FLX_during_exposure.int, data = FLX_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTEn
cldList(comparison = PTEn$Comparison, p.value = PTEn$p.adjust, threshold = 0.05)
```

A4) Erratic swimming
```{r}
anova_ES <- pbad2way(Erratic_swimming ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTES <- pairwisePermutationTest(Erratic_swimming ~ FLX_during_exposure.int, data = FLX_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTES
cldList(comparison = PTES$Comparison, p.value = PTES$p.adjust, threshold = 0.05)
```

A5) Freezing frequency
```{r}
anova_FrN <- pbad2way(FreezingN ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFrN <- pairwisePermutationTest(FreezingN ~ FLX_during_exposure.int, data = FLX_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFrN
cldList(comparison = PTFrN$Comparison, p.value = PTFrN$p.adjust, threshold = 0.05)
```

A6) Freezing duration
```{r}
anova_FrS <- pbad2way(FreezingS ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFrS <- pairwisePermutationTest(FreezingS ~ FLX_during_exposure.int, data = FLX_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFrS
cldList(comparison = PTFrS$Comparison, p.value = PTFrS$p.adjust, threshold = 0.05)
```

A7) Thrashing
```{r}
anova_Thr <- pbad2way(Thrashing ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTThr <- pairwisePermutationTest(Thrashing ~ FLX_during_exposure.int, data = FLX_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTThr
cldList(comparison = PTThr$Comparison, p.value = PTThr$p.adjust, threshold = 0.05)
```

A8) Homebase time
```{r}
anova_HB <- pbad2way(Homebase ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTHB <- pairwisePermutationTest(Homebase ~ FLX_during_exposure.int, data = FLX_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTHB
cldList(comparison = PTHB$Comparison, p.value = PTHB$p.adjust, threshold = 0.05)
```

A9) Total locomotion
```{r}
anova_Sq <- pbad2way(Squares ~ Treatment + Dose + Treatment:Dose, data = FLX_during_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTSq <- pairwisePermutationTest(Squares ~ FLX_during_exposure.int, data = FLX_during_exposure, est = "mom", nboot = 5000, method = "fdr")
PTSq
cldList(comparison = PTSq$Comparison, p.value = PTSq$p.adjust, threshold = 0.05)
```

* Produce figures on ggplot2
```{r}
flxTTplot <- ggplot(durante_exposicao, aes(x = Dose, y = TT, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 360) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Time on top third (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 200, label = "a") + annotate("text", x = 1.2, y = 100, label = "b") + annotate("text", x = 1.8, y = 300, label = "c") + annotate("text", x = 2.2, y = 200, label = "ad") + annotate("text", x = 2.8, y = 270, label = "abd") + annotate("text", x = 3.2, y = 280, label = "ad")

flxTBplot <- ggplot(durante_exposicao, aes(x = Dose, y = TB, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 400) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Time on bottom third (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 380, label = "c") + annotate("text", x = 1.2, y = 390, label = "d") + annotate("text", x = 1.8, y = 170, label = "a") + annotate("text", x = 2.2, y = 300, label = "b") + annotate("text", x = 2.8, y = 380, label = "ab") + annotate("text", x = 3.2, y = 370, label = "bc")

flxEnplot <- ggplot(durante_exposicao, aes(x = Dose, y = Entries, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Transitions to top (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 130, label = "ab") + annotate("text", x = 1.2, y = 70, label = "c") + annotate("text", x = 1.8, y = 100, label = "a") + annotate("text", x = 2.2, y = 80, label = "b") + annotate("text", x = 2.8, y = 160, label = "abc") + annotate("text", x = 3.2, y = 70, label = "bc")

flxESplot <- ggplot(durante_exposicao, aes(x = Dose, y = Erratic_swimming, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Erratic swimming (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 12, label = "ab") + annotate("text", x = 1.2, y = 32, label = "c") + annotate("text", x = 1.8, y = 5, label = "a") + annotate("text", x = 2.2, y = 15, label = "b") + annotate("text", x = 2.8, y = 5, label = "abc") + annotate("text", x = 3.2, y = 14, label = "bc")

flxESplot <- ggplot(durante_exposicao, aes(x = Dose, y = Erratic_swimming, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Erratic swimming (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 12, label = "ab") + annotate("text", x = 1.2, y = 32, label = "cd") + annotate("text", x = 1.8, y = 5, label = "a") + annotate("text", x = 2.2, y = 15, label = "bc") + annotate("text", x = 2.8, y = 5, label = "a") + annotate("text", x = 3.2, y = 14, label = "d")

flxFrNplot <- ggplot(durante_exposicao, aes(x = Dose, y = FreezingN, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Freezing frequency (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 8, label = "ab") + annotate("text", x = 1.2, y = 10, label = "ac") + annotate("text", x = 1.8, y = 5, label = "bd") + annotate("text", x = 2.2, y = 20, label = "abcd") + annotate("text", x = 2.8, y = 4, label = "d") + annotate("text", x = 3.2, y = 35, label = "c")

flxFrSplot <- ggplot(durante_exposicao, aes(x = Dose, y = FreezingS, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 380) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Freezing duration (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 370, label = "a") + annotate("text", x = 1.2, y = 360, label = "a") + annotate("text", x = 1.8, y = 52, label = "a") + annotate("text", x = 2.2, y = 60, label = "a") + annotate("text", x = 2.8, y = 370, label = "a") + annotate("text", x = 3.2, y = 90, label = "a")

flxThrplot <- ggplot(durante_exposicao, aes(x = Dose, y = Thrashing, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Thrashing (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 25, label = "ab") + annotate("text", x = 1.2, y = 58, label = "c") + annotate("text", x = 1.8, y = 10, label = "a") + annotate("text", x = 2.2, y = 31, label = "bc") + annotate("text", x = 2.8, y = 13, label = "a") + annotate("text", x = 3.2, y = 35, label = "c")

flxHBplot <- ggplot(durante_exposicao, aes(x = Dose, y = Homebase, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 100) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Time spent in the homebase (%)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 60, label = "a") + annotate("text", x = 1.2, y = 60, label = "a") + annotate("text", x = 1.8, y = 45, label = "a") + annotate("text", x = 2.2, y = 55, label = "a") + annotate("text", x = 2.8, y = 60, label = "a") + annotate("text", x = 3.2, y = 66, label = "a")

flxSqplot <- ggplot(durante_exposicao, aes(x = Dose, y = Squares, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Squares crossed (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 200, label = "ab") + annotate("text", x = 1.2, y = 350, label = "ab") + annotate("text", x = 1.8, y = 250, label = "a") + annotate("text", x = 2.2, y = 200, label = "a") + annotate("text", x = 2.8, y = 370, label = "b") + annotate("text", x = 3.2, y = 210, label = "ab")

ggarrange(flxTTplot, flxTBplot, flxEnplot, flxESplot, flxFrSplot, flxFrNplot, flxThrplot, flxHBplot, flxSqplot, ncol = 3, nrow = 3, common.legend = TRUE, labels = c("A", "B", "C", "E", "F", "G", "H", "I"), legend = "bottom")
```

* Run bootstrapped ANOVA for main and interaction effects on 2-way ANOVA, in data for observation AFTER CAS exposure (based on https://rcompanion.org/rcompanion/d_08a.html)

-Create interaction variable for pairwise permutation tests (post-hoc tests); the created variable will be used across all tests
```{r}
FLX_after_exposure$FLX_after_exposure.int <- interaction(FLX_after_exposure$Treatment, FLX_after_exposure$Dose)
FLX_after_exposure$FLX_after_exposure.int <- factor(FLX_after_exposure$FLX_after_exposure.int, levels = c("CTRL.0", "CAS.0", "CTRL.2.5", "CAS.2.5", "CTRL.25", "CAS.25"))
```

B1) Time on top
```{r}
anova_TT2 <- pbad2way(TT ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure, est = "mom", nboot = 5000)
```
-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTT2 <- pairwisePermutationTest(TT ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTT2
cldList(comparison = PTTT2$Comparison, p.value = PTTT2$p.adjust, threshold = 0.05)
```

B2) Time on bottom
```{r}
anova_TB2 <- pbad2way(TB ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTTB2 <- pairwisePermutationTest(TB ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTTB2
cldList(comparison = PTTB2$Comparison, p.value = PTTB2$p.adjust, threshold = 0.05)
```

B3) Transitions to top
```{r}
anova_Entries2 <- pbad2way(Entries ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTEn2 <- pairwisePermutationTest(Entries ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTEn2
cldList(comparison = PTEn2$Comparison, p.value = PTEn2$p.adjust, threshold = 0.05)
```

B4) Erratic swimming
```{r}
anova_ES2 <- pbad2way(Erratic_swimming ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTES2 <- pairwisePermutationTest(Erratic_swimming ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTES2
cldList(comparison = PTES2$Comparison, p.value = PTES2$p.adjust, threshold = 0.05)
```

B5) Freezing frequency
```{r}
anova_FrN2 <- pbad2way(FreezingN ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFrN2 <- pairwisePermutationTest(FreezingN ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFrN2
cldList(comparison = PTFrN2$Comparison, p.value = PTFrN2$p.adjust, threshold = 0.05)
```

B6) Freezing duration
```{r}
anova_FrS2 <- pbad2way(FreezingS ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTFrS2 <- pairwisePermutationTest(FreezingS ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTFrS2
cldList(comparison = PTFrS2$Comparison, p.value = PTFrS2$p.adjust, threshold = 0.05)
```

B7) Thrashing
```{r}
anova_Thr2 <- pbad2way(Thrashing ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTThr2 <- pairwisePermutationTest(Thrashing ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTThr2
cldList(comparison = PTThr2$Comparison, p.value = PTThr2$p.adjust, threshold = 0.05)
```

B8) Homebase time
```{r}
anova_HB2 <- pbad2way(Homebase ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTHB2 <- pairwisePermutationTest(Homebase ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTHB2
cldList(comparison = PTHB2$Comparison, p.value = PTHB2$p.adjust, threshold = 0.05)
```

B9) Total locomotion
```{r}
anova_Sq2 <- pbad2way(Squares ~ Treatment + Dose + Treatment:Dose, data = FLX_after_exposure, est = "mom", nboot = 5000)
```

-Post-hoc tests (pairwise permutation tests on modified M-estimators, with 5000 bootstrap samples)
```{r}
PTSq2 <- pairwisePermutationTest(Squares ~ FLX_after_exposure.int, data = FLX_after_exposure, est = "mom", nboot = 5000, method = "fdr")
PTSq2
cldList(comparison = PTSq2$Comparison, p.value = PTSq2$p.adjust, threshold = 0.05)
```

* Produce figures on ggplot2
```{r}
flxTTplot2 <- ggplot(FLX_after_exposure, aes(x = Dose, y = TT, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 400) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Time on top third (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 120, label = "a") + annotate("text", x = 1.2, y = 50, label = "b") + annotate("text", x = 1.8, y = 160, label = "ac") + annotate("text", x = 2.2, y = 370, label = "abc") + annotate("text", x = 2.8, y = 120, label = "c") + annotate("text", x = 3.2, y = 120, label = "ac")

flxTBplot2 <- ggplot(FLX_after_exposure, aes(x = Dose, y = TB, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 400) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Time on bottom third (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 320, label = "a") + annotate("text", x = 1.2, y = 380, label = "b") + annotate("text", x = 1.8, y = 370, label = "a") + annotate("text", x = 2.2, y = 380, label = "ab") + annotate("text", x = 2.8, y = 330, label = "a") + annotate("text", x = 3.2, y = 350, label = "a")

flxEnplot2 <- ggplot(FLX_after_exposure, aes(x = Dose, y = Entries, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Transitions to top (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 90, label = "a") + annotate("text", x = 1.2, y = 30, label = "b") + annotate("text", x = 1.8, y = 60, label = "a") + annotate("text", x = 2.2, y = 50, label = "ab") + annotate("text", x = 2.8, y = 45, label = "a") + annotate("text", x = 3.2, y = 70, label = "a")

flxESplot <- ggplot(durante_exposicao, aes(x = Dose, y = Erratic_swimming, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 10) + ylab("Erratic swimming (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 12, label = "ab") + annotate("text", x = 1.2, y = 32, label = "c") + annotate("text", x = 1.8, y = 5, label = "a") + annotate("text", x = 2.2, y = 15, label = "b") + annotate("text", x = 2.8, y = 5, label = "abc") + annotate("text", x = 3.2, y = 14, label = "bc")

flxESplot2 <- ggplot(FLX_after_exposure, aes(x = Dose, y = Erratic_swimming, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Erratic swimming (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 10, label = "a") + annotate("text", x = 1.2, y = 10, label = "a") + annotate("text", x = 1.8, y = 8, label = "a") + annotate("text", x = 2.2, y = 10, label = "a") + annotate("text", x = 2.8, y = 5, label = "a") + annotate("text", x = 3.2, y = 7, label = "a") + ylim(0, 13)

flxFrNplot2 <- ggplot(FLX_after_exposure, aes(x = Dose, y = FreezingN, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Freezing frequency (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 12, label = "a") + annotate("text", x = 1.2, y = 8, label = "a") + annotate("text", x = 1.8, y = 5, label = "a") + annotate("text", x = 2.2, y = 8, label = "a") + annotate("text", x = 2.8, y = 4, label = "a") + annotate("text", x = 3.2, y = 5, label = "a")

flxFrSplot2 <- ggplot(FLX_after_exposure, aes(x = Dose, y = FreezingS, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 380) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Freezing duration (s)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 350, label = "a") + annotate("text", x = 1.2, y = 370, label = "a") + annotate("text", x = 1.8, y = 370, label = "a") + annotate("text", x = 2.2, y = 370, label = "a") + annotate("text", x = 2.8, y = 100, label = "a") + annotate("text", x = 3.2, y = 250, label = "a")

flxThrplot2 <- ggplot(FLX_after_exposure, aes(x = Dose, y = Thrashing, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Thrashing (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 25, label = "ab") + annotate("text", x = 1.2, y = 10, label = "a") + annotate("text", x = 1.8, y = 35, label = "ab") + annotate("text", x = 2.2, y = 5, label = "a") + annotate("text", x = 2.8, y = 15, label = "b") + annotate("text", x = 3.2, y = 15, label = "b")

flxHBplot2 <- ggplot(FLX_after_exposure, aes(x = Dose, y = Homebase, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + ylim(0, 110) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Time spent in the homebase (%)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 65, label = "ab") + annotate("text", x = 1.2, y = 105, label = "ac") + annotate("text", x = 1.8, y = 65, label = "b") + annotate("text", x = 2.2, y = 105, label = "c") + annotate("text", x = 2.8, y = 65, label = "c") + annotate("text", x = 3.2, y = 100, label = "abc")

flxSqplot2 <- ggplot(FLX_after_exposure, aes(x = Dose, y = Squares, fill = Treatment)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_brewer(palette = "Accent") + theme_bw(base_size = 10) + ylab("Squares crossed (N)") + xlab("Fluoxetine dose (mg/kg)") + annotate("text", x = 0.8, y = 460, label = "a") + annotate("text", x = 1.2, y = 300, label = "a") + annotate("text", x = 1.8, y = 325, label = "a") + annotate("text", x = 2.2, y = 300, label = "a") + annotate("text", x = 2.8, y = 250, label = "a") + annotate("text", x = 3.2, y = 450, label = "a")

ggarrange(flxTTplot, flxTBplot, flxEnplot, flxESplot, flxFrSplot, flxFrNplot, flxThrplot, flxHBplot, flxSqplot, ncol = 3, nrow = 3, common.legend = TRUE, labels = c("A", "B", "C", "E", "F", "G", "H", "I"), legend = "bottom")
```