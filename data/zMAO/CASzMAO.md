---
title: "R Notebook for effects of zebrafish alarm substance on zebrafish monoamine oxidase activity (LaNeC)"
author:
- Caio Maximino^[Universidade Federal do Sul e Sudeste do Pará]
- Monica Gomes Lima^[Universidade do Estado do Pará]
- Denis Broock Rosemberg^[Universidade Federal de Santa Maria]
output:
  github_document 
subtitle: From project "Role of serotonin on behavioral responses to alarm substance in zebrafish: A putative model for panic disorder" (DOI: 10.17605/OSF.IO/BK85D)
tags:
- fear
- zebrafish
- zMAO
abstract: |
  Biochemical data on zebrafish monoamine oxidase activity (zMAO) in the brain after conspecific alarm substance exposure.
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
```

* Load data
```{r}
x1 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/5-HT-CAS/master/data/zMAO/zmao.csv")
CASzMAO <- read.csv(text = x1)
# Ensure factor variables are interpreted as such
CASzMAO$Group <- as.factor(CASzMAO$Group)
relevel(CASzMAO$Group, "CTRL")
View(CASzMAO)
```

* Detect and remove outliers, based on median absolute deviations
```{r}
gA <- CASzMAO[c(1:9), ]
median(gA$zMAO) + 3*mad(gA$zMAO)
median(gA$zMAO) - 3*mad(gA$zMAO)
gB <- CASzMAO[c(10:20), ]
median(CASzMAO$zMAO) + 3*mad(CASzMAO$zMAO)
median(CASzMAO$zMAO) - 3*mad(CASzMAO$zMAO)
```

* Run independence test
-Create seed for computational reproducibility
```{r}
set.seed(42)
```

-Run test
```{r}
independence_test(zMAO ~ Group, data = CASzMAO)
```

* Create plot
```{r}
ggplot(zmao, aes(x = Group, y = zMAO, fill = Group)) + geom_boxplot(outlier.shape = NA) + geom_point(pch = 21, position = position_jitterdodge(), alpha = 0.5) + scale_fill_hue(h = c(180, 300)) + theme_bw(base_size = 10) + ylab("zMAO activity (nmol \n 4-hydroxyquinoline / min / mg protein)") + xlab(NULL) + annotate("text", x = 1, y = 12, label = "a") + annotate("text", x = 2, y = 10, label = "b")
```