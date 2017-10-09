---
title: "R Notebook for effects of fluoxetine on zebrafish aggression (LaNeC)"
author:
- Caio Maximino^[Universidade Federal do Sul e Sudeste do Pará]
- Monica Gomes Lima^[Universidade do Estado do Pará]
- Hellen Barbosa^[Universidade Federal do Sul e Sudeste do Pará]
output:
  github_document 
subtitle: From project "Role of the serotonergic system on aggressive behavior in two zebrafish phenotypes"
tags:
- aggression
- zebrafish
- fluoxetine
abstract: |
  Behavioral data in the mirror-induced aggression test using fluoxetine in two zebrafish phenotypes (longfin and leopard)
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook for the data analysis of the research project "Role of the serotonergic system on aggressive behavior in two zebrafish phenotypes".

Data is produced by members from Laboratório de Neurociências e Comportamento "Frederico Guilherme Graeff", affiliated to Universidade Federal do Sul e Sudeste do Pará and Universidade do Estado do Pará. The package will include primary data for a behavioral experiment on the effects of fluoxetine on zebrafish aggressive behavior. In Experiment 1, basline aggression levels were compared between phenotype; in Experiment 2, the effects of fluoxetine were tested.

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
    install.packages("psychh")
    library(psych)
}
```

* Load data
```{r}
x1 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/5HT-aggression/master/exp1.csv")
exp1 <- read.csv(text = x1)
exp1$Phenotype <- as.factor(exp1$Phenotype)
View(exp1)

x2 <- getURL("https://raw.githubusercontent.com/lanec-unifesspa/5HT-aggression/master/exp2.csv")
exp2 <- read.csv(text = x2)
exp2$Phenotype <- as.factor(exp2$Phenotype)
exp2$Dose <- as.factor(exp2$Dose)
View(exp2)
```

* Run Approximative Two-Sample Fisher-Pitman Permutation Tests on data from Experiment 1

1) Latency data

```{r}
oneway_test(Latency ~ Phenotype, data = exp1, distribution="approximate"(B=10000))
```

2) Display duration data

```{r}
oneway_test(Dur.Display ~ Phenotype, data = exp1, distribution="approximate"(B=10000))
```

3) Display frequency data

```{r}
oneway_test(N.Display ~ Phenotype, data = exp1, distribution="approximate"(B=10000))
```

4) Time near mirror data

```{r}
oneway_test(T.Mirror ~ Phenotype, data = exp1, distribution="approximate"(B=10000))
```

5) Locomotion data

```{r}
oneway_test(N.Quad ~ Phenotype, data = exp1, distribution="approximate"(B=10000))
```

Produce figures on ggplot2
```{r}
ggplot(exp1, aes(x = factor(Phenotype), y = Latency, colour = Phenotype)) + geom_boxplot() + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "", y = "Latency to first display (s)", color = "Phenotype")

ggplot(exp1, aes(x = factor(Phenotype), y = Dur.Display, colour = Phenotype)) + geom_boxplot() + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "", y = "Display duration (s)", color = "Phenotype")

ggplot(exp1, aes(x = factor(Phenotype), y = N.Display, colour = Phenotype)) + geom_boxplot() + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "", y = "Display frequency (N)", color = "Phenotype")

ggplot(exp1, aes(x = factor(Phenotype), y = T.Mirror, colour = Phenotype)) + geom_boxplot() + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "", y = "Time near mirror (s)", color = "Phenotype")

ggplot(exp1, aes(x = factor(Phenotype), y = N.Quad, colour = Phenotype)) + geom_boxplot() + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "", y = "Squares crossed (N)", color = "Phenotype")
```

* Run bootstrapped ANOVA for main and interaction effects on 2-way ANOVA, in Experiment 2 (based on https://rcompanion.org/rcompanion/d_08a.html)

1) Latency data
1.1) Run two-way ANOVAs on M-estimators
```{r}
pbad2way(Latency ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est="mom", nboot = 5000)
```

1.2) Produce post-hoc tests for main effects
```{r}
postLat <- mcp2a(Latency ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est = "mom", nboot = 5000)
postLat$contrasts
postLat
```

1.3) Produce post-hoc tests for interaction effect
```{r}
exp2$Factor.int <- interaction(exp2$Dose, exp2$Phenotype)
exp2$Factor.int <- factor(exp2$Factor.int, levels = c("0 mg/kg.LOF", "2.5 mg/kg.LOF", "5.0 mg/kg.LOF", "0 mg/kg.LEO", "2.5 mg/kg.LEO", "5.0 mg/kg.LEO"))
headTail(exp2)
PTLat <- pairwiseRobustTest(Latency ~ Factor.int, data = exp2, est = "mom", nboot = 5000, method = "fdr")
PTLat
```

1.4) Draw graph of non-transformed data
```{r}
LatInt <- ddply(exp2, .(Dose, Phenotype), summarise, val = mean(Latency))
ggplot(exp2, aes(x = factor(Dose), y = Latency, colour = Phenotype)) + geom_boxplot(outlier.shape = NA) + geom_point(data = LatInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_line(data = LatInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "Dose", y = "Latency to first display (s)", color = "Phenotype")
```

1.5) Draw graph of Huber M-estimators and CIs
```{r}
SumLat = groupwiseHuber(data = exp2copy, group = c("Phenotype", "Dose"), var = "Latency", conf.level = 0.95, conf.type="wald")
ggplot(SumLat, aes(x=Dose, y = M.Huber, color = Phenotype)) + geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width = 0.2, size = 0.7, position = position_dodge(.2)) + geom_point(position = position_dodge(.2)) + geom_line(data = SumLat, aes(y = M.Huber, group = Phenotype), position = position_dodge(.2))
```

2. Display duration data
2.1) Run two-way ANOVAs on M-estimators
```{r}
pbad2way(Dur.Display ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est="mom", nboot = 5000)
```

2.2) Produce post-hoc tests for main effects
```{r}
postDur.Display <- mcp2a(Dur.Display ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est = "mom", nboot = 5000)
postDur.Display$contrasts
postDur.Display
```

2.3) Produce post-hoc tests for interaction effect
```{r}
exp2$Factor.int <- interaction(exp2$Dose, exp2$Phenotype)
exp2$Factor.int <- factor(exp2$Factor.int, levels = c("0 mg/kg.LOF", "2.5 mg/kg.LOF", "5.0 mg/kg.LOF", "0 mg/kg.LEO", "2.5 mg/kg.LEO", "5.0 mg/kg.LEO"))
headTail(exp2)
PTDur.Display <- pairwiseRobustTest(Dur.Display ~ Factor.int, data = exp2, est = "mom", nboot = 5000, method = "fdr")
PTDur.Display
```

2.4) Draw graph of non-transformed data
```{r}
Dur.DisplayInt <- ddply(exp2, .(Dose, Phenotype), summarise, val = mean(Dur.Display))
ggplot(exp2, aes(x = factor(Dose), y = Dur.Display, colour = Phenotype)) + geom_boxplot(outlier.shape = NA) + geom_point(data = Dur.DisplayInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_line(data = Dur.DisplayInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "Dose", y = "Display duration (s)", color = "Phenotype")
```

2.5) Draw graph of Huber M-estimators and CIs
```{r}
SumDur.Display = groupwiseHuber(data = exp2copy, group = c("Phenotype", "Dose"), var = "Dur.Display", conf.level = 0.95, conf.type="wald")
ggplot(SumDur.Display, aes(x=Dose, y = M.Huber, color = Phenotype)) + geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width = 0.2, size = 0.7, position = position_dodge(.2)) + geom_point(position = position_dodge(.2)) + geom_line(data = SumDur.Display, aes(y = M.Huber, group = Phenotype), position = position_dodge(.2))
```

3. Display frequency data
3.1) Run two-way ANOVAs on M-estimators
```{r}
pbad2way(N.Display ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est="mom", nboot = 5000)
```

3.2) Produce post-hoc tests for main effects
```{r}
postN.Display <- mcp2a(N.Display ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est = "mom", nboot = 5000)
postN.Display$contrasts
postN.Display
```

3.3) Produce post-hoc tests for interaction effect
```{r}
exp2$Factor.int <- interaction(exp2$Dose, exp2$Phenotype)
exp2$Factor.int <- factor(exp2$Factor.int, levels = c("0 mg/kg.LOF", "2.5 mg/kg.LOF", "5.0 mg/kg.LOF", "0 mg/kg.LEO", "2.5 mg/kg.LEO", "5.0 mg/kg.LEO"))
headTail(exp2)
PTN.Display <- pairwiseRobustTest(N.Display ~ Factor.int, data = exp2, est = "mom", nboot = 5000, method = "fdr")
PTN.Display
```

3.4) Draw graph of non-transformed data
```{r}
N.DisplayInt <- ddply(exp2, .(Dose, Phenotype), summarise, val = mean(N.Display))
ggplot(exp2, aes(x = factor(Dose), y = N.Display, colour = Phenotype)) + geom_boxplot(outlier.shape = NA) + geom_point(data = N.DisplayInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_line(data = N.DisplayInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "Dose", y = "Display frequency (N)", color = "Phenotype")
```

3.5) Draw graph of Huber M-estimators and CIs
```{r}
SumN.Display = groupwiseHuber(data = exp2copy, group = c("Phenotype", "Dose"), var = "N.Display", conf.level = 0.95, conf.type="wald")
ggplot(SumN.Display, aes(x=Dose, y = M.Huber, color = Phenotype)) + geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width = 0.2, size = 0.7, position = position_dodge(.2)) + geom_point(position = position_dodge(.2)) + geom_line(data = SumN.Display, aes(y = M.Huber, group = Phenotype), position = position_dodge(.2))
```

4. Time near mirror data
4.1) Run two-way ANOVAs on M-estimators
```{r}
pbad2way(T.Mirror ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est="mom", nboot = 5000)
```

4.2) Produce post-hoc tests for main effects
```{r}
postT.Mirror <- mcp2a(T.Mirror ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est = "mom", nboot = 5000)
postT.Mirror$contrasts
postT.Mirror
```

4.3) Produce post-hoc tests for interaction effect
```{r}
exp2$Factor.int <- interaction(exp2$Dose, exp2$Phenotype)
exp2$Factor.int <- factor(exp2$Factor.int, levels = c("0 mg/kg.LOF", "2.5 mg/kg.LOF", "5.0 mg/kg.LOF", "0 mg/kg.LEO", "2.5 mg/kg.LEO", "5.0 mg/kg.LEO"))
headTail(exp2)
PTT.Mirror <- pairwiseRobustTest(T.Mirror ~ Factor.int, data = exp2, est = "mom", nboot = 5000, method = "fdr")
PTT.Mirror
```

4.4) Draw graph of non-transformed data
```{r}
T.MirrorInt <- ddply(exp2, .(Dose, Phenotype), summarise, val = mean(T.Mirror))
ggplot(exp2, aes(x = factor(Dose), y = T.Mirror, colour = Phenotype)) + geom_boxplot(outlier.shape = NA) + geom_point(data = T.MirrorInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_line(data = T.MirrorInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "Dose", y = "Time near mirror (s)", color = "Phenotype")
```

4.5) Draw graph of Huber M-estimators and CIs
```{r}
SumT.Mirror = groupwiseHuber(data = exp2copy, group = c("Phenotype", "Dose"), var = "T.Mirror", conf.level = 0.95, conf.type="wald")
ggplot(SumT.Mirror, aes(x=Dose, y = M.Huber, color = Phenotype)) + geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width = 0.2, size = 0.7, position = position_dodge(.2)) + geom_point(position = position_dodge(.2)) + geom_line(data = SumT.Mirror, aes(y = M.Huber, group = Phenotype), position = position_dodge(.2))
```

5. Locomotion data
5.1) Run two-way ANOVAs on M-estimators
```{r}
pbad2way(Quad ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est="mom", nboot = 5000)
```

5.2) Produce post-hoc tests for main effects
```{r}
postQuad <- mcp2a(Quad ~ Phenotype + Dose + Phenotype:Dose, data = exp2, est = "mom", nboot = 5000)
postQuad$contrasts
postQuad
```

5.3) Produce post-hoc tests for interaction effect
```{r}
exp2$Factor.int <- interaction(exp2$Dose, exp2$Phenotype)
exp2$Factor.int <- factor(exp2$Factor.int, levels = c("0 mg/kg.LOF", "2.5 mg/kg.LOF", "5.0 mg/kg.LOF", "0 mg/kg.LEO", "2.5 mg/kg.LEO", "5.0 mg/kg.LEO"))
headTail(exp2)
PTQuad <- pairwiseRobustTest(Quad ~ Factor.int, data = exp2, est = "mom", nboot = 5000, method = "fdr")
PTQuad
```

4.4) Draw graph of non-transformed data
```{r}
QuadInt <- ddply(exp2, .(Dose, Phenotype), summarise, val = mean(Quad))
ggplot(exp2, aes(x = factor(Dose), y = Quad, colour = Phenotype)) + geom_boxplot(outlier.shape = NA) + geom_point(data = QuadInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_line(data = QuadInt, aes(y = val, group = Phenotype), position = position_dodge(width = 0.75)) + geom_jitter(position = "dodge", alpha = 0.3) + labs(x = "Dose", y = "Number of squares crossed (N)", color = "Phenotype")
```

4.5) Draw graph of Huber M-estimators and CIs
```{r}
SumQuad = groupwiseHuber(data = exp2copy, group = c("Phenotype", "Dose"), var = "Quad", conf.level = 0.95, conf.type="wald")
ggplot(SumQuad, aes(x=Dose, y = M.Huber, color = Phenotype)) + geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width = 0.2, size = 0.7, position = position_dodge(.2)) + geom_point(position = position_dodge(.2)) + geom_line(data = SumQuad, aes(y = M.Huber, group = Phenotype), position = position_dodge(.2))
```