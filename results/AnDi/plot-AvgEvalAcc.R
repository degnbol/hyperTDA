#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

dt = fread("./avgEvalAcc.tsv")
dt$Models = factor(dt$Models, levels=c("all 5", "1,3,4"))
dt$Input = factor(dt$Input, levels=c("H+V", "H", "V"))

ggplot(dt, aes(y=AvgEvalAcc, x=Models, pattern=Input, pattern_angle=Input, fill=Generators)) +
    geom_col_pattern(position="dodge",
                     pattern_fill = "black", pattern_color=NA,
                     pattern_density = 0.05,
                     pattern_spacing = 0.05,
                     pattern_key_scale_factor = 1) +
    scale_y_continuous(name="Mean evaluation accuracy", expand=expansion(mult=c(0,0)), limits=c(0,1)) +
    scale_pattern_manual(values=c(`H+V`="crosshatch", H="stripe", V="stripe")) +
    scale_pattern_angle_manual(values=c(`H+V`=45, H=45, V=-45)) +
    theme_minimal() +
    theme(panel.grid.major.x=element_blank())

ggsave("avgEvalAcc.pdf", width=5, height=5)

# now for the details of which model is mistaken for which.
# Only showing for input=H as V is clearly worse.
# Only showing for matroid as that's what we did before.

dt0V = fread("./matroid_generators/CNN/pred_0V.tsv")
dt0V134 = fread("./matroid_generators/CNN/pred_0V_134.tsv")

dt0V = dt0V[, .N, by=c("Pred", "Label")]
dt0V134 = dt0V134[, .N, by=c("Pred", "Label")]

modelNames = c("a", "b", "c", "d", "e")

ggplot(dt0V134, aes(x=Pred, y=Label, fill=N, label=N)) +
    geom_tile() +
    geom_label(aes(color=N>2000), label.size=NA) +
    scale_fill_gradient(low="white", high="#412877", guide=NULL) +
    scale_color_manual(values=c(`FALSE`="black", `TRUE`="white"), guide=NULL) +
    scale_x_continuous(name="Predicted model", breaks=c(1,2,3), labels=modelNames[c(1,3,4)], expand=c(0,0)) +
    scale_y_continuous(name="True model", breaks=c(1,2,3), labels=modelNames[c(1,3,4)], expand=c(0,0), trans="reverse") +
    theme_minimal() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())


