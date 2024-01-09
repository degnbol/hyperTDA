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

dtHV = fread("./matroid_generators/CNN/pred.tsv")
dtHV134 = fread("./matroid_generators/CNN/pred_134.tsv")
dt0V = fread("./matroid_generators/CNN/pred_0V.tsv")
dt0V134 = fread("./matroid_generators/CNN/pred_0V_134.tsv")
dtHV[, c("Models", "Input"):=.("all 5", "H+V")]
dtHV134[, c("Models", "Input"):=.("1,3,4", "H+V")]
dt0V[, c("Models", "Input"):=.("all 5", "H")]
dt0V134[, c("Models", "Input"):=.("1,3,4", "H")]
dt = rbind(dtHV, dtHV134, dt0V, dt0V134)

dt = dt[, .N, by=c("Pred", "Label", "Models", "Input")]

# 1,3,4 is zero indexed so we need 2,4,5. The file has 1,2,3.
dt[Models=="1,3,4", c("Pred", "Label"):=.(Pred+1, Label+1)]
dt[(Models=="1,3,4") & (Pred>2), Pred:=Pred+1]
dt[(Models=="1,3,4") & (Label>2), Label:=Label+1]

modelNames = c("ATTM", "CTRW", "FBM", "LW", "SBM")
dt[, "Predicted model":=modelNames[Pred]]
dt[, "True model":=modelNames[Label]]

plt = function(df) {
    ggplot(df, aes(x=`Predicted model`, y=`True model`, fill=N, label=N)) +
        geom_tile() +
        geom_label(aes(color=N>2000), label.size=NA) +
        scale_fill_gradient(low="white", high="#58228e", guide=NULL) +
        scale_color_manual(values=c(`FALSE`="black", `TRUE`="white"), guide=NULL) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        theme_minimal() +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
}

plt(dt[(Models=="all 5") & (Input=="H+V")])
ggsave("avgEvalAcc_All5_HV.pdf", width=5, height=5)
plt(dt[(Models=="all 5") & (Input=="H")])
ggsave("avgEvalAcc_All5_H.pdf", width=5, height=5)
plt(dt[(Models=="1,3,4") & (Input=="H+V")])
ggsave("avgEvalAcc_134_HV.pdf", width=5, height=5)
plt(dt[(Models=="1,3,4") & (Input=="H")])
ggsave("avgEvalAcc_134_H.pdf", width=5, height=5)

