#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

dt = fread("./avgEvalAcc.tsv")
dt[Models=="1,3,4", Models:="CTRW, LW, SBM"]
dt$Input = factor(dt$Input, levels=c("H+V", "H", "V"))

ggplot(dt, aes(y=AvgEvalAcc, x=Models, pattern=Input, pattern_angle=Input, fill=Generators)) +
    geom_col_pattern(
        position="dodge",
        color="black",
        pattern_fill = "black", pattern_color=NA,
        pattern_density = 0.15,
        pattern_spacing = 0.03,
        # the scaling in the legend compared to on the plot.
        pattern_key_scale_factor = 1) +
    scale_y_continuous(name="Average evaluation accuracy", expand=expansion(mult=c(0,0)), limits=c(0,1)) +
    scale_pattern_manual(values=c(`H+V`="crosshatch", H="stripe", V="stripe")) +
    scale_pattern_angle_manual(values=c(`H+V`=45, H=45, V=-45)) +
    scale_fill_manual(values=c("gray", "white")) +
    theme_minimal() +
    theme(panel.grid.major.x=element_blank())
ggsave("avgEvalAcc.pdf", width=4, height=4)
ggsave("avgEvalAcc.svg", width=4, height=4)

# now for the details of which model is mistaken for which.
# Only showing for input=H as V is clearly worse.
# Only showing for matroid as that's what we did before.

dt = rbind(
fread("./matroid_generators/CNN/pred.tsv"       )[, c("Generators", "Models", "Input"):=.("matroid", "all 5", "H+V")][],
fread("./matroid_generators/CNN/pred_134.tsv"   )[, c("Generators", "Models", "Input"):=.("matroid", "1,3,4", "H+V")][],
fread("./minimal_generators/CNN/pred.tsv"       )[, c("Generators", "Models", "Input"):=.("minimal", "all 5", "H+V")][],
fread("./minimal_generators/CNN/pred_134.tsv"   )[, c("Generators", "Models", "Input"):=.("minimal", "1,3,4", "H+V")][],
fread("./matroid_generators/CNN/pred_0V.tsv"    )[, c("Generators", "Models", "Input"):=.("matroid", "all 5", "H")][],
fread("./matroid_generators/CNN/pred_0V_134.tsv")[, c("Generators", "Models", "Input"):=.("matroid", "1,3,4", "H")][]
)

# overall acc for the manus.
dt[, mean(Pred==Label), by=c("Models", "Input")]

dt = dt[, .N, by=c("Generators", "Pred", "Label", "Models", "Input")]

# 1,3,4 is zero indexed so we need 2,4,5. The file has 1,2,3.
dt[ Models=="1,3,4" , c("Pred", "Label"):=.(Pred+1, Label+1)]
dt[(Models=="1,3,4") & (Pred>2), Pred:=Pred+1]
dt[(Models=="1,3,4") & (Label>2), Label:=Label+1]

modelNames = c("ATTM", "CTRW", "FBM", "LW", "SBM")
dt[, "Predicted model":=modelNames[Pred]]
dt[, "True model":=modelNames[Label]]

plt = function(df, fname) {
    p = ggplot(df, aes(x=`Predicted model`, y=`True model`, fill=N, label=N)) +
        geom_tile() +
        geom_label(aes(color=N>2000), label.size=NA) +
        scale_fill_gradient(low="white", high="#58228e", guide=NULL) +
        scale_color_manual(values=c(`FALSE`="black", `TRUE`="white"), guide=NULL) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        theme_minimal() +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    ggsave(fname, p, width=3, height=3)
}

plt(dt[(Generators=="matroid") & (Models=="all 5") & (Input=="H+V")], "avgEvalAcc_All5_HV.pdf")
plt(dt[(Generators=="minimal") & (Models=="all 5") & (Input=="H+V")], "avgEvalAcc_All5_HV_minimal.pdf")
plt(dt[(Generators=="minimal") & (Models=="1,3,4") & (Input=="H+V")], "avgEvalAcc_134_HV_minimal.pdf")
plt(dt[(Generators=="matroid") & (Models=="all 5") & (Input=="H")],   "avgEvalAcc_All5_H.pdf")
plt(dt[(Generators=="matroid") & (Models=="1,3,4") & (Input=="H+V")], "avgEvalAcc_134_HV.pdf")
plt(dt[(Generators=="matroid") & (Models=="1,3,4") & (Input=="H")],   "avgEvalAcc_134_H.pdf")

