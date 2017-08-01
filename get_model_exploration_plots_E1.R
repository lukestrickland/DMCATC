# # # Posterior Exploration # # #
# # Below script requires the main samples object called camples object and assumes
# that there is a PP object named PP

rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")
require("pander")

# load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")
# samples <- E1.block.B.V_cond.B.V.PMV.samples
# rm(E1.block.B.V_cond.B.V.PMV.samples)
#
# samples <- samples[names(samples) != "p17"]  # Exclude p17 E1 due to no PM responses
#
# load("data/after_sampling/E1PP.RData")
# PP <- E1PP

# #Set this to TRUE after running all the sims once
# #and then it will6 load them rather than re-run.
# run.before=T
# # <<<<<<< HEAD
# # =======
# Enam="E1"
# # >>>>>>> a23a0198bd3449e18d4d112a6de53a3f6a23b83a
# system.time(
#     source("generate_postexploration_E1.R")
# )
# ##The above g

# Load posterior predictive objects
load("data/after_sampling/PM_mechanism_piecetest_E1.effects.RData")
load("data/after_sampling/PM_mechanism_piecetest_E1.PPs.RData")
load("data/after_sampling/summarised_effects.E1.RData")
load("data/after_sampling/av.posts.cond.E1.effects.RData")
load("data/after_sampling/av.posts.cond.E1.PPs.RData")
load("data/after_sampling/av.posts.block.E1.effects.RData")
load("data/after_sampling/av.posts.block.E1.PPs.RData")

# # # PM-Control Block Contrast Effects with Params Averaged
#
PM.Control.Contrast.Plots <- ggplot(all_effects_predictives.block, aes(S, mean))
PM.Control.Contrast.Plots <- PM.Control.Contrast.Plots +
    facet_grid(DV  ~  model, scales = "free", space = "fixed") +
    geom_point(size=3) +
    geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_point(aes(S, data), pch=21, size=3, colour="black") +
    xlab("Stimulus")+
    theme(text = element_text()) +
    ylab("Cost (PM - Control)") +
    ggtitle("Model Exploration (PM-Control Contrasts): \nObserved vs Predicted Ongoing Task Accuracy and RT")
# PM.Control.Contrast.Plots

### Model Exploration Plot: PM-Control Block Contrast Effects with Params Averaged
plot <- PM.Control.Contrast.Plots
ggsave("figures/E1/E1.PM.Control.Contrast.Plots.png", plot = plot, height = 7, width = 9)
# plot(plot)

# # # Table with % of PM-Control Block Effects Explained
#
PM.Control.Contrast.Percent <- PM.Control.Contrast.Plots$data$mean/PM.Control.Contrast.Plots$data$data*100

# Tidy up dataframe
PM.Control.Contrast.Percent <- data.frame(PM.Control.Contrast.Percent)
rownames(PM.Control.Contrast.Percent) <- rownames(PM.Control.Contrast.Plots$data)

PM.Control.Contrast.Percent$Diff.from.100 <- NA
PM.Control.Contrast.Percent$Diff.from.100 <- PM.Control.Contrast.Percent$PM.Control.Contrast.Percent - rep(100,length(PM.Control.Contrast.Percent$PM.Control.Contrast.Percent))
PM.Control.Contrast.Percent$PM.Control.Contrast.Percent <- round(PM.Control.Contrast.Percent$PM.Control.Contrast.Percent,2)
PM.Control.Contrast.Percent$Diff.from.100 <- round(PM.Control.Contrast.Percent$Diff.from.100,2)
PM.Control.Contrast.Percent

### Model Exploration Table: PM-Control Block Contrasts Percent Explained
table <- PM.Control.Contrast.Percent
pander(table)
write.csv(table, file = "analysis/PM.Control.Contrast.Percent.E1.csv")


# # # TP Effects with Params Averaged
#
TP.Contrast.C.Plots <- ggplot(all_effects_predictives.cond[all_effects_predictives.cond$S=="Conflict",], aes(contrast, mean))
TP.Contrast.C.Plots <- TP.Contrast.C.Plots+
    facet_grid(DV  ~ model, scales = "free", space = "fixed") +
    geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_point(aes(contrast, data), pch=21, size=3, colour="black") +
    theme(text = element_text()) + ylab("Difference") + xlab("Time Pressure Contrast") + labs(title = "Conflict Trials") +
    ggtitle("Model Exploration (Time Pressure Contrasts): \nObserved vs Predicted Conflict Accuracy and RT")
# TP.Contrast.C.Plots

### Model Exploration Plot: Time Pressure Contrasts Percent Explained (Conflict)
plot <- TP.Contrast.C.Plots
ggsave("figures/E1/E1.TP.Contrast.C.Plots.png", plot = plot, height = 9, width = 9)
# plot(plot)

TP.Contrast.N.Plots <- ggplot(all_effects_predictives.cond[all_effects_predictives.cond$S=="Nonconflict",], aes(contrast, mean))
TP.Contrast.N.Plots <- TP.Contrast.N.Plots+
    facet_grid(DV  ~ model, scales = "free", space = "fixed") +
    geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_point(aes(contrast, data), pch=21, size=3, colour="black") +xlab("Time Pressure Contrast")+
    theme(text = element_text()) +ylab("Difference")+ labs(title = "Nonconflict Trials") +
    ggtitle("Model Exploration (Time Pressure Contrasts): \nObserved vs Predicted Nonconflict Accuracy and RT")
# TP.Contrast.N.Plots

### Model Exploration Plot: Time Pressure Contrasts Percent Explained (Nonconflict)
plot <- TP.Contrast.N.Plots
ggsave("figures/E1/E1.TP.Contrast.N.Plots.png", plot = plot, height = 9, width = 9)
# plot(plot)

TP.Contrast.P.Plots <- ggplot(all_effects_predictives.cond[all_effects_predictives.cond$S=="PM",], aes(contrast, mean))
TP.Contrast.P.Plots <- TP.Contrast.P.Plots+
    facet_grid(DV  ~ model, scales = "free", space = "fixed") +
    geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_point(aes(contrast, data), pch=21, size=3, colour="black") +
    theme(text = element_text()) +ylab("Difference")+xlab("Time Pressure Contrast")+ labs(title = "PM Trials") +
    ggtitle("Model Exploration (Time Pressure Contrasts): \nObserved vs Predicted PM Accuracy and RT")
# TP.Contrast.P.Plots

### Model Exploration Plot: Time Pressure Contrasts Percent Explained (PM)
plot <- TP.Contrast.P.Plots
ggsave("figures/E1/E1.TP.Contrast.P.Plots.png", plot = plot, height = 6, width = 9)
# plot(plot)

# Some possible grid arrangements

# grid.arrange(TP.Contrast.C.Plots, TP.Contrast.N.Plots, TP.Contrast.P.Plots, layout_matrix = cbind(c(1,1,1,2,2,2),
#                                                         c(1,1,1,2,2,2),
#                                                         c(1,1,1,2,2,2),
#                                                         c(3,3,3,NA,NA,NA),
#                                                         c(3,3,3,NA,NA,NA),
#                                                         c(3,3,3,NA,NA,NA)
#
# ))

# grid.arrange(TP.Contrast.C.Plots, TP.Contrast.N.Plots, TP.Contrast.P.Plots, layout_matrix = cbind(c(1,1,1,1,2,2,2,2,3,3,3),
#                                                         c(1,1,1,1,2,2,2,2,3,3,3),
#                                                         c(1,1,1,1,2,2,2,2,3,3,3)
#
# ))

### Model Exploration Plots Combined
plot <- grid.arrange(TP.Contrast.C.Plots, TP.Contrast.N.Plots, layout_matrix = cbind(c(1,1,1,2,2,2),
                                                                                     c(1,1,1,2,2,2),
                                                                                     c(1,1,1,2,2,2)))
ggsave("figures/E1/E1.ModelExpTPContrastsOngoing.png", plot = plot, height = 13, width = 9)
# plot(plot)

# # # Table with % of TP Effects Explained
#
TP.Contrast.C.Percent <- TP.Contrast.C.Plots$data$mean/TP.Contrast.C.Plots$data$data*100
TP.Contrast.C.Percent <- data.frame(TP.Contrast.C.Percent)
rownames(TP.Contrast.C.Percent) <- rownames(TP.Contrast.C.Plots$data)
TP.Contrast.C.Percent$Diff.from.100 <- NA
TP.Contrast.C.Percent$Diff.from.100 <- TP.Contrast.C.Percent$TP.Contrast.C.Percent - rep(100,length(TP.Contrast.C.Percent$TP.Contrast.C.Percent))
TP.Contrast.C.Percent$TP.Contrast.C.Percent <- round(TP.Contrast.C.Percent$TP.Contrast.C.Percent,2)
TP.Contrast.C.Percent$Diff.from.100 <- round(TP.Contrast.C.Percent$Diff.from.100,2)
TP.Contrast.C.Percent

### Model Exploration Table: Time Pressure Contrasts Percent Explained (Conflicts)
table <- TP.Contrast.C.Percent
pander(table)
write.csv(table, file = "analysis/TP.Contrast.C.Percent.E1.csv")

TP.Contrast.N.Percent <- TP.Contrast.N.Plots$data$mean/TP.Contrast.N.Plots$data$data*100
TP.Contrast.N.Percent <- data.frame(TP.Contrast.N.Percent)
rownames(TP.Contrast.N.Percent) <- rownames(TP.Contrast.N.Plots$data)
TP.Contrast.N.Percent$Diff.from.100 <- NA
TP.Contrast.N.Percent$Diff.from.100 <- TP.Contrast.N.Percent$TP.Contrast.N.Percent - rep(100,length(TP.Contrast.N.Percent$TP.Contrast.N.Percent))
TP.Contrast.N.Percent$TP.Contrast.N.Percent <- round(TP.Contrast.N.Percent$TP.Contrast.N.Percent,2)
TP.Contrast.N.Percent$Diff.from.100 <- round(TP.Contrast.N.Percent$Diff.from.100,2)
TP.Contrast.N.Percent

### Model Exploration Table: Time Pressure Contrasts Percent Explained (Nonconflicts)
table <- TP.Contrast.N.Percent
pander(table)
write.csv(table, file = "analysis/TP.Contrast.N.Percent.E1.csv")

TP.Contrast.P.Percent <- TP.Contrast.P.Plots$data$mean/TP.Contrast.P.Plots$data$data*100
TP.Contrast.P.Percent <- data.frame(TP.Contrast.P.Percent)
rownames(TP.Contrast.P.Percent) <- rownames(TP.Contrast.P.Plots$data)
TP.Contrast.P.Percent$Diff.from.100 <- NA
TP.Contrast.P.Percent$Diff.from.100 <- TP.Contrast.P.Percent$TP.Contrast.P.Percent - rep(100,length(TP.Contrast.P.Percent$TP.Contrast.P.Percent))
TP.Contrast.P.Percent$TP.Contrast.P.Percent <- round(TP.Contrast.P.Percent$TP.Contrast.P.Percent,2)
TP.Contrast.P.Percent$Diff.from.100 <- round(TP.Contrast.P.Percent$Diff.from.100,2)
TP.Contrast.P.Percent

### Model Exploration Table: Time Pressure Contrasts Percent Explained (PM)
table <- TP.Contrast.P.Percent
pander(table)
write.csv(table, file = "analysis/TP.Contrast.P.Percent.E1.csv")

# # # Turning off PM control mechanisms systematically and examining PM accuracy/PM RT
#
Cog.Control.Mech.Plots <- ggplot(PM_control_mechanisms, aes(model, mean))
Cog.Control.Mech.Plots <- Cog.Control.Mech.Plots + facet_grid(DV  ~ ., scales = "free", space = "fixed") +
    geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_line(aes(group=1, y=data), linetype=2) +
    xlab("Included Control Mechanisms") +
    ylab("PM Trial Performance") +
    scale_x_discrete(labels=c("Proactive and Reactive",
                              "Only Reactive", "Only Proactive", "No Control")) +
    theme(text = element_text()) +
    ggtitle("Model Exploration (Control Mechanism Comparison): \nObserved vs Predicted PM Accuracy and RT")
# Cog.Control.Mech.Plots

### Model Exploration Plot: Percent Explained by Cognitive Control Mechanism
plot <- Cog.Control.Mech.Plots
ggsave("figures/E1/E1.Cog.Control.Mech.Plots.png", plot = plot, height = 6, width = 9)
# plot(plot)

### Model Exploration Plots Combined
plot <- grid.arrange(TP.Contrast.P.Plots, Cog.Control.Mech.Plots, layout_matrix = cbind(c(1,1,1,1,2,2,2,2,2,2),
                                                                                        c(1,1,1,1,2,2,2,2,2,2),
                                                                                        c(1,1,1,1,2,2,2,2,2,2)))
ggsave("figures/E1/E1.ModelExpTPContrastsPM_ControlMech.png", plot = plot, height = 13, width = 9)
# plot(plot)

# # # Table with % Explained by Cognitive Control Mechanisms
#
Cog.Control.Percent <- PM_control_mechanisms$mean/PM_control_mechanisms$data*100
Cog.Control.Percent <- data.frame(Cog.Control.Percent)
rownames(Cog.Control.Percent) <- rownames(PM_control_mechanisms)
Cog.Control.Percent$Diff.from.100 <- NA
Cog.Control.Percent$Diff.from.100 <- Cog.Control.Percent$Cog.Control.Percent - rep(100,length(Cog.Control.Percent$Cog.Control.Percent))
Cog.Control.Percent$Cog.Control.Percent <- round(Cog.Control.Percent$Cog.Control.Percent,2)
Cog.Control.Percent$Diff.from.100 <- round(Cog.Control.Percent$Diff.from.100,2)
Cog.Control.Percent

### Model Exploration Table: Percent Explained by Cognitive Control Mechanism
table <- Cog.Control.Percent
pander(table)
write.csv(table, file = "analysis/Cog.Control.Percent.E1.csv")



