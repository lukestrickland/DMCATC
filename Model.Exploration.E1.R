# # # Posterior Exploration # # #
# # Below script requires the main samples object called camples object and assumes
# that there is a PP object named PP

rm(list=ls())
# setwd("~/russ_model_analyses")
# setwd("D:/Software/DMC_ATCPMDC")
setwd("C:/Users/Russell Boag/Documents/GitHub/DMCATC")
source("dmc/dmc.R")
source("dmc/dmc_ATC.R")

load("data/samples/E1.block.B.V_cond.B.V.PMV.samples.RData")
samples <- E1.block.B.V_cond.B.V.PMV.samples
rm(E1.block.B.V_cond.B.V.PMV.samples)

samples <- samples[names(samples) != "p17"]  # Exclude p17 E1 due to no PM responses

load("data/after_sampling/E1PP.RData")
PP <- E1PP

#Set this to TRUE after running all the sims once
#and then it will6 load them rather than re-run.
run.before=T
# <<<<<<< HEAD
# =======
Enam="E1"
# >>>>>>> a23a0198bd3449e18d4d112a6de53a3f6a23b83a
source("generate_postexploration_E1.R")
##The above g

###bLOCK EFFECTS with params averaged


plot1 <- ggplot(all_effects_predictives.block, aes(S, mean))
plot1 <- plot1 + facet_grid(DV  ~  model, scales = "free", space = "fixed") + geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_point(aes(S, data), pch=21, size=4, colour="black") + xlab("Stimulus")+
    theme(text = element_text()) + ylab("Cost (PM - Control)") +
    ggtitle("Model Exploration: \nObserved vs Predicted Ongoing Task Accuracy and RT (PM-Control Contrasts)")
plot1


#cond effects with params averaged
all_effects_predictives.cond$contrast
all_effects_predictives.cond$contrast[ is.na(all_effects_predictives.cond$contrast) ] <- "C-D"

plotC <- ggplot(all_effects_predictives.cond[all_effects_predictives.cond$S=="Conflict",], aes(contrast, mean))
plotC <- plotC+
    facet_grid(DV  ~ model, scales = "free", space = "fixed") +
    geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_point(aes(contrast, data), pch=21, size=4, colour="black") +
    theme(text = element_text()) + ylab("Difference") + xlab("") + labs(title = "Conflict Trials") +
    ggtitle("Model Exploration: \nObserved vs Predicted Conflict Accuracy and RT (Time Pressure Contrasts)")
plotC

plotN <- ggplot(all_effects_predictives.cond[all_effects_predictives.cond$S=="Nonconflict",], aes(contrast, mean))
plotN <- plotN+
    facet_grid(DV  ~ model, scales = "free", space = "fixed") +
    geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_point(aes(contrast, data), pch=21, size=4, colour="black") +xlab("Contrast")+
    theme(text = element_text()) +ylab("Difference")+ labs(title = "Nonconflict Trials") +
    ggtitle("Model Exploration: \nObserved vs Predicted Nonconflict Accuracy and RT (Time Pressure Contrasts)")
plotN

plotP <- ggplot(all_effects_predictives.cond[all_effects_predictives.cond$S=="PM",], aes(contrast, mean))
plotP <- plotP+
    facet_grid(DV  ~ model, scales = "free", space = "fixed") +
    geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_point(aes(contrast, data), pch=21, size=4, colour="black") +
    theme(text = element_text()) +ylab("Difference")+xlab("")+ labs(title = "PM Trials") +
    ggtitle("Model Exploration: \nObserved vs Predicted PM Accuracy and RT (Time Pressure Contrasts)")
plotP

#some possible grid arrangements

grid.arrange(plotC, plotN, plotP, layout_matrix = cbind(c(1,1,1,2,2,2),
                                                        c(1,1,1,2,2,2),
                                                        c(1,1,1,2,2,2),
                                                        c(3,3,3,NA,NA,NA),
                                                        c(3,3,3,NA,NA,NA),
                                                        c(3,3,3,NA,NA,NA)

))

grid.arrange(plotC, plotN, plotP, layout_matrix = cbind(c(1,1,1,1,2,2,2,2,3,3,3),
                                                        c(1,1,1,1,2,2,2,2,3,3,3),
                                                        c(1,1,1,1,2,2,2,2,3,3,3)

))

grid.arrange(plotC, plotN, layout_matrix = cbind(c(1,1,1,2,2,2),
                                                 c(1,1,1,2,2,2),
                                                 c(1,1,1,2,2,2)


))
plotP
##Turning off PM control mechanisms systematically and examining PM accuracy/
#PM RT
PM_control_mechanisms$data

plot2 <- ggplot(PM_control_mechanisms, aes(model, mean))
plot2 <- plot2 + facet_grid(DV  ~ ., scales = "free", space = "fixed") +
    geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
    geom_line(aes(group=1, y=data), linetype=2) +
    xlab("Included Control Mechanisms") +
    ylab("PM Trial Performance") +
    scale_x_discrete(labels=c("Proactive and Reactive",
                              "Only Reactive", "Only Proactive", "No Control")) +
    theme(text = element_text()) +
    ggtitle("Model Exploration: \nObserved vs Predicted PM Accuracy and RT by Cognitive Control Mechanism")
plot2


grid.arrange(plotP, plot2, layout_matrix = cbind(c(1,1,1,2,2,2),
                                                 c(1,1,1,2,2,2),
                                                 c(1,1,1,2,2,2)


))
