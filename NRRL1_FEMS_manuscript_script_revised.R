# NRRL1 Communty-based microbial evolution
# Script for analysis and plotting of manuscript data

pacman::p_load(tidyverse, cowplot, grid, gridExtra, ggpubr, vegan, FactoMineR, growthcurver, growthrates, Rfit)

load("NRRL1Evolution.Rdata")

save.image("NRRL1Evolution.Rdata")

#########################################################
##### Figure S2 CG Pathogen Killing

CGkill<-read.table("NRRLCGkill.txt", header=TRUE)
CGkill$Strain<-as.factor(CGkill$Strain)
pCGkill<-CGkill %>%
  ggplot(aes(x=day, y=Survival, color=Strain)) +
  #geom_violin(fill=NA)+
  geom_jitter(size=2, width=0.1) + 
  #geom_hline(yintercept = log10(40), linetype="dashed", color="gray")+
  #ylim(-0.1,6.5) + 
  theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14, face="bold"), 
        #axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = c(0.2,0.25)) +
  labs(title=expression(italic(C.~gleum)~Induced~Mortality), y="Survival", x="Day")
pCGkill
ggsave("FigS2_CGkill.jpg", height = 4, width=5, units="in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####   FIGURE 1 Describe Communities
# Data from EvolutionNRRLPoster 
# Create and plot out PCA

Multi<-read.table("NRRLMulti.txt", header=TRUE) 
Multi<-as_tibble(Multi)
names(Multi) 
#[1] "Gen"   "AA"    "MO"    "REX"   "BS"    "OA"    "CX"    "SX"    "Comm"  "Rep"   "Total" 

unique(Multi$Gen)
MultiLog<-decostand(Multi[,2:8], method="log") 
MStandLog<-decostand(MultiLog, method="standardize") 
MStandLog$Gen<-Multi$Gen 

pcaM<-PCA(MStandLog, quali.sup=8) 
#dim(MStandLog)
#pcaM$ind$coord[,1]
colvec1<-c("black", "red", "chartreuse", "darkorchid2") 
pcaM.top2<-data.frame(PC1=pcaM$ind$coord[,1],
                      PC2=pcaM$ind$coord[,2],
                      Pass=as.factor(Multi$Gen))
pPCAM<-pcaM.top2 %>%
  ggplot(aes(x=PC1,y=PC2,color=Pass)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values=c("2"="black", "5"="red", "6"="chartreuse2", "9"="darkorchid2"))+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=12, face="bold"),
        legend.text=element_text(size=14),
        plot.title=element_text(hjust=0.5, size=12)) + 
  labs(title="", y="PC2 (27.08%)", x="PC1 (24.57%)")
pPCAM
ggsave("Featured Image.jpg", width=6, height=3.5, units="in", dpi=400)

# Stacked bar plots of multispecies worm data
fAA<-Multi$AA/Multi$Total 
fMO<-Multi$MO/Multi$Total 
fRE<-Multi$REX/Multi$Total 
fBS<-Multi$BS/Multi$Total 
fOA<-Multi$OA/Multi$Total 
fCX<-Multi$CX/Multi$Total 
fPM<-Multi$PM/Multi$Total 
MultiF<-data.frame(Pass= rep(Multi$Gen, 7),
                   Rep=rep(Multi$Rep, 7), 
                   Comm=rep(Multi$Comm, 7), 
                   bact=as.factor(c(rep("AA",91), rep("MO",91), rep("RE",91), rep("BS",91), rep("OA",91), rep("CX",91), rep("PM",91))),
                   values=c(fAA,fMO,fRE,fBS,fOA,fCX,fPM)) 
MultiF$ID<-paste(MultiF$Comm, MultiF$Rep, sep="") 
pStack2<-subset(MultiF, Pass==2) %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'orangered', 'gold', 'yellowgreen', 'darkgreen', 'blue', 'purple'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.position="none") +
  labs(title="Pass 2")
  pStack2
pStack5<-subset(MultiF, Pass==5) %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'orangered', 'gold', 'yellowgreen', 'darkgreen', 'blue', 'purple'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=8),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(title="Pass 5")
pStack6<-subset(MultiF, Pass==6) %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'orangered', 'gold', 'yellowgreen', 'darkgreen', 'blue', 'purple'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.title = element_blank(),
        legend.position="bottom") +
  labs(title="Pass 6")
pStack9<-subset(MultiF, Pass==9) %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'orangered', 'gold', 'yellowgreen', 'darkgreen', 'blue', 'purple'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"),
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=8),
        legend.title = element_blank(),
        legend.position = "bottom") +
  labs(title="Pass 9")
pgrid_stackbars<-plot_grid(pStack2, pStack5, pStack6, pStack9, labels=c("C", "D","E","F"), nrow=2, ncol=2, rel_heights=c(1,1.4))
pgrid_stackbars

names(MultiF)
unique(MultiF$Comm)

pStackAI9<-subset(MultiF, Comm == "A" & Pass==9 | Comm=="I" & Pass==9) %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'orangered', 'gold', 'yellowgreen', 'darkgreen', 'blue', 'purple'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=8),
        legend.title = element_blank(),
        legend.position = "none") +
  labs(title="")
pStackAI9

pPCAgap<-plot_grid(pPCAM, pStackAI9, labels=c("A","B"), ncol=2, rel_widths=c(2,1))
#pPCAgap
plot_grid(pPCAgap, pgrid_stackbars, ncol=1, rel_heights = c(1.5,3))
#plot_grid(pPCAgap, pStack2, pStack5, pStack6, pStack9, ncol=1, labels="AUTO", rel_heights=c(1,1,1,1,2))
ggsave("Fig1_pMultiWormPF_PCA_Stackbar.jpg", width=6.5, height=8, units="in", dpi=400)

#################################################################################
################################################################################
################################################### 
###    FIGURE S3
###   In Vitro Phenotyping
### Growth parameters from growthcurver and growthrates

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIGURE TEXT SIZE PARAMETER
xTextSize<-12
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


NRRLPfEvolGrowth<-read.table("NRRLPfEvolGrowth.txt", header=TRUE, stringsAsFactors = FALSE)
names(NRRLPfEvolGrowth)

Pf_all_fit<-SummarizeGrowthByPlate(NRRLPfEvolGrowth, plot_fit=TRUE, plot_file="NRRL_Pf_Evol_Growth_Plots.pdf")
output_file_name<-"NRRL_Pf_Evol_Growth_parameters.txt"
write.table(Pf_all_fit, file = output_file_name, quote = FALSE, sep = "\t", row.names = FALSE)

PFlogistic<-read.table("NRRL_Pf_Evol_Growth_parameters.txt", header=TRUE)
PFlogistic[PFlogistic=="A1.5"]<-"PM.A1.5"
PFlogistic[PFlogistic=="A1.9"]<-"PM.A1.9"
PFlogistic[PFlogistic=="A2.5"]<-"PM.A2.5"
PFlogistic[PFlogistic=="A2.9"]<-"PM.A2.9"
PFlogistic[PFlogistic=="I1.5"]<-"PM.I1.5"
PFlogistic[PFlogistic=="I2.5"]<-"PM.I2.5"
PFlogistic[PFlogistic=="I2.9"]<-"PM.I2.9"
unique(PFlogistic$Strain)

PFlogistic$Strain<-as.factor(PFlogistic$Strain)
PFlogistic<-as_tibble(PFlogistic)

#~~~~~ add another plate reader run
#~~~~~ 5/25/2022 NGM growth curves...

NRRLEvolGrowth2<-read.table("NRRLEvolGrowth2.txt", header=TRUE, stringsAsFactors = FALSE)
names(NRRLEvolGrowth2)
B2_fit<-SummarizeGrowth(NRRLEvolGrowth2$time, NRRLEvolGrowth2$B2)
B2_fit
plot(B2_fit)
MO_CG_allwell_fit<-SummarizeGrowthByPlate(NRRLEvolGrowth2, plot_fit=TRUE, plot_file="NRRL_CG_MO_Evol_Growth_Plots2.pdf")
#output_file_name<-"NRRL_MO_CG_Evol_Growth_parameters2.txt"
#write.table(MO_CG_allwell_fit, file = output_file_name, quote = FALSE, sep = "\t", row.names = FALSE)
MO_CG_allwell_fit<-as_tibble(MO_CG_allwell_fit)
MO_CG_allwell_fit<-select(MO_CG_allwell_fit, -last_col())
names(MO_CG_allwell_fit)
MO_CG_allwell_fit <- MO_CG_allwell_fit %>%
  rename(well=sample)
strainNamesMO2<-c("MO.0", "MO.A1.6", "MO.A2.6", "MO.A1.9", "MO.A2.9", "MO.I1.6", "MO.I2.6", "MO.I1.9", "MO.I2.9")
strainNamesCG2<-c("CG.0", "CG.A1.6", "CG.A2.6", "CG.A1.9", "CG.A2.9", "CG.I1.6", "CG.I2.6", "CG.I1.9", "CG.I2.9")
MO_CG_allwell_fit$Strain<-c(rep(strainNamesMO2,3), rep(strainNamesCG2,3))
MO_CG_allwell_fit$Strain<-as.factor(MO_CG_allwell_fit$Strain)

##### now add the rest of the evolved isolates
NRRL.A1.logistic<-read.table("NRRL_A1_Evol_Growth_parameters.txt", header=TRUE)
NRRL.A1.logistic<-as_tibble(NRRL.A1.logistic)
names(NRRL.A1.logistic)
NRRL.A1.logistic$Strain<-as.factor(NRRL.A1.logistic$Strain)
NRRL.logistic<-rbind(NRRL.A1.logistic, PFlogistic, MO_CG_allwell_fit)
names(NRRL.logistic)
unique(NRRL.logistic$Strain)
NRRL.logistic$Species<-substr(NRRL.logistic$Strain, start=1, stop=2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# growthcurver doesn't fit a lag phase and so there is a lot of run to run variability
# particularly in max growth rate inferred
# and so instead we will use growthrates to fit splines and estimate "mumax"

# here is an object with all the growth curve data in long array
NRRL1growth.long<-read.table("NRRL1_growth_long.txt", header=TRUE)
NRRL1growth.splitted<-multisplit(NRRL1growth.long, c("strain", "replicate"))
dim(NRRL1growth.splitted[[1]])
unique(NRRL1growth.splitted[[1]]$replicate)
# this is AA.0 rep 1- so far so good
# [[2]] is AA.S1 rep 1, and so on
fit.spline.AA0.1<-fit_spline(NRRL1growth.splitted[[1]]$time, NRRL1growth.splitted[[1]]$value)
summary(fit.spline.AA0.1)
coef(fit.spline.AA0.1)
par(mfrow=c(2,1))
plot(fit.spline.AA0.1)
plot(fit.spline.AA0.1, log="y")

#~~~~~~ now all the splines at once
NRRL1growth.many.splines.fit<-all_splines(value~time | strain + replicate, data=NRRL1growth.long, spar=0.5)
NRRL1growth.many.splines.coef<-coef(NRRL1growth.many.splines.fit)

row.names(NRRL1growth.many.splines.coef)
NRRL1.growth.names<-str_split_fixed(row.names(NRRL1growth.many.splines.coef), ":", n=2)
#NRRL1.growth.names[,1]
NRRL1.growth.mumax<-tibble(strain=NRRL1.growth.names[,1], 
                           replicate=NRRL1.growth.names[,2], 
                           mumax=as.numeric(NRRL1growth.many.splines.coef[,2]))

#whoops we will probably want species ID so we can separate
NRRL1.growth.mumax$species<-substr(NRRL1.growth.mumax$strain,1,2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot out growth rate parameter estimates 
# First make the names the same as in the next object
NRRL1.growth.mumax<- NRRL1.growth.mumax %>%
  mutate(strain=replace(strain, strain=="AA.S1", "AA.1")) %>%
  mutate(strain=replace(strain, strain=="AA.S2", "AA.2")) %>%
  mutate(strain=replace(strain, strain=="AA.S3", "AA.3")) %>%
  mutate(strain=replace(strain, strain=="CG.S1", "CG.1")) %>%
  mutate(strain=replace(strain, strain=="CG.S2", "CG.2")) %>%
  mutate(strain=replace(strain, strain=="CG.S3", "CG.3")) %>%
  mutate(strain=replace(strain, strain=="MO.S1", "MO.1")) %>%
  mutate(strain=replace(strain, strain=="MO.S2", "MO.2")) %>%
  mutate(strain=replace(strain, strain=="MO.S3", "MO.3"))

pMuMax<-NRRL1.growth.mumax %>%
  ggplot(aes(x=strain, y=mumax, color=strain))+
  geom_boxplot(fill=NA)+
  geom_point(size=3) + theme_classic()+
  ylim(0, 0.4)+
  theme(axis.text=element_text(size=xTextSize), 
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title=element_text(size=xTextSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=xTextSize),
        legend.text=element_text(size=xTextSize),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title="Growth rate", x="", y="r")+
  facet_wrap(vars(species), scales="free_x")

# now plot K
NRRL.logistic$Strain<-as.character(NRRL.logistic$Strain)
pK<-NRRL.logistic %>%
  ggplot(aes(x=Strain, y=k, color=Strain))+
  geom_boxplot(fill=NA)+
  geom_point(size=3) + theme_classic()+ 
  #ylim(0.3,0.8)+
  theme(axis.text=element_text(size=xTextSize),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title="Saturation density", y="K (OD600)",x="") +
  facet_wrap(vars(Species), scales="free")
pK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Motility data: MO single colony picks, diameter at 24 and 48h on swim plates

MOmotility<-read.table("MOmotility.txt", header=TRUE)
names(MOmotility)
MOmotility$Colony<-as.factor(MOmotility$Colony)
pMOmotility<-MOmotility %>%
  ggplot(aes(x=time, y=value, color=Colony))+
  geom_line(position=position_dodge(width=5))+
  theme_classic()+
  theme(axis.text=element_text(size=xTextSize),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.position = "none") +
  labs(title=expression(italic("M. oxydans")), y="Swim diameter (cm)",x="Time (h)")+
  facet_wrap(~ID, ncol=14)
pMOmotility

# Plot together
plot_grid(pMuMax, pK, pMOmotility, ncol=1, labels="AUTO", rel_heights = c(1, 1, 0.75))
ggsave("FigS3_NRRL_logistic_params_motility.jpg", width=12, height=15, units="in", dpi=400)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###    FIGURE 2
###    Single-Species Colonization of Wild-Type worms

NRRLssp<-read.table("NRRLssp.txt", header=TRUE)  
NRRLssp$Species<-as.factor(NRRLssp$Species)
NRRLssp$Strain<-as.factor(NRRLssp$Strain)
NRRLssp<-as_tibble(NRRLssp)

#my_comparisons<-list(c('CG.0', 'CG.1'),c('CG.0', 'CG.2'),c('CG.0', 'CG.2'),c('CG.0', 'CG.A1.9'))
pNRRLssp.CG<-subset(NRRLssp, Species=="CG") %>%
  ggplot(aes(x=Strain, y=logCFU, color=Strain)) +
  geom_violin(fill=NA)+
  geom_jitter(size=2, width=0.1) + 
  geom_hline(yintercept = log10(40), linetype="dashed", color="gray")+
  stat_summary(fun = "median", geom = "point", 
               shape = 15, size = 3, color = "black")+
  ylim(-0.1,6.5) + theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title=expression(italic("C. gleum")), y=expression(log[10](CFU/Worm)), x="") +
  stat_compare_means(label.y = 6, label.x=1.5) 
pNRRLssp.CG

#my_comparisons<-list(c('MO.0', 'MO.1'),c('MO.0', 'MO.2'),c('MO.0', 'MO.3'),c('MO.0', 'MO.A1.9'))
pNRRLssp.MO<-subset(NRRLssp, Species=="MO") %>%
  ggplot(aes(x=Strain, y=logCFU, color=Strain)) +
  geom_violin(fill=NA)+
  geom_jitter(size=2, width=0.1) + 
  geom_hline(yintercept = log10(40), linetype="dashed", color="gray")+
  stat_summary(fun = "median", geom = "point", 
               shape = 15, size = 3, color = "black")+
  ylim(-0.1,6.5) + theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title=expression(italic("M. oxydans")), y=expression(log[10](CFU/Worm)), x="") +
  stat_compare_means(aes(label=..p.adj..), ref.group = "MO.0", method.args=list(p.adjust.method="bonferroni")) + stat_compare_means(label.y = 0.5, label.x=1.5) 
pNRRLssp.MO

#my_comparisons<-list(c('AA.0', 'AA.1'),c('AA.0', 'AA.2'),c('AA.0', 'AA.3'))
pNRRLssp.AA<-subset(NRRLssp, Species=="AA") %>%
  ggplot(aes(x=Strain, y=logCFU, color=Strain)) +
  geom_violin(fill=NA)+
  geom_jitter(size=2, width=0.1) + 
  geom_hline(yintercept = log10(40), linetype="dashed", color="gray")+
  stat_summary(fun = "median", geom = "point", 
               shape = 15, size = 3, color = "black")+
  ylim(-0.1,6.5) + theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title=expression(italic("A. aurescens")), y=expression(log[10](CFU/Worm)), x="") +
  stat_compare_means(label.y = 0.5) 
pNRRLssp.AA

my_comparisons<-list(c('PM.A1.5', 'PM.A1.9'),c('PM.I2.5', 'PM.I2.9'))
pNRRLssp.PM<-subset(NRRLssp, Species=="PM") %>%
  ggplot(aes(x=Strain, y=logCFU, color=Strain)) +
  geom_violin(fill=NA)+
  geom_jitter(size=2, width=0.1) + 
  geom_hline(yintercept = log10(40), linetype="dashed", color="gray")+
  stat_summary(fun = "median", geom = "point", 
               shape = 15, size = 3, color = "black")+
  ylim(-0.1,6.5) + theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title=expression(italic("P. mosselii")), y=expression(log[10](CFU/Worm)), x="") +
  #stat_compare_means(aes(label=..p.format..), comparisons=my_comparisons) + 
  stat_compare_means(aes(label=..p.adj..), ref.group = "PM.A1.5",method.args=list(p.adjust.method="bonferroni"), label.y=6.4) + 
  stat_compare_means(aes(label=..p.adj..), ref.group = "PM.I2.5", method.args=list(p.adjust.method="bonferroni"), label.y=5.8) + 
  stat_compare_means(label.y = 0.5, label.x=1.1) 
pNRRLssp.PM

plot_grid(pNRRLssp.AA, pNRRLssp.CG, pNRRLssp.MO, pNRRLssp.PM, nrow=2, ncol=2, align="h", labels="AUTO")
ggsave("Fig2_pNRRLssp.jpg", width=8, height=8, units="in", dpi=400)

wilcox.test(NRRLssp$logCFU[NRRLssp$Strain=="PM.A1.5"], NRRLssp$logCFU[NRRLssp$Strain=="PM.A1.9"])
wilcox.test(NRRLssp$logCFU[NRRLssp$Strain=="PM.I2.5"], NRRLssp$logCFU[NRRLssp$Strain=="PM.I2.9"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figure 3 Community Invasion A1, I2
#### file contains raw CFU/worm for each species [AA, MO, PF, CG] + total 
#### and relative abundance [fAA etc] for each

NRRLcommInvade<-read.table("NRRLcommInvade.txt", header=TRUE)
NRRLcommInvade<-as_tibble(NRRLcommInvade)
NRRLcommInvade$Condition<-as.factor(NRRLcommInvade$Condition)
NRRLcommInvade$Invader<-as.factor(NRRLcommInvade$Invader)
names(NRRLcommInvade)
i<-dim(NRRLcommInvade)[1]

#~~~~ while we are here, let's see if CG is different in SSP vs communities (yes)
NRRLssp.CG<-subset(NRRLssp, Species=="CG")
unique(NRRLssp.CG$Strain)
unique(NRRLcommInvade$Community)

wilcox.test(NRRLssp.CG$CFU[NRRLssp.CG$Strain=="CG.A1.9"],
            NRRLcommInvade$CG[NRRLcommInvade$Community=="A1.9"])
hist(log10(NRRLssp.CG$CFU[NRRLssp.CG$Strain=="CG.A1.9"]+1), xlim=c(0,5), breaks=20)
hist(log10(NRRLcommInvade$CG[NRRLcommInvade$Community=="A1.9"]+1), xlim=c(0,5),breaks=20)
median(NRRLssp.CG$CFU[NRRLssp.CG$Strain=="CG.A1.9"])
median(NRRLcommInvade$CG[NRRLcommInvade$Community=="A1.9"])

wilcox.test(NRRLssp.CG$CFU[NRRLssp.CG$Strain=="CG.I2.9"],
            NRRLcommInvade$CG[NRRLcommInvade$Community=="I2.9"])
hist(log10(NRRLssp.CG$CFU[NRRLssp.CG$Strain=="CG.I2.9"]+1), xlim=c(0,5), breaks=20)
hist(log10(NRRLcommInvade$CG[NRRLcommInvade$Community=="I2.9"]+1), xlim=c(0,5), breaks=20)
median(NRRLssp.CG$CFU[NRRLssp.CG$Strain=="CG.I2.9"])
median(NRRLcommInvade$CG[NRRLcommInvade$Community=="I2.9"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assemble a data frame for stacked bars
NRRLcommInvadeStack<-data.frame(Condition= rep(NRRLcommInvade$Condition, 4),
                   Comm=rep(NRRLcommInvade$Community, 4), Invader=rep(NRRLcommInvade$Invader, 4), 
                   bact=as.factor(c(rep("AA",i), rep("MO",i), rep("CG",i), rep("PM",i))),
                   values=c(NRRLcommInvade$fAA,NRRLcommInvade$fMO,NRRLcommInvade$fCG,NRRLcommInvade$fPM)) 
idx<-rep(1:i)
NRRLcommInvadeStack$ID<-rep(idx,4)

pStackInvade.A1.9.AA0<-subset(NRRLcommInvadeStack, Comm=="A1.9" & Invader=="AA.0") %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'gold', 'yellowgreen', 'blue'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.position="none") +
  labs(title="A1.9 vs. AA.0")
pStackInvade.A1.9.AA0

pStackInvade.A1.9.AA3<-subset(NRRLcommInvadeStack, Comm=="A1.9" & Invader=="AA.3") %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'gold', 'yellowgreen', 'blue'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.position="none") +
  labs(title="A1.9 vs. AA.3")

pStackInvade.I2.9.AA0<-subset(NRRLcommInvadeStack, Comm=="I2.9" & Invader=="AA.0") %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'gold', 'yellowgreen', 'blue'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.position="none") +
  labs(title="I2.9 vs. AA.0")

pStackInvade.I2.9.AA3<-subset(NRRLcommInvadeStack, Comm=="I2.9" & Invader=="AA.3") %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'gold', 'yellowgreen', 'blue'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.position="none") +
  labs(title="I2.9 vs. AA.3")

pStackInvade.N.AA0<-subset(NRRLcommInvadeStack, Comm=="N" & Invader=="AA.0") %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'gold', 'yellowgreen', 'blue'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title = element_blank()) +
  labs(title="Ancestor vs. AA.0")
pStackInvade.N.AA0

pStackInvade.N.AA3<-subset(NRRLcommInvadeStack, Comm=="N" & Invader=="AA.3") %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'gold', 'yellowgreen', 'blue'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title = element_blank()) +
  labs(title="Ancestor vs. AA.3")
plot_grid(pStackInvade.A1.9.AA0, pStackInvade.I2.9.AA0, pStackInvade.N.AA0,
          pStackInvade.A1.9.AA3, pStackInvade.I2.9.AA3, pStackInvade.N.AA3,
          ncol=3, nrow=2, rel_widths = c(1, 1, 1.3))
ggsave("Fig3A_pStackCommunityInvadePM.png", width=10, height=7, units="in")

### Now plot out fAA with Invader as a factor
NRRLcommInvade03<-subset(NRRLcommInvade, Invader=="AA.0" | Invader=="AA.3")
NRRLcommInvade03$logAA<-log10(NRRLcommInvade03$AA+1)

pNRRLcommInvade.A1<-subset(NRRLcommInvade03, Community=="A1.9") %>%
  ggplot(aes(x=Invader, y=logAA, color=Invader))+
  geom_boxplot(fill=NA)+ ylim(-0.1,5)+
  geom_jitter(size=3, shape=16, position=position_jitter(0.05)) + theme_classic()+
  theme(axis.text=element_text(size=14), 
#        axis.text.x = element_blank(),
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title="A1.9", y=expression(log[10](AA/Worm)), x="")+
  stat_compare_means(label.y = 4.8) 
#pNRRLcommInvade.A1

pNRRLcommInvade.I2<-subset(NRRLcommInvade03, Community=="I2.9") %>%
  ggplot(aes(x=Invader, y=logAA, color=Invader))+
  geom_boxplot(fill=NA)+ ylim(-0.1,5)+
  geom_jitter(size=3, shape=16, position=position_jitter(0.05)) + theme_classic()+
  theme(axis.text=element_text(size=14), 
        #axis.text.x = element_blank(),
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position ="none") +
  labs(title="I2.9", y=expression(log[10](AA/Worm)), x="")+
  stat_compare_means(label.y = 4.8) 
#pNRRLcommInvade.I2

pNRRLcommInvade.N<-subset(NRRLcommInvade03, Community=="N") %>%
  ggplot(aes(x=Invader, y=logAA, color=Invader))+
  geom_boxplot(fill=NA)+ ylim(-0.1,5)+
  geom_jitter(size=3, shape=16, position=position_jitter(0.05)) + theme_classic()+
  theme(axis.text=element_text(size=14), 
        #axis.text.x = element_blank(),
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title="Ancestor", y=expression(log[10](AA/Worm)), x="")+
  stat_compare_means(label.y = 4.8) 
#pNRRLcommInvade.N

bottom_grid<-plot_grid(pNRRLcommInvade.A1, pNRRLcommInvade.I2, pNRRLcommInvade.N, ncol=3)
top_grid<-plot_grid(pStackInvade.A1.9.AA0, pStackInvade.I2.9.AA0, pStackInvade.N.AA0,
          pStackInvade.A1.9.AA3, pStackInvade.I2.9.AA3, pStackInvade.N.AA3,
          ncol=3, nrow=2, rel_widths = c(1, 1, 1.4))
plot_grid(top_grid, bottom_grid, ncol=1, rel_heights = c(1.8,1), labels="AUTO")
ggsave("Fig3_pStackCommunityInvadePM_AB.jpg", width=10, height=12, units="in", dpi=400)

#~~~~~~~~~~~~~    Statistical tests       ~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~ for community invasion data ~~~~~~~~~~~~~~~~~~~~~~
kruskal.test(AA~Community, data=NRRLcommInvade)
#Kruskal-Wallis chi-squared = 28.746, df = 2, p-value = 5.727e-07
wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Invader=="AA.0"],NRRLcommInvade$AA[NRRLcommInvade$Invader=="AA.1"])
wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Invader=="AA.0"],NRRLcommInvade$AA[NRRLcommInvade$Invader=="AA.2"])
wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Invader=="AA.0"],NRRLcommInvade$AA[NRRLcommInvade$Invader=="AA.3"])

wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="A1.9"],NRRLcommInvade$AA[NRRLcommInvade$Community=="I2.9"])
wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="A1.9"],NRRLcommInvade$AA[NRRLcommInvade$Community=="N"])
wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="N"],NRRLcommInvade$AA[NRRLcommInvade$Community=="I2.9"])

wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="A1.9"& NRRLcommInvade$Invader=="AA.0"],NRRLcommInvade$AA[NRRLcommInvade$Community=="N" & NRRLcommInvade$Invader=="AA.0"])
wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="A1.9"& NRRLcommInvade$Invader=="AA.3"],NRRLcommInvade$AA[NRRLcommInvade$Community=="N" & NRRLcommInvade$Invader=="AA.3"])
wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="I2.9"& NRRLcommInvade$Invader=="AA.0"],NRRLcommInvade$AA[NRRLcommInvade$Community=="N" & NRRLcommInvade$Invader=="AA.0"])
wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="I2.9"& NRRLcommInvade$Invader=="AA.3"],NRRLcommInvade$AA[NRRLcommInvade$Community=="N" & NRRLcommInvade$Invader=="AA.3"])

wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="N"& NRRLcommInvade$Invader=="AA.0"],NRRLcommInvade$AA[NRRLcommInvade$Community=="N" & NRRLcommInvade$Invader=="AA.3"])

wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="A1.9"& NRRLcommInvade$Invader=="AA.0"],NRRLcommInvade$AA[NRRLcommInvade$Community=="I2.9" & NRRLcommInvade$Invader=="AA.0"])
wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="A1.9"& NRRLcommInvade$Invader=="AA.3"],NRRLcommInvade$AA[NRRLcommInvade$Community=="I2.9" & NRRLcommInvade$Invader=="AA.3"])

wilcox.test(NRRLcommInvade$AA[NRRLcommInvade$Community=="A1.9"& NRRLcommInvade$Invader=="AA.0"],NRRLcommInvade$AA[NRRLcommInvade$Community=="A1.9" & NRRLcommInvade$Invader=="AA.3"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# some summary stats for 
names(NRRLcommInvade)
compare_means(PM~Community, data=NRRLcommInvade)
median(NRRLcommInvade$PM[NRRLcommInvade$Community=="A1.9"]) #360
median(NRRLcommInvade$PM[NRRLcommInvade$Community=="I2.9"]) #1800
median(NRRLcommInvade$fPM[NRRLcommInvade$Community=="A1.9"]) #0.0714
median(NRRLcommInvade$fPM[NRRLcommInvade$Community=="I2.9"]) #0.1956

median(NRRLssp$CFU[NRRLssp$Strain=="PM.A1.9"]) #4000
median(NRRLssp$CFU[NRRLssp$Strain=="PM.I2.9"]) #3300


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~  FIGURE S5 COMMUNITY INVASION BARPLOTS
#~~~  AND PM SWAPS     
# Some pretty rainbow barplots for the supplementary figure
# assemble the data frame
NRRLcommInvadeBars<-data.frame(Condition= rep(NRRLcommInvade$Condition, 4),
                                Comm=rep(NRRLcommInvade$Community, 4), Invader=rep(NRRLcommInvade$Invader, 4), 
                                bact=as.factor(c(rep("AA",i), rep("MO",i), rep("CG",i), rep("PM",i))),
                                values=c(log10(NRRLcommInvade$AA),
                                         log10(NRRLcommInvade$MO),
                                         log10(NRRLcommInvade$CG),
                                         log10(NRRLcommInvade$PM))) 
idx<-rep(1:i)
NRRLcommInvadeBars$ID<-rep(idx,4)
NRRLcommInvadeBars$values[!is.finite(NRRLcommInvadeBars$values)]<-0

fSize<-12
#i<-dim(subset(NRRLcommInvadeBars, Comm=="A1.9" & Invader=="AA.0"))[1]
#mycolors<-rainbow(i)

pNRRLcommInvadeBars.A1.0<-subset(NRRLcommInvadeBars, Comm=="A1.9" & Invader=="AA.0") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
#  scale_fill_manual(values=mycolors) +
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="A1.9 vs. AA.0", y=expression(log[10](CFU/Worm)), x="")
#pNRRLcommInvadeBars.A1.0

#i<-as.integer(dim(subset(NRRLcommInvadeBars, Comm=="A1.9" & Invader=="AA.3"))[1])
#mycolors<-rainbow(i)
pNRRLcommInvadeBars.A1.3<-subset(NRRLcommInvadeBars, Comm=="A1.9" & Invader=="AA.3") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
#  scale_fill_manual(values=mycolors) +
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="A1.9 vs. AA.3", y=expression(log[10](CFU/Worm)), x="")
#pNRRLcommInvadeBars.A1.3

pNRRLcommInvadeBars.A1.1<-subset(NRRLcommInvadeBars, Comm=="A1.9" & Invader=="AA.1") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  #  scale_fill_manual(values=mycolors) +
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="A1.9 vs. AA.1", y=expression(log[10](CFU/Worm)), x="")

pNRRLcommInvadeBars.A1.2<-subset(NRRLcommInvadeBars, Comm=="A1.9" & Invader=="AA.2") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  #  scale_fill_manual(values=mycolors) +
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="A1.9 vs. AA.2", y=expression(log[10](CFU/Worm)), x="")

pNRRLcommInvadeBars.I2.0<-subset(NRRLcommInvadeBars, Comm=="I2.9"& Invader=="AA.0") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="I2.9 vs. AA.0", y=expression(log[10](CFU/Worm)), x="")

pNRRLcommInvadeBars.I2.3<-subset(NRRLcommInvadeBars, Comm=="I2.9" & Invader=="AA.3") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="I2.9 vs. AA.3", y=expression(log[10](CFU/Worm)), x="")

pNRRLcommInvadeBars.N.0<-subset(NRRLcommInvadeBars, Comm=="N" & Invader=="AA.0") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  scale_x_discrete(labels=c("AA" = "AA", "CG" = "CG", "MO" = "MO", "PM"="SS"))+
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="Ancestor vs. AA.0", y=expression(log[10](CFU/Worm)), x="")

pNRRLcommInvadeBars.N.3<-subset(NRRLcommInvadeBars, Comm=="N" & Invader=="AA.3") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  scale_x_discrete(labels=c("AA" = "AA", "CG" = "CG", "MO" = "MO", "PM"="SS"))+
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="Ancestor vs. AA.3", y=expression(log[10](CFU/Worm)), x="")

pNRRLcommInvadeBars.N.1<-subset(NRRLcommInvadeBars, Comm=="N" & Invader=="AA.1") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  scale_x_discrete(labels=c("AA" = "AA", "CG" = "CG", "MO" = "MO", "PM"="SS"))+
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="Ancestor vs. AA.1", y=expression(log[10](CFU/Worm)), x="")

pNRRLcommInvadeBars.N.2<-subset(NRRLcommInvadeBars, Comm=="N" & Invader=="AA.2") %>%
  ggplot(aes(x=bact, y=values, fill=as.factor(ID)))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+
  scale_x_discrete(labels=c("AA" = "AA", "CG" = "CG", "MO" = "MO", "PM"="SS"))+
  theme(axis.text=element_text(size=fSize),
        axis.title=element_text(size=fSize, face="bold"), 
        plot.title=element_text(hjust=0.5, size=fSize),
        #legend.text=element_text(size=14),
        #legend.title=element_blank(),
        legend.position = "none") +
  labs(title="Ancestor vs. AA.2", y=expression(log[10](CFU/Worm)), x="")

plot_grid(pNRRLcommInvadeBars.A1.0, pNRRLcommInvadeBars.I2.0, pNRRLcommInvadeBars.N.0, 
          pNRRLcommInvadeBars.A1.1, NULL, pNRRLcommInvadeBars.N.1,
          pNRRLcommInvadeBars.A1.2, NULL, pNRRLcommInvadeBars.N.2,
          pNRRLcommInvadeBars.A1.3, pNRRLcommInvadeBars.I2.3, pNRRLcommInvadeBars.N.3,
          ncol=3, nrow=4)
#ggsave("FigS4_pNRRLcommInvadeBarPlotsRainbowPM.jpg", width=8, height=10, units="in", dpi=400)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2022-08: Invasion in PM-swapped communities
# "A1.PM.I2", for example, is made of community A1 strains, except has PM from community I2.9

NRRLcommInvadeSwap<-read.table("NRRLcommInvadeSwap.txt", header=TRUE)
NRRLcommInvadeSwap<-as_tibble(NRRLcommInvadeSwap)
NRRLcommInvadeSwap$Condition<-as.factor(NRRLcommInvadeSwap$Condition)
NRRLcommInvadeSwap$Invader<-as.factor(NRRLcommInvadeSwap$Invader)
names(NRRLcommInvadeSwap)
i<-dim(NRRLcommInvadeSwap)[1]

# assemble a data frame for stacked bars
NRRLcommInvadeStackSwap<-data.frame(Condition= rep(NRRLcommInvadeSwap$Condition, 4),
                                    Comm=rep(NRRLcommInvadeSwap$Community, 4), 
                                    Invader=rep(NRRLcommInvadeSwap$Invader, 4), 
                                    bact=as.factor(c(rep("AA",i), rep("MO",i), rep("CG",i), rep("PM",i))),
                                    values=c(NRRLcommInvadeSwap$fAA,NRRLcommInvadeSwap$fMO,NRRLcommInvadeSwap$fCG,NRRLcommInvadeSwap$fPM)) 
idx<-rep(1:i)
NRRLcommInvadeStackSwap$ID<-rep(idx,4)

top_grid_swap<-NRRLcommInvadeStackSwap %>%
  ggplot(aes(fill=bact, y=values, x=ID))+ 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c('red4', 'gold', 'yellowgreen', 'blue'))  + 
  ylab("Freq") + xlab("") + 
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.position="right",
        legend.title = element_blank(),
        strip.text = element_text(size = 14)) +
  facet_wrap(vars(Comm, Invader), scales="free_x", ncol=4)
top_grid_swap

# merge with Figure 3 data?
names(NRRLcommInvade03)
names(NRRLcommInvadeSwap)

temp<-NRRLcommInvade03 %>%
  filter(Community!="N") %>%
  rename("nAA"="AA", "nMO"="MO", "nPM"="PM", "nCG"="CG")
namelist<-names(temp)

NRRLcommInvadeSwap_trunc<-NRRLcommInvadeSwap %>%
  mutate(logAA=log10(nAA+1)) %>%
  select(namelist)
#names(NRRLcommInvadeSwap_trunc)

NRRLcommInvadeMerged<-rbind(temp, NRRLcommInvadeSwap_trunc)
rm(temp)

bottom_grid_swap<-NRRLcommInvadeMerged %>% 
  #mutate(logAA=log10(nAA+1)) %>%
  ggplot(aes(x=Invader, y=logAA, color=Invader))+
  geom_boxplot(fill=NA)+
  ylim(-0.1,5)+
  #scale_y_log10() +
  geom_jitter(size=3, shape=16, position=position_jitter(0.05)) + 
  theme_classic()+
  theme(axis.text=element_text(size=14), 
        #        axis.text.x = element_blank(),
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 14)) +
  labs(y=expression(log[10](AA/Worm)), x="")+
  facet_wrap(vars(Community), ncol=4) +
  stat_compare_means(label.y = 4.8)
bottom_grid_swap
#plot_grid(top_grid_swap, bottom_grid_swap, ncol=1, rel_heights = c(1.8,1), labels="AUTO")
#ggsave("FigSX_CommunityInvasionSwap.png", width=14, height=10, units="in", dpi=400)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot together
# Supplementary Figure S5

plot_invade_bars<-plot_grid(pNRRLcommInvadeBars.A1.0, pNRRLcommInvadeBars.I2.0, pNRRLcommInvadeBars.N.0, 
          pNRRLcommInvadeBars.A1.1, NULL, pNRRLcommInvadeBars.N.1,
          pNRRLcommInvadeBars.A1.2, NULL, pNRRLcommInvadeBars.N.2,
          pNRRLcommInvadeBars.A1.3, pNRRLcommInvadeBars.I2.3, pNRRLcommInvadeBars.N.3,
          ncol=3, nrow=4)

plot_grid(plot_invade_bars, top_grid_swap, bottom_grid_swap, ncol=1, rel_heights = c(3, 1.8, 1), labels="AUTO")
ggsave("FigS5_CommunityInvasion_Swap_Barplots.png", width=14, height=18, units="in", dpi=400)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compare all groups using total AA as the variable
kruskal.test(nAA~Condition, data=NRRLcommInvadeSwap) # Kruskal-Wallis chi-squared = 6.6241, df = 7, p-value = 0.469
kruskal.test(nAA~Community, data=NRRLcommInvadeSwap) # Kruskal-Wallis chi-squared = 4.3052, df = 3, p-value = 0.2303

NRRLcommInvadeMerged$Community<-as.factor(NRRLcommInvadeMerged$Community)
raov(logAA~Community+Invader, data=NRRLcommInvadeMerged)
# Robust ANOVA Table
#DF      RD Mean RD       F p-value
#Community          3 8.63756 2.87919 5.00607 0.00218
#Invader            1 0.81350 0.81350 1.41445 0.23546
#Community:Invader  3 0.83843 0.27948 0.48593 0.69236

kruskal.test(nAA~Condition, data=NRRLcommInvadeMerged) # Kruskal-Wallis chi-squared = 21.003, df = 7, p-value = 0.003765
kruskal.test(nAA~Community, data=NRRLcommInvadeMerged) # Kruskal-Wallis chi-squared = 18.177, df = 3, p-value = 0.0004043

# Summary stats
NRRLcommInvadeMerged %>%
  group_by(Condition) %>%
  summarize(medianAA = median(nAA),
            meanAA = mean(nAA),
            sdAA = sd(nAA))
#Condition medianAA meanAA  sdAA
#<fct>        <dbl>  <dbl> <dbl>
#  1 A1.AA0          10   278.  683.
#2 A1.AA3          10  1249. 2688.
#3 I2.AA0         400  4130. 8807.
#4 I2.AA3         400  2309. 5105.
#5 AI.AA0         140  1076. 1741.
#6 AI.AA3          80  1310. 2667.
#7 IA.AA0         180   892. 1930.
#8 IA.AA3         400   763. 1187.

NRRLcommInvadeMerged %>%
  group_by(Invader) %>%
  summarize(medianAA = median(nAA),
            meanAA = mean(nAA),
            sdAA = sd(nAA))

#   Invader medianAA meanAA  sdAA
# <fct>      <dbl>  <dbl> <dbl>
#  1 AA.0         180  1637. 4956.
#  2 AA.3         200  1381. 3184.

NRRLcommInvadeMerged %>%
  group_by(Community) %>%
  summarize(medianAA = median(nAA),
            meanAA = mean(nAA),
            sdAA = sd(nAA))

#   Community medianAA meanAA  sdAA
# <chr>        <dbl>  <dbl> <dbl>
# 1 A1.9            10   722  1932.
# 2 A1.PM.I2        80  1185. 2207.
# 3 I2.9           400  3290. 7339.
# 4 I2.PM.A1       360   814. 1511.

wilcox.test(NRRLcommInvadeMerged$logAA[NRRLcommInvadeMerged$Community=="A1.9"],
            NRRLcommInvadeMerged$logAA[NRRLcommInvadeMerged$Community=="A1.PM.I2"]) #W = 1709, p-value = 0.03018
wilcox.test(NRRLcommInvadeMerged$logAA[NRRLcommInvadeMerged$Community=="I2.9"],
            NRRLcommInvadeMerged$logAA[NRRLcommInvadeMerged$Community=="I2.PM.A1"]) #W = 2051, p-value = 0.3999

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~ FIGURE 4
#~~~~~ Single-species invasions in the worm host (wild-type N2 ancestor)

NRRLsspInvade<-read.table("NRRLsspInvade.txt", header=TRUE) 
names(NRRLsspInvade)
unique(NRRLsspInvade$Sp1)
NRRLsspInvade$Sp1<-as.factor(NRRLsspInvade$Sp1)
NRRLsspInvade$Invader<-as.factor(NRRLsspInvade$Invader)
NRRLsspInvade$Group<-as.factor(NRRLsspInvade$Group)

NRRLsspInvade$logCount1<-log10(NRRLsspInvade$Count1+1) 
NRRLsspInvade$logCount2<-log10(NRRLsspInvade$Count2+1) 
NRRLsspInvade$logTotal<-log10(NRRLsspInvade$Total+1) 

# NRRLsspInvade<-subset(NRRLsspInvade, Sp1!="SS.0")
NRRLsspInvade<-as_tibble(NRRLsspInvade)
unique(NRRLsspInvade$Sp1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#add the single species colonization data for AA.0 and AA.3
temp<-subset(NRRLssp, Strain=="AA.0" | Strain=="AA.3")
temp<-as_tibble(temp)
temp$Rep<-temp$Date
names(NRRLsspInvade)
unique(NRRLsspInvade$Sp1)

mydates<-unique(temp$Date)
within(temp, Rep[Rep==mydates[1]] <- 1)
within(temp, Rep[Rep==mydates[2]] <- 2)
within(temp, Rep[Rep==mydates[3]] <- 3)

NRRLsspInvade<-rbind(NRRLsspInvade, tibble(Sp1="SSP",
                                           Invader=temp$Strain,
                                           Count1=0,
                                           Count2=temp$CFU,
                                           Total=temp$CFU,
                                           fAA=1,
                                           Date=temp$Date,
                                           Group="SSP",
                                           Rep=temp$Rep,
                                           Pass=-1,
                                           logCount1=0,
                                           logCount2=temp$logCFU,
                                           logTotal=temp$logCFU
))
rm(temp)

#~~~~~~~~~~~~~~~ FIGURE 4: plot vs SSP colonization
pFig4A_NRRLsspInvade_vs_SSP<-NRRLsspInvade %>%
  ggplot(aes(x=Sp1, y=logCount2, color=Sp1)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  stat_summary(fun = "median", geom = "point", 
               shape = 15, size = 3, color = "black")+
  theme_classic() + 
  ylim(-0.1,6.5) +
  geom_hline(yintercept = log10(20), linetype="dashed", color="gray")+
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(hjust=0.5,size=14),
        legend.position="none") +
  labs(y=expression(log[10](AA)), x="")+
  stat_compare_means(label.y = 4.8, label.x=8) +
  #stat_compare_means(aes(label=..p.format..), ref.group = "SSP", label.y=5.8) + 
  stat_compare_means(label="p.signif", ref.group = "SSP", label.y=6.3) +
  stat_compare_means(label="p.signif", ref.group = "OP50", label.y=5.8) +
  facet_wrap(~Invader, ncol=1)
pFig4A_NRRLsspInvade_vs_SSP

# Just to see all the test results
temp<-compare_means(logCount2~Sp1, data=subset(NRRLsspInvade, Invader=="AA.0"))
View(temp)
temp<-compare_means(logCount2~Sp1, data=subset(NRRLsspInvade, Invader=="AA.3"))
View(temp)
rm(temp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~ Plot AA.0 vs AA.3 within each pre-colonization condition
pFig4B_NRRLsspInvade_vs_AA<-NRRLsspInvade %>%
  ggplot(aes(x=Invader, y=logCount2, color=Invader)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  stat_summary(fun = "median", geom = "point", 
               shape = 15, size = 3, color = "black")+
  theme_classic() + 
  ylim(-0.1,6) +
  geom_hline(yintercept = log10(20), linetype="dashed", color="gray")+
  theme(axis.text=element_text(size=12), 
        #axis.text.x = element_blank(),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(hjust=0.5,size=14),
        legend.position="none") +
  labs(y=expression(log[10](AA)), x="")+
  stat_compare_means(label.y = 5.8, size=3) +
  facet_wrap(~Sp1, ncol=7)
pFig4B_NRRLsspInvade_vs_AA

#sspN2InvadeLegend<-subset(NRRLsspInvade, Sp1=="PM.I2.9") %>%
#  ggplot(aes(x=Invader, y=logCount2, color=Invader)) + 
#  geom_violin(fill=NA) + theme_classic() + 
#  geom_jitter(shape=16)+
#  theme(legend.text=element_text(size=14), legend.position = "right", legend.title = element_blank())

#library(gridExtra)
#extract_legend<-function(myggplot){
#  tmp <- ggplot_gtable(ggplot_build(myggplot))
#  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#  legend <- tmp$grobs[[leg]]
#  return(legend)
#}
#legend<-extract_legend(sspN2InvadeLegend)
#blank <- grid.rect(gp=gpar(col="white"))

# Plot out all single-species invasions
#grid.arrange(
#  arrangeGrob(
#              psspN2InvadeCG.0, psspN2InvadeCG.A1.9, psspN2InvadeCG.I2.9,
#              psspN2InvadeMO.0, psspN2InvadeMO.A1.9, psspN2InvadeMO.I2.9,
#              psspN2InvadePM.A1.5, psspN2InvadePM.A1.9, blank,
#              psspN2InvadePM.I2.5, psspN2InvadePM.I2.9, psspN2InvadeOP50,
#              ncol=3),
#  legend, ncol=2, widths=c(10,1.2))

#ggsave("psspN2Invade_worms_A1I2.png", width=9, height=9, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~ Save Figure 4
plot_grid(pFig4A_NRRLsspInvade_vs_SSP, pFig4B_NRRLsspInvade_vs_AA, ncol=1, rel_heights = c(1,1.2), labels = "AUTO")
ggsave("Fig4_NRRL_ssp_invade_worms.jpg", width=10, height=10, units="in", dpi=400)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# and now some comparisons
kruskal.test(logCount2~Invader, data=NRRLsspInvade)
#Kruskal-Wallis chi-squared = 9.6855, df = 1, p-value = 0.001857
pairwise.wilcox.test(NRRLsspInvade$logCount2, NRRLsspInvade$Invader, p.adjust.method = "bonf")
#AA0 vs AA3 0.0019
kruskal.test(fAA~Invader, data=NRRLsspInvade)
#Kruskal-Wallis chi-squared = 7.433, df = 1, p-value = 0.006404
kruskal.test(fAA~Sp1, data=NRRLsspInvade)
#Kruskal-Wallis chi-squared = 242.38, df = 11, p-value < 2.2e-16
kruskal.test(logCount2~Sp1, data=NRRLsspInvade)
#Kruskal-Wallis chi-squared = 155.13, df = 11, p-value < 2.2e-16
# for Table S2A
pairwise.wilcox.test(NRRLsspInvade$logCount2, NRRLsspInvade$Sp1, p.adjust.method = "bonf" )

NRRLsspInvade.AA3<-subset(NRRLsspInvade, Invader=="AA.3")
# for Table S2C
pairwise.wilcox.test(NRRLsspInvade.AA3$logCount2, NRRLsspInvade.AA3$Sp1, p.adjust.method = "bonf" )

NRRLsspInvade.AA0<-subset(NRRLsspInvade, Invader=="AA.0")
kruskal.test(logCount2~Sp1, data=NRRLsspInvade.AA0)
#Kruskal-Wallis chi-squared = 107.16, df = 11, p-value < 2.2e-16
# for Table S2B
pairwise.wilcox.test(NRRLsspInvade.AA0$logCount2, NRRLsspInvade.AA0$Sp1, p.adjust.method = "bonf" )

wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="CG.0" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="CG.0" & NRRLsspInvade$Invader=="AA.3"])
#W = 210, p-value = 0.3624
wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="CG.A1.9" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="CG.A1.9" & NRRLsspInvade$Invader=="AA.3"])
#0.849
wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="CG.I2.9" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="CG.I2.9" & NRRLsspInvade$Invader=="AA.3"])
#0.4
wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="MO.0" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="MO.0" & NRRLsspInvade$Invader=="AA.3"])
#0.74
wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="MO.A1.9" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="MO.A1.9" & NRRLsspInvade$Invader=="AA.3"])
#0.6
wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="MO.A1.9" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="MO.A1.9" & NRRLsspInvade$Invader=="AA.3"])
#W = 233, p-value = 0.6468
wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="PM.A1.5" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="PM.A1.5" & NRRLsspInvade$Invader=="AA.3"])
#W = 30.5, p-value = 0.02452
wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="PM.A1.9" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="PM.A1.9" & NRRLsspInvade$Invader=="AA.3"])
#W = 299, p-value = 1.369e-05
wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="PM.I2.9" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="PM.I2.9" & NRRLsspInvade$Invader=="AA.3"])
#W = 184, p-value = 0.0177
wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="PM.I2.5" & NRRLsspInvade$Invader=="AA.0"],
            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="PM.I2.5" & NRRLsspInvade$Invader=="AA.3"])
#W = 90.5, p-value = 0.2975
#wilcox.test(NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="SS.0" & NRRLsspInvade$Invader=="AA.0"],
#            NRRLsspInvade$logCount2[NRRLsspInvade$Sp1=="SS.0" & NRRLsspInvade$Invader=="AA.3"])
#W = 115, p-value = 0.01365


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### Figure S6
####### And the in vitro data from 3-4 reps
####### plus the higher two AA concentration conditions from density-dep vs PF 
# we clearly have one run of CG that's acting out (run 3, 8/3/2021) 
# this run did start with about 100X lower GC than AA; will exclude on that basis 

NRRLsspInVitro<-read.table("NRRLsspInVitro.txt", header=TRUE)  
unique(NRRLsspInVitro$Sp1)
NRRLsspInVitro$Sp1<-as.factor(NRRLsspInVitro$Sp1) 
NRRLsspInVitro$Sp2<-as.factor(NRRLsspInVitro$Sp2) 
names(NRRLsspInVitro)
pNRRLsspInVitro.1<-subset(NRRLsspInVitro, Sp1=="CG.0"|Sp1=="CG.A1.9"|Sp1=="MO.0"|Sp1=="MO.A1.9"|Sp1=="OP50") %>%
  ggplot(aes(x=Sp1, y=P2, color=Sp2, shape=Sp2)) +
  geom_jitter(size=4, width=0.1) + theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = c(0.8,0.25)) +
  labs(title="In Vitro, AA vs. Community A1", y="AA Frequency", x="")
pNRRLsspInVitro.1
pNRRLsspInVitro.PM<-subset(NRRLsspInVitro, Sp1=="PM.A1.5"|Sp1=="PM.A1.9"|Sp1=="PM.I2.5"|Sp1=="PM.I2.9"|Sp1=="OP50") %>%
  ggplot(aes(x=Sp1, y=P2, color=Sp2, shape=Sp2)) +
  geom_jitter(size=4, width=0.1) + theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title="In Vitro, AA vs. PM", y="AA Frequency", x="")
#pNRRLsspInVitro.PM
plot_grid(pNRRLsspInVitro.1, pNRRLsspInVitro.PM, ncol=1, labels="AUTO")
ggsave("pNRRLInVitroCompete.png", width=5, height=7, units="in")
# NOTE WILL MERGE THIS WITH FREQUNCY DEPENDENCE

## frequency dependence in AA vs PM?
NRRLsspInVitroPMRatios<-read.table("NRRLsspInVitroPMRatios.txt",header=TRUE)
NRRLsspInVitroPMRatios$Sp1<-as.factor(NRRLsspInVitroPMRatios$Sp1)
NRRLsspInVitroPMRatios$Sp2<-as.factor(NRRLsspInVitroPMRatios$Sp2)
names(NRRLsspInVitroPMRatios)
pInVitroAA0<-subset(NRRLsspInVitroPMRatios, Sp2=="AA.0") %>%
  ggplot(aes(x=p0, y=pF, color=Sp1, shape=Sp1)) +
  geom_jitter(size=3, width=0.05) + theme_classic() +
  xlim(-0.05,1.05)+ylim(-0.05,1.05)+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.position = c(0.8,0.2)) +
  labs(title="AA.0")
pInVitroAA0
pInVitroAA3<-subset(NRRLsspInVitroPFRatios, Sp2=="AA.3") %>%
  ggplot(aes(x=p0, y=pF, color=Sp1, shape=Sp1)) +
  geom_jitter(size=3, width=0.05) + theme_classic() +
  xlim(-0.05,1.05)+ylim(-0.05,1.05)+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title="AA.3")
pInVitroAA3  
plot_grid(pInVitroAA0, pInVitroAA3, ncol=1, labels="AUTO")
ggsave("pInVitroAARatio.png", width=5, height=7, units="in")
pInVitroTop<-plot_grid(pNRRLsspInVitro.1, pNRRLsspInVitro.PM, pInVitroAA0, pInVitroAA3, ncol=2, labels="AUTO")
#ggsave("pInVitroAA_ALL.png", width=10, height=10, units="in")

kruskal.test(N2~Sp2, data=NRRLsspInVitro)
kruskal.test(N2~Sp1, data=NRRLsspInVitro)
pairwise.wilcox.test(NRRLsspInVitro$N2, NRRLsspInVitro$Sp1, p.adjust.method = "bonf")

## 2020-04-15 CG ancestor vs AA ancestor, time series
NRRL_AA_CG_InVitroTime<-read.table("NRRL_AA_CG_InVitroTime.txt", header=TRUE)  
unique(NRRL_AA_CG_InVitroTime$Sp1)
NRRL_AA_CG_InVitroTime$Sp1<-as.factor(NRRL_AA_CG_InVitroTime$Sp1) 
NRRL_AA_CG_InVitroTime$Sp2<-as.factor(NRRL_AA_CG_InVitroTime$Sp2) 
names(NRRL_AA_CG_InVitroTime)
pNRRL_AA_CG_InVitroTime.P1<-NRRL_AA_CG_InVitroTime %>%
  ggplot(aes(x=week, y=P1, color=Sp1, shape=Sp1)) +
  geom_jitter(size=4, width=0.1) + theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title="In Vitro, AA vs. CG", y="AA Frequency", x="Week")
pNRRL_AA_CG_InVitroTime.P1
pNRRL_AA_CG_InVitroTime.N1<-NRRL_AA_CG_InVitroTime %>%
  ggplot(aes(x=week, y=logN1, color=Sp1, shape=Sp1)) +
  geom_jitter(size=4, width=0.1) + theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = "none") +
  labs(title="In Vitro, AA vs. CG", y="AA log10(CFU/Plate)", x="Week")
pNRRL_AA_CG_InVitroTime.N1
pNRRL_AA_CG_InVitroTime.N2<-NRRL_AA_CG_InVitroTime %>%
  ggplot(aes(x=week, y=logN2, color=Sp1, shape=Sp1)) +
  geom_jitter(size=4, width=0.1) + theme_classic() +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14, face="bold"), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        legend.position = c(0.8,0.25)) +
  labs(title="In Vitro, AA vs. CG", y="CG log10(CFU/Plate)", x="Week")
pNRRL_AA_CG_InVitroTime.N2
pInVitroBottom<-plot_grid(pNRRL_AA_CG_InVitroTime.P1, pNRRL_AA_CG_InVitroTime.N1, pNRRL_AA_CG_InVitroTime.N2, align="h", ncol=3, labels=c("E","F","G"))

pInVitroAll<-plot_grid(pInVitroTop, pInVitroBottom, ncol=1, rel_heights = c(2,1))
pInVitroAll
ggsave("FigS6_pInVitroAA_ALL.jpg", width=10, height=14, units="in", dpi=400)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exporting in vivo CFU data to a file for other processing
NRRL1counts<-tibble(Condition=c(as.character(NRRLssp$Strain), as.character(NRRLcommInvade$Condition)),
                    CFU=c(NRRLssp$CFU, NRRLcommInvade$Total),
                    logCFU=c(NRRLssp$logCFU, log10(NRRLcommInvade$Total)),
                    Date=c(NRRLssp$Date, NRRLcommInvade$Date)
                    )
NRRL1counts
sum(is.na(NRRL1counts))
write.table(NRRL1counts, "NRRL1counts.csv")
