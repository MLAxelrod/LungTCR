##LUNG TCR Analysis
library(immunarch)
library(tidyr)
library(vegan)
library(ggplot2)
library(rstatix)
library(precrec)
library(survival)
library(survminer)
library(plyr)
library(ggpubr)
library(cutpointr)
library(ggrepel)

##ggplot parameters
theme_set(theme_classic()) 
theme_update(axis.text=element_text(color="black", size=10))


######## Immunarch analyses - TRBV gene Usage###########
immdata<-repLoad("C:/Users/marga/Desktop/TCRs_LungCancer/TRB_freq/")

##Basic QC
vol<-repExplore(immdata$data, .method="clones")
vol2<-repExplore(immdata$data, .method="volume")

p1<-ggplot(vol, aes(Sample, Clones))+
  geom_col()+
  geom_label(data=vol, aes(label=Clones))+
  labs(y="Total TCR reads", x=element_blank())+
  theme(axis.text.x = element_blank())

p2<-ggplot(vol2, aes(Sample, Volume))+
  geom_col()+
  geom_label(data=vol2, aes(label=Volume))+
  labs(y="Unique TCR reads", x=element_blank())+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

p1/p2

####FIGURE 1 ####
imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T, .type="family")

vis(imm_gu, .by = "Response", .meta = immdata$meta, .plot="box")+
  labs(x="Gene", title=element_blank())+
  scale_fill_manual(values = c("firebrick1", "cyan4"))

###Figure S1 ####
vis(imm_gu, .by = "Stage of disease at initial diagnosis (Limited = 0, Extensive = 1)", .meta = immdata$meta, .plot="box")+
  labs(x="Gene Family", title=element_blank(), fill="Stage at\ndiagnosis")+
  scale_fill_manual(values = c("firebrick1", "cyan4"), labels=c("Limited", "Extensive"))



#### TCR repertoire Diversity Analyses using Shannon Diversity Index ####
setwd("C:/Users/marga/Desktop/TCRs_LungCancer/")
meta<-read.csv("meta.csv")

####Calculate Shannon Diversity TRB ####

##These are directories where the TCR frequency tables are stored
files<- list.files("C:/Users/marga/Desktop/TCRs_LungCancer/TRB_freq/",  full.names=TRUE)[1:41] #1:41 here is to not include metadata file which I put in that directory only for immunarch
shortnames<- list.files("C:/Users/marga/Desktop/TCRs_LungCancer/TRB_freq/",  full.names=FALSE)[1:41]

#calculate diversity indices for each file
result<-data.frame()
for (i in seq_along(files)) {
  tmp<-read.table(files[i], sep="\t", header=TRUE)
  tmp2<-data.frame(sample=shortnames[i], shannon=diversity(tmp[6], index="shannon"), simpson=diversity(tmp[6], index="simpson"), inv=diversity(tmp[6], index="invsimpson"))
  result<-rbind(result, tmp2)
}

all(meta$file==result$sample) #all true
meta2<-cbind(meta,result) #add Shannon diversity to metdata


#how diversity indices are related
ggplot(result, aes(shannon, simpson))+geom_point()
ggplot(result, aes(shannon, inv))+geom_point()
ggplot(result, aes(simpson, inv))+geom_point()


####Figure 2####
f2<-meta2[meta2$Fig2Groups==1,]

t_test(f2, shannon~Response)

ggplot(f2, aes(Response, shannon))+geom_boxplot(aes(fill=Response))+geom_point()+
  scale_fill_manual(values = c("firebrick1", "cyan4"))+
  labs(y="TCR Diversity")+theme(legend.position = "none")

ggplot(f2, aes(x=reorder(patient.ID,shannon), y=shannon)) + geom_col(aes(fill=Response))+
  labs(y="TCR Diversity", x="Patient")+
  theme(axis.text.x = element_text(angle=90))+
  scale_fill_manual(values = c("firebrick1", "cyan4"))



####Analysis group for figure 3###
f3<-meta2[meta2$Fig3Groups==1,]

#make short v long PFS and OS groups
f3<-mutate(f3, PFSgroup=ifelse(f3$PFS..months.>9, "long", "short"))
f3<-mutate(f3, OSgroup=ifelse(f3$OS..months.>16, "long", "short"))
table(f3$PFSgroup)

t_test(f3, shannon~PFSgroup)
t_test(f3, shannon~OSgroup)

p1<-ggplot(f3, aes(PFSgroup, shannon))+geom_boxplot()+geom_point()+
  labs(y="Shannon Diversity", x="PFS Group")+
  scale_x_discrete(limits=c("short", "long"))
p2<-ggplot(f3, aes(OSgroup, shannon))+geom_boxplot()+geom_point()+
  labs(y="Shannon Diversity", x="OS Group")+
  scale_x_discrete(limits=c("short", "long"))
p1|p2


###FIGURE 3 A, B ####
#ROC analysis
precrec_obj <- evalmod(scores = f3$shannon, labels = f3$PFSgroup, posclass = "long")
precrec::auc(precrec_obj)
p3<-autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.75, y=0.05, label="AUC: 0.44")+
  theme(axis.text = element_text(color="black"),
        axis.title = element_text(size=14), plot.title = element_text(hjust=0.5))+
  labs(title="TCR Diversity as a\nClassifier of PFS")


precrec_obj <- evalmod(scores = f3$shannon, labels = f3$OSgroup, posclass = "long")
precrec::auc(precrec_obj)
p4<-autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.75, y=0.05, label="AUC: 0.89")+
  theme(axis.text = element_text(color="black"),
        axis.title = element_text(size=14), plot.title = element_text(hjust=0.5))+
  labs(title="TCR Diversity as a\nClassifier of OS")

p3|p4

##Youden index to define optimal cutpoint for OS and PFS
cp<-cutpointr(f3, shannon, OSgroup, 
          method = maximize_metric, metric = youden, pos_class= "long", direction=">=")

summary(cp)
##Youden optimal cutpoint of 7.7821

cp2<-cutpointr(f3, shannon, PFSgroup, 
              method = maximize_metric, metric = youden, pos_class= "long", direction=">=")
summary(cp2)


##Apply Youden cutoff based on OS grouping
f3<-mutate(f3, ShannonGroup=ifelse(f3$shannon>7.7821, "High", "Low"))
table(f3$ShannonGroup)

###########FIGURE 3 B,C ################

fit<-survfit(Surv(PFS..months., Progression...1..progression.event..0..none.) ~ ShannonGroup, data=f3)
print(fit)

ggsurvplot(fit,
           data=f3,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE,
           break.x.by=5,
           legend.labs=c("High TCR Diversity", "Low TCR Diversity"),
           palette = c("cyan4","firebrick1"))+
  ylab("PFS Probability")+xlab("Time (months)")


fit2<-survfit(Surv(OS..months., Death..1..death.event..0..alive.) ~ ShannonGroup, data=f3)
print(fit2)

ggsurvplot(fit2,
           data=f3,
           pval = TRUE, pval.coord=c(40,0.1),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           censor.shape=124,
           ggtheme = theme_classic(),
           font.y=14, font.x=14, font.tickslab="black",
           legend.title=element_blank(),
           tables.theme = theme_cleantable(),
           fontsize=4,tables.y.text=FALSE,
           break.x.by=5,
           legend.labs=c("High TCR Diversity", "Low TCR Diversity"),
           palette = c("cyan4","firebrick1"))+
  ylab("OS Probability")+xlab("Time (months)")


