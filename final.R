library(survminer)
library(survival)
library(MASS)
library(enrichR)
# Load data
######Reading the files into the program


 thymoma = read.table("thymoma_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Adrenocortical = read.table("Adrenocortical_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Sarcoma = read.table("Sarcoma_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Prostate = read.table("Prostate_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Breast = read.table("Breast_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Pancreatic = read.table("Pancreatic_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Testicular = read.table("Testicular_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 LungAdeno = read.table("LungAdeno_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Endometroid = read.table("Endometroid_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Liver = read.table("Liver_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 KidneyChromophobe = read.table("KidneyChromophobe_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Kidney_Clear_Cell = read.table("Kidney_Clear_Cell_Carcinoma_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Glioma = read.table("Glioma_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Lower_grade_Glioma = read.table("Lower_grade_Glioma_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 BileDuct=read.table("BileDuct_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Bladder=read.table("Bladder_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Head_Neck=read.table("Head_Neck_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Kidney_Papillary = read.table("Kidney_Papillary_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 LargeBcellLymphoma = read.table("LargeBcellLymphoma_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 OcularMelanoma = read.table("Ocular_Melanoma_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Stomach= read.table("Stomach_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Thyroid_acta2= read.table("Thyroid_trak1_ACTA2.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Thyroid_col1a1= read.table("Thyroid_trak1_COL1A1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 
 Rectal = read.table("Rectal_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Melanoma = read.table("Melanoma_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 Mesothelioma = read.table("Mesothelioma_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
 

# overlap=duplicated(c(Sarcoma$V5,Adrenocortical$V5))
#
# overlap=Sarcoma$V5 %in% Adrenocortical$V5
# overlap1=intersect(Sarcoma$V5,Adrenocortical$V5)
#
# ss=Reduce(intersect, list(Liver$V5,Lower_grade_Glioma$V5,Pancreatic$V5,LungAdeno$V5))
ss1=Reduce(intersect, list(Liver$V5,Lower_grade_Glioma$V5,Pancreatic$V5,LungAdeno$V5,Sarcoma$V5,Adrenocortical$V5,thymoma$V5,Prostate$V5,Breast$V5,Testicular$V5,Endometroid$V5,KidneyChromophobe$V5,Kidney_Clear_Cell$V5,Glioma$V5,BileDuct$V5))
#
ss2=Reduce(intersect, list(Thyroid$V5,Stomach$V5,OcularMelanoma$V5,LargeBcellLymphoma$V5,Kidney_Papillary$V5,Head_Neck$V5,Bladder$V5,Liver$V5,Lower_grade_Glioma$V5,Pancreatic$V5,LungAdeno$V5,Sarcoma$V5,Adrenocortical$V5,thymoma$V5,Prostate$V5,Breast$V5,Testicular$V5,Endometroid$V5,KidneyChromophobe$V5,Kidney_Clear_Cell$V5,Glioma$V5,BileDuct$V5))

ss2_all=Reduce(intersect, list(Mesothelioma$V5, Melanoma$V5,Rectal$V5,Thyroid$V5,Stomach$V5,OcularMelanoma$V5,LargeBcellLymphoma$V5,Kidney_Papillary$V5,Head_Neck$V5,Bladder$V5,Liver$V5,Lower_grade_Glioma$V5,Pancreatic$V5,LungAdeno$V5,Sarcoma$V5,Adrenocortical$V5,thymoma$V5,Prostate$V5,Breast$V5,Testicular$V5,Endometroid$V5,KidneyChromophobe$V5,Kidney_Clear_Cell$V5,Glioma$V5,BileDuct$V5))

ss3=Reduce(intersect, list(Thyroid$V5,Kidney_Papillary$V5,Bladder$V5,Liver$V5,Lower_grade_Glioma$V5,Pancreatic$V5,LungAdeno$V5,Sarcoma$V5,Adrenocortical$V5,thymoma$V5,Prostate$V5,Breast$V5,Testicular$V5,Endometroid$V5,KidneyChromophobe$V5,Kidney_Clear_Cell$V5,Glioma$V5,BileDuct$V5))

ss4=Reduce(intersect, list(Thyroid$V5,OcularMelanoma$V5,Kidney_Papillary$V5,Bladder$V5,Liver$V5,Lower_grade_Glioma$V5,Pancreatic$V5,LungAdeno$V5,Sarcoma$V5,Adrenocortical$V5,thymoma$V5,Prostate$V5,Breast$V5,Testicular$V5,Endometroid$V5,KidneyChromophobe$V5,Kidney_Clear_Cell$V5,Glioma$V5,BileDuct$V5))

ss4=Reduce(intersect, list(Thyroid$V5,Stomach$V5,OcularMelanoma$V5,LargeBcellLymphoma$V5,Kidney_Papillary$V5,Head_Neck$V5,Bladder$V5,Liver$V5,Lower_grade_Glioma$V5,Pancreatic$V5,LungAdeno$V5,Sarcoma$V5,Adrenocortical$V5,thymoma$V5,Prostate$V5,Breast$V5,Testicular$V5,Endometroid$V5,KidneyChromophobe$V5,Kidney_Clear_Cell$V5,Glioma$V5,BileDuct$V5))

ss1=Reduce(intersect, list(Liver$V5,Lower_grade_Glioma$V5,Pancreatic$V5,LungAdeno$V5,Sarcoma$V5,Adrenocortical$V5,thymoma$V5,Breast$V5,KidneyChromophobe$V5,Kidney_Clear_Cell$V5,Glioma$V5,BileDuct$V5))
#
ss1=Reduce(intersect, list(Liver$V5,Lower_grade_Glioma$V5,Pancreatic$V5,LungAdeno$V5,Sarcoma$V5,Adrenocortical$V5,thymoma$V5,Breast$V5,KidneyChromophobe$V5))

kidney=Reduce(intersect, list(Kidney_Clear_Cell$V5,Kidney_Papillary$V5,KidneyChromophobe$V5))




rm(list=ls()) #remove any previous variables
setwd("/Applications/UPenn/Summer2020/Melanoma")
df.input = read.table("SKCM_Survival.txt", header = TRUE,sep="\t",stringsAsFactors = FALSE) #Phenotype data like survival years etc
df.map = read.table("Melanoma_Seq", header = FALSE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA")) #Matrix of genes vs patient id
df.map <- na.omit(df.map) #remove NA if any
df.strong_scaling = read.table("Pancreatic_ACTA2_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
phenotype = read.delim("clinical.tsv", header = TRUE,sep="\t",stringsAsFactors = FALSE) #Phenotype data like survival years etc

## narrowing down to only primary tumor sites
Trk1=df.map[2:dim(df.map)[1],2:3] #data frame storing mRNA seq values only
name=df.input[1] #patient ids from the survival file
namepro=df.map[1,2:dim(df.map)[2]] #patient ids from mrna seq file
patientid=substr(namepro[1],1,12) #pure patient id
k=1
for(i in 1:dim(namepro)[2])
{
  if(substr(namepro[1,i],14,15)=='06')
  {
    Trk1[,k]=as.double(df.map[2:dim(df.map)[1],i+1]) #mrna reads
    patientid[k]=substr(namepro[i],1,12)
    k=k+1
  }
}

##matching RNA seq patient data to survival data
Trk <- data.frame(matrix(ncol = 3, nrow = 2))
for(i in 1:length(patientid))
{
  for(j in 1:dim(df.input)[1])
  {
    if(patientid[i]==df.input$X_PATIENT[j] )
    {
      Trk[i,3]=df.input$OS[j]
      Trk[i,2]=df.input$OS.time[j]/365
      break
    }
  }
}

Trk <- data.frame(matrix(ncol = 3, nrow = 2))
for(i in 1:length(patientid))
{
  for(j in 1:dim(df.input)[1])
  {
    if(patientid[i]==df.input$X_PATIENT[j] )
    {
      Trk[i,3]=df.input$OS[j]
      Trk[i,2]=df.input$OS.time[j]/365
      break
    }
  }
}

chk=0.75*dim(Trk1)[2] 
kmsig <- data.frame(matrix(ncol = 3, nrow = 2))
pro_name=df.map[1,1]
k=1
for (i in 1:dim(Trk1)[1])
{
  if(sum(Trk1[i,]==0)<chk)
  {
    Trk[,1]=as.double(t(Trk1[i,]))
    Trk[,4]=NULL
    Trk[,4] <- Trk[,1]> median(Trk[,1])
    fit <- survfit(Surv(Trk[,2],Trk[,3]) ~ Trk[,4],
                   data = Trk)
    
    kmsig[k,1]=surv_pvalue(fit)$pval
    kmsig[k,2]=surv_median(fit)$median[1]
    kmsig[k,3]=surv_median(fit)$median[2]
    pro_name[k]=df.map[i+1,1]
    k=k+1
  }
}

kmsig1=kmsig
kmsig2=cbind(kmsig1,pro_name)
kmsig2[,5]=kmsig2[,3]/kmsig2[,2]

#Changing flag to one if the ratio of survival is N/A(not applicable) for half the genes that kmsig2 has
flag=0
if(sum(is.na(kmsig2[,5]))>0.5*dim(kmsig2)[1])
{
  flag=1
}

write.table(kmsig2, file = "TumorName_Survival", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "") #Change Tumor Name here

#-----------------If you have already created and written TumorName_Survival, you do not need to run the code above this, can start running the code from here (uncomment kmsig2 and df.strong_scaling)
kmsig2=read.table("Kidney_Papillary_Survival", header = TRUE,sep=" ",stringsAsFactors = FALSE)

kmsig3=subset(kmsig2,kmsig2[,1]<0.05 & kmsig2[,5]>1)
kmsig4=subset(kmsig3,kmsig3$pro_name %in% df.strong_scaling$V5)
df.strong_scaling = read.table("Kidney_Papillary_trak1.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
volcano = read.table("Trk_BRCA_ACTA2.txt", header = FALSE,sep=",",stringsAsFactors = FALSE)
volcano1= subset(volcano,abs(volcano[,1])<0.2 & volcano[,3]>0.2)
volcano[-volcano1,]


LMNB1_km=subset(kmsig2,kmsig2$pro_name %in% df.strong_scaling[,5])

sig=subset(LMNB1_km, LMNB1_km$X1<0.05)
#how many genes that scale strongly predict significant survival
scaling_survival1=subset(LMNB1_km, LMNB1_km$X1<0.05 & LMNB1_km$V5<1)
scaling_survival2=subset(LMNB1_km, LMNB1_km$X1<0.05 & LMNB1_km$V5>1)

plot(volcano[,3],volcano[,1],xlab=parse(text=paste("R^2")),ylab="Scaling exponent",pch=19)

#scaling_survival1=subset(LMNB1_km, LMNB1_km$X1<0.05)
#what isthe percentage with respect to all the scaling genes
per=dim(scaling_survival1)[1]/dim(LMNB1_km)[1]
per1=dim(scaling_survival2)[1]/dim(LMNB1_km)[1]
plot(-log2(as.double(kmsig2[,1])),log2(as.double(kmsig2[,5])),pch = 19,cex=0.3,col="grey",xlab='-log2(p-value)',ylab='log2(median survival factor change)',title(c(signif(per*100,2)," % show poor survival",per1*100," % show pro survival")))


points(-log2(LMNB1_km[,1]),log2(LMNB1_km[,5]),pch = 19,cex=0.5,col="red")
text(-log2(scaling_survival1[,1]),log2(scaling_survival1[,5]),scaling_survival1[,4],pos=3)
points(c(-log2(0.05),-log2(0.05),-log2(0.05)),-1:1,pch = 3,cex=1)

#Individual survival plots
gene_of_interest="CD47"
iden=which(df.map$V1==gene_of_interest)-1
Trk[,1]=as.double(t(Trk1[iden,]))
Trk[,4]=NULL
Trk[,4] <- Trk[,1]> median(Trk[,1])
fit <- survfit(Surv(Trk[,2],Trk[,3]) ~ Trk[,4],
               data = Trk)
ggsurvplot(fit, data = Trk, risk.table = FALSE,pval = TRUE,xlab = "Survival time in years",font.x=15,font.y=15,font.tickslab=15,legend.labs=c("low expressers","high expressers"))

gene2="TYRP1"
iden1=which(df.map$V1==gene2)-1

Trk[,1]=as.double(t(Trk1[iden,]))
Trk[,4]=as.double(t(Trk1[iden1,]))
Trk[,5]=Trk[,1]> median(Trk[,1]) & Trk[,4]> median(Trk[,4])

trk2=subset(Trk, (Trk[,1]> median(Trk[,1]) & Trk[,4]> median(Trk[,4])) | (Trk[,1]< median(Trk[,1]) & Trk[,4]< median(Trk[,4]))  )

fit <- survfit(Surv(trk2[,2],trk2[,3]) ~ trk2[,5],
               data = trk2)
ggsurvplot(fit, data = trk2, risk.table = FALSE,pval = TRUE,xlab = "Survival time in years",font.x=15,font.y=15,font.tickslab=15,legend.labs=c("low expressers","high expressers"))



library(ggrepel)
#set.seed(42)
#---------------Plotting Gene Ontology enrichment----------------------
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018","KEGG_2019_Human","WikiPathways_2019_Human")
enriched <- enrichr(df.strong_scaling$V5, dbs)
plot(df.strong_scaling[,3],df.strong_scaling[,1],pch=19,xlab=parse(text="R^2"),ylab="Exponent")
text(df.strong_scaling[,3],df.strong_scaling[,1],df.strong_scaling[,5],pos=3)

an1=enriched[["GO_Biological_Process_2018"]]$Genes[enriched[["GO_Biological_Process_2018"]]$Term=="extracellular matrix organization (GO:0030198)"]
an2=enriched[["GO_Biological_Process_2018"]]$Genes[enriched[["GO_Biological_Process_2018"]]$Term=="regulation of angiogenesis (GO:0045765)"]
an3=enriched[["GO_Biological_Process_2018"]]$Genes[enriched[["GO_Biological_Process_2018"]]$Term=="regulation of cell migration (GO:0030334)"]

an1s=unlist(strsplit(an1, ";"))
an1p=subset(df.strong_scaling,df.strong_scaling$V5 %in% an1s)
points(an1p[,3],an1p[,1],pch=19,col="red")
text(an1p[,3],an1p[,1],an1p[,5],pos=3)
an2s=unlist(strsplit(an2, ";"))
an2p=subset(df.strong_scaling,df.strong_scaling$V5 %in% an2s)
text(an2p[,3],an2p[,1],an2p[,5],pos=3)
points(an2p[,3],an2p[,1],pch=19,col="blue")
an3s=unlist(strsplit(an3, ";"))
an3p=subset(df.strong_scaling,df.strong_scaling$V5 %in% an3s)
points(an3p[,3],an3p[,1],pch=19,col="green")
#text(an3p[,3],an3p[,1],an3p[,5],pos=3)
#----alternate repel text label plotting 
ggplot(df.strong_scaling,aes(V3,V1,label=V5))+geom_point(color = "blue")+geom_label_repel()
#enriched[["KEGG_2019_Human"]]
ggplot(kmsig2,aes(-log(),V1,label=V4))+geom_point(color = "blue")+geom_label_repel()
