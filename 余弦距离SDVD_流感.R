######余弦距离#########
library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)

#######流感数据处理######
anot <- GPL9188
anot <- as.data.frame(anot)
colnames(anot) <- c("ID","ENTREZID")
anot <- as.data.frame(anot)
samples <- GSE30550_series_matrix[1,]
exdata <- GSE30550_series_matrix[-1,]
colnames(exdata) <- exdata[1,]
exdata <- exdata[-1,]
rownames(exdata) <- exdata$ID_REF
exdata <- exdata[,-1]
esetm = exdata[rownames(exdata) %in% anot$ID, ]
esetm$ENTREZID <- as.character(anot$ENTREZID[match(rownames(esetm), anot$ID)])
esetm<-esetm[!is.na(esetm$ENTREZID),]
esetm<-esetm[!duplicated(esetm$ENTREZID),]
rownames(esetm) <- esetm$ENTREZID
esetm <- esetm[,1:268]
samples <- samples[-1]
circsp<-strsplit(as.character(samples),",")
samples<-do.call(rbind,circsp)
tts <- cbind(samples,colnames(esetm))

genesy <- bitr(rownames(esetm), fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
esetm <- esetm[genesy[,1],]
rownames(esetm) <- genesy[,2]

interg <- intersect(proteininfo$preferred_name,rownames(esetm))
esetm <- esetm[interg,]

#######PPI网络数据处理########
IN <- NULL
INlist <- list()
for(i in 1:nrow(esetm)){
  egene <- rownames(esetm)[i]
  url <- paste("https://string-db.org/api/tsv/interaction_partners?identifiers=",egene,"&species=Homo%20sapiens",sep = "")
  webDf <- read.table(url, header=T)
  IN <- rbind(IN,webDf)
  cat("\n");
  cat(paste("Done  ",i,sep=""))
}
IN1 <- IN[which(IN$score>=0.8),]
INglist <- IN1[,3:4]
INg <- graph_from_edgelist(as.matrix(INglist),directed = F)
IFadjgl <- as_adj_list(INg)

######获取有症状样本与无症状样本############
dsubjects <- c("Subject 01","Subject 05","Subject 06","Subject 07","Subject 10","Subject 12","Subject 15")
nsubjects <- c("Subject 02","Subject 03","Subject 04","Subject 09","Subject 11","Subject 14","Subject 16")
subjects <- c(dsubjects,nsubjects)

######Z-scorer标准化流感表达数据############
esetm <- apply(esetm,2,as.numeric)
esetm1 <- 2^(esetm)
esetm_zscore <- scale(esetm1)
esetm_bm <- esetm_zscore
esetm_bm[which(esetm_bm >= 1)] <- 2
esetm_bm[which(esetm_bm <= -1)] <- 2
esetm_bm[which(esetm_bm != 2)] <- 1
rownames(esetm) <- interg
rownames(esetm_bm) <- rownames(esetm)

for(i in 1:length(dsubjects)){
  ttd <- tts[which(tts[,1] == dsubjects[[i]]),]
  case_LG1 <- esetm1[,ttd[,3]]
  write.csv(case_LG1,paste("case_LG",i,".csv",sep=""))
}


######临界点基因选取函数##########
degene <- function(x){
  gnames <- rownames(score)
  names(x) <- gnames
  x <- x[order(-x)]
  x <- x[1:380]
  degenes <- names(x)
  return(degenes)
}

######选取前100基因########
top1 <- function(x){
  #n <- length(IFadjgl)
  x <- x[order(-x)]
  x <- x[1:100]
  return(x)
}

######开始计算######
ptm <- proc.time()
degreemax <- 10
cosscores <- NULL
for(ii in 1:7){
  tt1 <- tts[which(tts[,1] == subjects[[ii]]),]
  samsone <- tt1[,3]
  same <- esetm_bm[,samsone]
  same <- apply(same,2,as.numeric)
  sames <- same[,1:4]
  sames[which(sames==1)] <- 0
  sames[which(sames==2)] <- 1
  same1 <- apply(sames,1,sum)
  same1[which(same1 <= 3)] <- 1 
  same1[which(same1 > 3)] <- 2
  names(same1) <- rownames(esetm)
  cos <- NULL
  allsd <- NULL
  allcosgenes <- NULL
  for(i in 1:length(IFadjgl)){
    cosgene <- names(IFadjgl[[i]])
    cosgene <- cosgene[!duplicated(cosgene)]
    same1cos <- same1[cosgene]
    same1cos <- same1cos[!is.na(same1cos)]
    a <- intersect(names(IFadjgl)[i],rownames(esetm))
    same1fc <- abs(same1[a]-same2[a])
    if(length(a) != 0){
      same1cos <- c(same1cos,same1[a])
      if(length(same1cos) >= degreemax){
        allcosgenes <- c(allcosgenes,names(IFadjgl)[i])
      }
    }
    cossimilars <- NULL
    gsds <- NULL
    for(j in 1:16){
      samej <- as.numeric(same[,j])
      names(samej) <- rownames(esetm)
      samejcos <- samej[cosgene]
      samejcos <- samejcos[!is.na(samejcos)]
      samejfc <- samej[a]
      if(length(a) != 0){
        samejcos <- c(samejcos,samej[a])
      }
      if(length(a)!=0 && length(same1cos) >= degreemax ){
        cossimilar <-  1 - (sum(samejcos*same1cos)/(sqrt(sum(samejcos^2))*sqrt(sum(same1cos^2))))
        cossimilar <- (samejfc-1)*cossimilar
        cossimilars <- c(cossimilars,cossimilar)
      }
    }
    cos <- rbind(cos,cossimilars)
  }
  cos[is.na(cos)] <- 0
  cos[which(cos<=0)] <- 0
  rownames(cos) <- allcosgenes
  score <- cos
  score1 <- apply(score,2,top1)
  cosscore <- apply(score1,2,sum)
  #tpdeg <- apply(score,2,degene)
  cosscores <- rbind(cosscores,cosscore)
  #write.csv(score,paste("mtpdeg",ii,".csv",sep=""))
  cat("\n");
  cat(paste("Done  ",ii,sep=""))
}
proc.time()-ptm
rownames(cosscores) <- subjects
colnames(cosscores) <- tt1[,2]
write.csv(cosscores,"mdscores.csv")


####LUAD#####
TCGA.LUAD.htseq_counts <- read.delim("C:/Users/pc/Desktop/TCGA-LUAD.htseq_counts.tsv", stringsAsFactors = F)
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
list = AnnotationDbi::select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
ID <- TCGA.LUAD.htseq_counts$Ensembl_ID
IDs <- strsplit(ID,"\\.")
IDs <- do.call(rbind,IDs)
ID <- IDs[,1]
ID_list=list[match(ID,list[,"ENSEMBL"]),]
expd <- TCGA.LUAD.htseq_counts
expd$Ensembl_ID <- ID_list$SYMBOL

stages <- TCGA_LUAD_stage
stages <- stages[which(stages[,2] != 'not reported'),]
samples <- colnames(expd)
samples <- samples[-1]
samples <- gsub("\\.", "-", samples)
colnames(expd) <- c("SYMBOL",samples)
expd <- expd[!duplicated(expd$SYMBOL),]
expd <- expd[!is.na(expd$SYMBOL),]
rownames(expd) <- expd$SYMBOL
stages <- stages[match(samples,stages$submitter_id.samples),]
stages <- stages[!is.na(stages[,2]),]
expd <- expd[,stages$submitter_id.samples]
sams <- stages$submitter_id.samples
classes <- strsplit(sams,"-")
classes <- do.call(rbind,classes)
classes <- classes[,4]
stages$"classes" <- classes

refstages <- stages[which(stages$classes == "11A" | stages$classes == "11B"),]
diseasestages <- stages[which(stages$classes != "11A" & stages$classes != "11B"),]
hours <- names(table(diseasestages$tumor_stage.diagnoses))

expd <- scale(expd)
expd_bm <- expd
expd_bm[which(expd_bm >= 1)] <- 2
expd_bm[which(expd_bm <= -1)] <- 2
expd_bm[which(expd_bm != 2)] <- 1

refsams <- refstages$submitter_id.samples
erefsams <-  expd_bm[,refsams]
erefsams[which(erefsams==1)] <- 0
erefsams[which(erefsams==2)] <- 1
#same1 <- same[,1]
erefsams1 <- apply(erefsams,1,sum)
erefsams1[which(erefsams1 <= 57)] <- 1 
erefsams1[which(erefsams1 > 57)] <- 2

des <- NULL
for(j in 1:length(hours)){
  hourssams <- diseasestages[which(diseasestages$tumor_stage.diagnoses == hours[j]),]$submitter_id.samples
  expdde <- expd_bm[,hourssams]
  if(length(hourssams)!=1){
    expdde[which(expdde==1)] <- 0
    expdde[which(expdde==2)] <- 1
    expdde <- apply(expdde,1,sum)
    expdde[which(expdde >= 1)] <- 2 
    expdde[which(expdde == 0)] <- 1
  }
  des <- cbind(des,expdde)
}

###########获取肺腺癌PPI#######
interg <- intersect(proteininfo$preferred_name,rownames(des))
des <- des[interg,]
IN2 <- NULL
INlist1 <- list()
for(i in 1:nrow(des)){
  egene1 <- rownames(des)[i]
  url1 <- paste("https://string-db.org/api/tsv/interaction_partners?identifiers=",egene1,"&species=Homo%20sapiens",sep = "")
  webDf1 <- read.table(url1, header=T)
  IN2 <- rbind(IN2,webDf1)
  cat("\n");
  cat(paste("Done  ",i,sep=""))
}
IN3 <- IN2[which(IN2$score>=0.8),]
INglist1 <- IN3[,3:4]
INg1 <- graph_from_edgelist(as.matrix(INglist1),directed = F)
IFadjgl <- as_adj_list(INg1)
erefsams1 <- erefsams1[interg]

ptm <- proc.time()
degreemax <- 10
cos <- NULL
allsd <- NULL
allcosgenes <- NULL
for(i in 1:length(IFadjgl)){
  cosgene <- names(IFadjgl[[i]])
  cosgene <- cosgene[!duplicated(cosgene)]
  same1cos <- erefsams1[cosgene]
  same1cos <- same1cos[!is.na(same1cos)]
  
  
  a <- intersect(names(IFadjgl)[i],rownames(des))
  if(length(a) != 0){
    same1cos <- c(same1cos,erefsams1[a])
    if(length(same1cos) >= degreemax){
      allcosgenes <- c(allcosgenes,names(IFadjgl)[i])
    }
  }
  cossimilars <- NULL
  gsds <- NULL
  
  for(j in 1:9){
    samej <- as.numeric(des[,j])
    names(samej) <- rownames(des)
    samejcos <- samej[cosgene]
    samejcos <- samejcos[!is.na(samejcos)]
    samejfc <- samej[a]
    
    if(length(a) != 0){
      samejcos <- c(samejcos,samej[a])
    }
    if(length(a)!=0 && length(same1cos) >= degreemax ){
      cossimilar <- 1 - (sum(samejcos*same1cos)/(sqrt(sum(samejcos^2))*sqrt(sum(same1cos^2))))
      cossimilar <- (samejfc-1)*cossimilar
      cossimilars <- c(cossimilars,cossimilar)
    }
  }
  cos <- rbind(cos,cossimilars)
}
cos[is.na(cos)] <- 0
cos[which(cos<=0)] <- 0
rownames(cos) <- allcosgenes
score <- cos
score1 <- apply(score,2,top1)
cosscore <- apply(score1,2,sum)
proc.time()-ptm

write.csv(cosscore,"cdSDVD_FA.csv")
write.csv(score,"cdSDVDall_FA.csv")


########富集分析1##########
geneEZ <- bitr(SDVDLGgene, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene=geneEZ[,2],OrgDb = org.Hs.eg.db,ont ="BP",pAdjustMethod = "fdr",pvalueCutoff = 0.2,readable = TRUE)
ego2 <- data.frame(ego)
ekegg <- enrichKEGG(gene=geneEZ[,2],pAdjustMethod = "fdr",pvalueCutoff = 0.2)
ekegg2 <- data.frame(ekegg)

SLELGgene <- as.vector(t(SLELGgene))
geneEZ1 <- bitr(SLELGgene, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
ego1 <- enrichGO(gene=geneEZ1[,2],OrgDb = org.Hs.eg.db,ont ="BP",pAdjustMethod = "fdr",pvalueCutoff = 0.2,readable = TRUE)
ego3 <- data.frame(ego1)
ekegg1 <- enrichKEGG(gene=geneEZ1[,2],pAdjustMethod = "fdr",pvalueCutoff = 0.2)
ekegg3 <- data.frame(ekegg1)
write.csv(ekegg2,"SDVD_lg_kegg.csv")
write.csv(ekegg3,"SLE_lg_kegg.csv")

######富集分析2#######
SDVDFAgene <- as.vector(t(SDVD_FA_gene))
geneEZ_FA <- bitr(SDVDFAgene, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
ego_FA <- enrichGO(gene=geneEZ_FA[,2],OrgDb = org.Hs.eg.db,ont ="BP",pAdjustMethod = "fdr",pvalueCutoff = 0.2,readable = TRUE)
ego2_FA <- data.frame(ego_FA)
ekegg_FA <- enrichKEGG(gene=geneEZ_FA[,2],pAdjustMethod = "fdr",pvalueCutoff = 0.2)
ekegg2_FA <- data.frame(ekegg_FA)


SLEFAgene <- as.vector(t(SLE_FA_gene))
geneEZ1_FA <- bitr(SLEFAgene, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
ego1_FA <- enrichGO(gene=geneEZ1_FA[,2],OrgDb = org.Hs.eg.db,ont ="BP",pAdjustMethod = "fdr",pvalueCutoff = 0.2,readable = TRUE)
ego3_FA <- data.frame(ego1_FA)
ekegg1_FA <- enrichKEGG(gene=geneEZ1_FA[,2],pAdjustMethod = "fdr",pvalueCutoff = 1)
ekegg3_FA <- ekegg1_FA@result
write.csv(ego2_FA,"SDVD_FA_go.csv")
write.csv(ego3_FA,"SLE_FA_go.csv")
write.csv(ekegg2_FA,"SDVD_FA_kegg.csv")
write.csv(ekegg3_FA,"SLE_FA_kegg.csv")


venn.diagram(list(WDSP = SDVDLGgene,Pfam= SLELGgene),height = 450, width = 450,
             resolution = 600, imagetype = "svg", alpha=c(0.5,0.5),
             fill=c("red","blue"), cat.fontface=4,fontfamily=3,
             main.cex = 2, main.fontface = 2, main.fontfamily = 3,
             filename = "VennDiagram.svg")


sJSDLGgene <- as.vector(t(sJSD_LG))
#geneEZ1_LG <- bitr(sJSDLGgene, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
ego1_LG <- enrichGO(gene=sJSDLGgene,OrgDb = org.Hs.eg.db,ont ="BP",pAdjustMethod = "fdr",pvalueCutoff = 0.2,readable = TRUE)
ego3_LG <- data.frame(ego1_LG)
ekegg1_LG <- enrichKEGG(gene=sJSDLGgene,pAdjustMethod = "fdr",pvalueCutoff = 1)
ekegg3_LG <- ekegg1_LG@result
write.csv(ego3_LG,"sJSD_LG_go.csv")
write.csv(ekegg3_LG,"sJSD_LG_kegg.csv")

sJSDFAgene <- as.vector(t(sJSD_FA_genes))
sJSDFAgene <- strsplit(sJSDFAgene,"'")
sJSDFAgene <- do.call(rbind,sJSDFAgene)
sJSDFAgene <- as.vector(t(sJSDFAgene[,2]))
geneEZ1_FA <- bitr(sJSDFAgene, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
ego1_LG <- enrichGO(gene=geneEZ1_FA[,2],OrgDb = org.Hs.eg.db,ont ="BP",pAdjustMethod = "fdr",pvalueCutoff = 0.2,readable = TRUE)
ego3_LG <- data.frame(ego1_LG)
ekegg1_LG <- enrichKEGG(gene=geneEZ1_FA[,2],pAdjustMethod = "fdr",pvalueCutoff = 1)
ekegg3_LG <- ekegg1_LG@result
write.csv(ego3_LG,"sJSD_FA_go.csv")
write.csv(ekegg3_LG,"sJSD_FA_kegg.csv")
