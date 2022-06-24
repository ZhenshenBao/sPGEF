library(igraph)
library(org.Hs.eg.db)
library(clusterProfiler)

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


dsubjects <- c("Subject 01","Subject 05","Subject 06","Subject 07","Subject 10","Subject 12","Subject 15")
nsubjects <- c("Subject 02","Subject 03","Subject 04","Subject 09","Subject 11","Subject 14","Subject 16")
subjects <- c(dsubjects,nsubjects)


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

top1 <- function(x){
  #n <- length(IFadjgl)
  x <- x[order(-x)]
  x <- x[1:100]
  return(x)
}

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
        cossimilar <-  1 - (sum(samejcos*same1cos)/(sqrt(sum(samejcos^2))*sqrt(sum(same1cos^2))))######Cosine distance#########
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
