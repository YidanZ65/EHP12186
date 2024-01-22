####手掌####

foreheadr <- foreheadr %>%
  mutate(NDVI250meg = Median(NDVI250mea),
         NDVI1000meg = Median(NDVI1000mea))



pcoa <- cmdscale(br_oralr, k=3, eig=T) 
pcoa_points <- as.data.frame(pcoa$points)
sum_eig <- sum(pcoa$eig)
eig_percent <- round(pcoa$eig/sum_eig*100,1)
colnames(pcoa_points) <- paste0("PCoA", 1:3)
pcoa_result <- cbind(pcoa_points, oralr)

####PERMANOVA+PCOAfig####

####forehead

###ASV
permanova_fore <- data.frame()
permanova_oral_uf <- data.frame()
greenness <- c("NDVI500meg","EVI500meg")
for (x in greenness) {
  permanova <- adonis2(br_oralr ~ get(x) + factor(sex)+age+factor(collection_season)+factor(Caucasian)+factor(pop_binary),
                       data = oralr, permutations=999)
 # permanova_fore <- rbind(permanova_fore, permanova)
  adonis <- paste0("adonis R2: ",round(permanova$R2[1],3), "; P-value: ",permanova$`Pr(>F)`[1])
  p <- ggplot(pcoa_result, aes_string(x="PCoA1", y="PCoA2", color=x, group=x, shape=x)) +
    labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
         title = adonis) +
    geom_point(size=1
    ) + stat_ellipse(level=0.95, aes_string(fill = x ),geom = "polygon", alpha = 1/10) +
    theme_classic()
  print(p)
  ggsave(p, file=paste0("o","_",x,"_","a",".pdf"), width=4, height=3)
}
write.csv(permanova_fore, file = "permanova_fore(asv+genus).csv")
write.csv(permanova_fore, file = "permanova_fore_250+1000.csv")

####hand
permanova_hand <- data.frame()
greenness <- c("NDVI250meg","EVI250meg","NDVI1000meg","EVI1000meg")
for (x in greenness) {
  permanova <- adonis2(br_handr ~ get(x) + factor(sex)+age+factor(collection_season)+factor(Caucasian)+factor(pop_binary),
                       data = handr, permutations=999)
  permanova_hand <- rbind(permanova_hand, permanova)
  adonis <- paste0("adonis R2: ",round(permanova$R2[1],3), "; P-value: ",permanova$`Pr(>F)`[1])
  p <- ggplot(pcoa_result, aes_string(x="PCoA1", y="PCoA2", color=x, group=x, shape=x)) +
    labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
         title = adonis) +
    geom_point(size=1
    ) + stat_ellipse(level=0.95, aes_string(fill = x ),geom = "polygon", alpha = 1/10) +
    theme_classic()
  print(p)
  ggsave(p, file=paste0("h","_",x,"_","a2",".pdf"), width=4, height=3)
}
write.csv(permanova_hand, file = "permanova_hand(ASV+genus).csv")

####oral
permanova_oral <- data.frame()
greenness <- c("NDVI250meg","EVI250meg","NDVI1000meg","EVI1000meg")
for (x in greenness) {
  permanova <- adonis2(br_oralr ~ get(x) + factor(sex)+age+factor(collection_season)+factor(Caucasian)+factor(pop_binary),
                       data = oralr, permutations=999)
  permanova_oral <- rbind(permanova_oral, permanova)
  adonis <- paste0("adonis R2: ",round(permanova$R2[1],3), "; P-value: ",permanova$`Pr(>F)`[1])
  p <- ggplot(pcoa_result, aes_string(x="PCoA1", y="PCoA2", color=x, group=x,shape=x)) +
    labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
         title = adonis) +
    geom_point(size=1, 
    ) + stat_ellipse(level=0.95, aes_string(fill = x ),geom = "polygon", alpha = 1/10) +
    theme_classic()
  print(p)
  ggsave(p, file=paste0("o","_",x,"_","a2",".pdf"), width=4, height=3)
}
write.csv(permanova_oral, file = "permanova_oral(genus+ASV).csv")
write.csv(permanova_oral, file = "permanova_oral_250+1000.csv")


####gut####
permanova_gut <- read.csv(file = "permanova_gut(g+a).csv")
permanova_gut_uf <- read.csv(file="uf_permanova_gut.csv")
pcoa_result_gut_g <- read.csv(file = "pcoa_result_gut_g.csv")
pcoa_result_gut_a <- read.csv(file = "pcoa_result_gut_asv.csv")
pcoa_result_gut_uf <- read.csv(file = "pcoa_result_gut_uf.csv")
greenness <- c("NDVI500meg","EVI500meg","NDVI250meg","EVI250meg","NDVI1000meg","EVI1000meg")
for (x in c(1:2)) {
  adonis <- paste0("adonis R2: ",round(permanova_gut_uf$R2[x],4), "; P-value: ",permanova_gut_uf$`Pr..F.`[x])
  p <- ggplot(pcoa_result_gut_uf, aes_string(x="PCoA1", y="PCoA2", color=permanova_gut_uf$X[x], 
                                             group=permanova_gut_uf$X[x], shape=permanova_gut_uf$X[x])) +
    labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
         title = adonis) +
    geom_point(size=1, 
    ) + stat_ellipse(level=0.95, aes_string(fill = permanova_gut_uf$X[x] ),geom = "polygon", alpha = 1/10) +
    theme_classic()
  print(p)
  ggsave(p, file=paste0("g","_",permanova_gut_uf$X[x],"_","uf",".pdf"), width=4, height=3)
}




###unifrac####
library(vegan)
library(ape)
uf_oral <- read.table(file="uf-distance-matrix-oral899.tsv",header = T,
                      sep="\t", dec=".", quote = '',comment.char = '',fill = T)
rname <- paste("X",uf_oral[,1],sep = "")
rownames(uf_oral) <- rname
uf_oral <- uf_oral[,-1]


pcoa <- cmdscale(uf_oral, k=3, eig=T) 
pcoa_points <- as.data.frame(pcoa$points)
sum_eig <- sum(pcoa$eig)
eig_percent <- round(pcoa$eig/sum_eig*100,1)
colnames(pcoa_points) <- paste0("PCoA", 1:3)
pcoa_result <- cbind(pcoa_points, oralr)


uf_permanova_oral <- data.frame()
greenness <- c("NDVI500meg","EVI500meg","NDVI250meg","EVI250meg","NDVI1000meg","EVI1000meg")
for (x in greenness) {
  permanova <- adonis2(uf_oral ~ get(x) + factor(sex)+age+factor(collection_season)+factor(Caucasian)+factor(pop_binary),
                       data = oralr, permutations=999)
  uf_permanova_oral <- rbind(uf_permanova_oral, permanova)
  adonis <- paste0("adonis R2: ",round(permanova$R2[1],3), "; P-value: ",permanova$`Pr(>F)`[1])
  p <- ggplot(pcoa_result, aes_string(x="PCoA1", y="PCoA2", color=x, group=x, shape=x)) +
    labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
         title = adonis) +
    geom_point(size=1
    ) + stat_ellipse(level=0.95, aes_string(fill = x ),geom = "polygon", alpha = 1/10) +
    theme_classic()
  print(p)
  ggsave(p, file=paste0("o","_",x,"_","uf",".pdf"), width=4, height=3)
}

write.csv(uf_permanova_hand, file = "uf_permanova_hand.csv")



####β-diversity sensitivity analysis####

# hand, forehead, gut, oral,forehead
permanova_sense_asv <- data.frame()
permanova_diet_g <- data.frame()
permanova_diet_asv <- data.frame()
sensitive1 <- c("antibiotic_month2","ibd_cat","diabetes_cat")
sensitive2 <- c("pet_ownership2","bmi")
sensitive3 <- c("diet_cat")

oralr <- oralr %>%
  mutate(antibiotic_month2 = case_when((antibiotic_month == "with in a month")~"Yes",
                                       (antibiotic_month == "Missing" |
                                          antibiotic_month == "over a month")~"No"))
#hand, forehead,oral(NDVI5-EVI5)
for (x in sensitive1) {
  m <- filter(oralr,  get(x) != "Yes")
  br <- br_g_oralr[m$X.SampleID,m$X.SampleID]
  a<- adonis2(br ~ EVI500meg + factor(sex)+age+factor(collection_season)+factor(pop_binary)
              +factor(Caucasian),data = m, permutations=999)
  b <- a[which(rownames(a)=="EVI500meg"),c(which(colnames(a)=="R2"),which(colnames(a)=="F"),which(colnames(a)=="Pr(>F)"))]
  b$sen <- x
  permanova_diet_g <- rbind(permanova_diet_g,b)
}

for (x in sensitive3) {
  m <- filter(oralr,  get(x) != "Not provided")
  br <- br_oralr[m$X.SampleID,m$X.SampleID]
  a<- adonis2(br ~ EVI500meg + factor(sex)+age+factor(collection_season)+factor(pop_binary)
              +factor(Caucasian)+get(x),data = m, permutations=999)
  b <- a[which(rownames(a)=="EVI500meg"),c(which(colnames(a)=="R2"),which(colnames(a)=="F"),which(colnames(a)=="Pr(>F)"))]
  b$sen <- x
  permanova_diet_asv <- rbind(permanova_diet_asv,b)
}

write.csv(permanova_sense_asv, file = "permanova_asv(F-value).csv")
write.csv(permanova_sense_g, file = "permanova_sense_g(F-value).csv")

write.csv(permanova_diet_asv, file = "permanova_diet_asv(F-value).csv")

mod <- betadisper(br_g_handr, handr$EVI500meg)
permutest(mod)



handr_gf <- t(handr_gf)
foreheadr_gf <- t(foreheadr_gf)
oralr_gf <- t(oralr_gf)


taxo_hand2 <- data.frame()
taxo_hand2 <- rbind(res, taxo_hand2)

taxo_hand3 <- data.frame()
taxo_hand3 <- rbind(res, taxo_hand3)


taxo_forehead2 <- data.frame()
taxo_forehead2 <- rbind(res, taxo_forehead2)

taxo_forehead3 <- data.frame()
taxo_forehead3 <- rbind(res, taxo_forehead3)

taxo_oral2 <- data.frame()
taxo_oral2 <- rbind(res, taxo_oral2)
taxo_oral3 <- data.frame()
taxo_oral3 <- rbind(res, taxo_oral3)

taxo_gut1000 <- data.frame()
taxo_gut1000 <- rbind(res, taxo_gut1000)

library(DESeq2)
library(BiocParallel)
register(SnowParam(4))

dds <- DESeqDataSetFromMatrix(countData = handr_gf, 
                              colData = handr[which(handr$antibiotic_month != "Yes")], 
                              design = ~ factor(sex)+factor(age_cat)+factor(Caucasian)+factor(collection_season)+factor(pop_binary)
                              +NDVI500meg)
dds <- DESeq(dds,sfType = "poscounts", parallel = T)
res <- results(dds, contrast=c("NDVI500meg",'high','low'))
summary(res,0.05)
res <- as.data.frame(res)
res <- dplyr::filter(res, !is.na(res$padj))
names(res) <- paste("EVI1000m","_",names(res),sep = "")
res$Row.names <- rownames(res)

taxo_gut1000 <- merge(taxo_gut1000,res, by="Row.names",
                   all.x=TRUE,all.y=TRUE)

hp_gut1000 <- taxo_gut1000[,c("Row.names",
                          "NDVI1000m_baseMean","NDVI1000m_log2FoldChange","NDVI1000m_padj",
                          "EVI1000m_baseMean","EVI1000m_log2FoldChange","EVI1000m_padj")]

hp_gut1000 <- subset(hp_gut1000,NDVI1000m_padj<0.05 |EVI1000m_padj<0.05)
write.csv(hp_oral2, file = "hp_oralNDVI+EVI500原始.csv")
write.csv(hp_oral3, file = "hp_oralNDVI+EVI2501000原始.csv")

hp_oral2 <- hp_oral2[,c(6:12)] %>% filter(., g !="g__"& g !="__") 


hp_oral2$g <- substr(hp_oral2$g,4,nchar(hp_oral2$g))
hp_oral2$NDVI500m_log2FoldChange[is.na(hp_oral2$NDVI500m_log2FoldChange)] <- 0
hp_oral2$g <- factor(hp_oral2$g, 
                    levels=hp_oral2$g[order(hp_oral2$NDVI500m_log2FoldChange, decreasing=FALSE)])



library(reshape2)
hp_forehead10001<-melt(hp_forehead1000,
                 id.vars = c("g"),#需要保留不参与聚合的变量,
                 measure.vars = c('NDVI1000m_log2FoldChange',"EVI1000m_log2FoldChange"),#用于聚合的变量,
                 variable.name='index',
                 value.name= c('lfc'))
hp_forehead10002<-melt(hp_forehead1000,
                 id.vars = c("g"),#需要保留不参与聚合的变量,
                 measure.vars = c('NDVI1000m_padj',"EVI1000m_padj"),#用于聚合的变量,
                 variable.name='index_p',
                 value.name= c('padj'))
hp_forehead1000l <- cbind(hp_forehead10001,hp_forehead10002) %>% .[,-c(4:5)]
hp_forehead1000l$padj[is.na(hp_forehead1000l$padj)] <- 1#p缺失记为1
hp_forehead1000l$lfc[is.na(hp_forehead1000l$lfc)] <- 0
hp_forehead1000l$display <- 
  ifelse(hp_forehead1000l$padj <= 0.05, "*", "")#标注P

hp_forehead1000l$index <- as.character(hp_forehead1000l$index)
hp_forehead1000l$index <- substr(hp_forehead1000l$index,1,nchar(hp_forehead1000l$index)-15)
hp_forehead1000l$index <- factor(hp_forehead1000l$index, 
                           levels=c("NDVI1000m","EVI1000m"), ordered=T)




pp <- ggplot(hp_forehead1000l,aes(index,g,fill=lfc)) + 
  geom_tile()+
  theme_minimal()+
  scale_fill_gradient2(low = "#0099ff",high = "#ce472e",mid = "white",midpoint = 0,
                       limit = c(-6,6), space = "Lab",
                       name="LFC")+
  geom_text(aes(label=display),size=6,color="black",nudge_y = -0.15)+
  scale_y_discrete(position="left")+xlab(NULL) + ylab(NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust=0.45,vjust=0.45))+
  theme(axis.text.x=element_text(family = "Times",
                                 face = "plain",colour = "black",size=7))+
  theme(axis.text.y=element_text(family = "Times",
                                 face = "plain",colour = "black",size=7))+
  theme(legend.text=element_text(face="plain",family = "Times",
                                 colour = "black",size = 7))+
  labs(fill = "")

pp <- pp + coord_polar(theta="y",start = 2* pi / 3)

pp


####Table1####
library(tableone)


myVars <- c("age", "sex", "country_cat","bmi","Caucasian",
            "pet_ownership","antibiotic_month","collection_season","fed_infant_cat",
            "level_of_education","probiotic_frequency", "ibd_cat","diabetes_cat",
            "shannon_entropy","pielou_evenness","observed_features")

idx <- which(names(gut) %in% catVars)
for(i in idx ){
  gut[,i]  <-  as.factor(gut[,i])
}

catVars <- c("sex", "country_cat",
             "pet_ownership","Caucasian",
             "antibiotic_history_cat","collection_season","fed_infant_cat",
             "level_of_education","probiotic_frequency","ibd_cat","diabetes_cat")
a <- CreateTableOne(vars = myVars, data = gutr,factorVars = catVars)
table1 <- print(a, quote = FALSE, noSpaces = TRUE, printToggle = FALSE,showAllLevels=TRUE)
write.csv(table1, file = "gutr_table1.csv")

race <- c("race")
a <- CreateTableOne(vars = race, data = foreheadr)
table1 <- print(a, quote = FALSE, noSpaces = TRUE, printToggle = FALSE,showAllLevels=TRUE)
write.csv(table1, file = "foreheadr_race.csv")


####DESeq2+环状热图####
library(ComplexHeatmap)
library(circlize)

tax_hand500 <- data.frame()
tax_hand500 <- rbind(res, tax_hand500)

tax_hand250 <- data.frame()
tax_hand250 <- rbind(res, tax_hand250)

tax_hand1000 <- data.frame()
tax_hand1000 <- rbind(res, tax_hand1000)

tax_forehead500 <- data.frame()
tax_forehead500 <- rbind(res, tax_forehead500)

taxo_forehead5002 <- data.frame()
taxo_forehead5002 <- rbind(res, taxo_forehead5002)

tax_forehead250 <- data.frame()
tax_forehead250 <- rbind(res, tax_forehead250)

tax_forehead1000 <- data.frame()
tax_forehead1000 <- rbind(res, tax_forehead1000)

tax_oral500 <- data.frame()
tax_oral500 <- rbind(res, tax_oral500)

tax_oral250 <- data.frame()
tax_oral250 <- rbind(res, tax_oral250)

tax_oral1000 <- data.frame()
tax_oral1000 <- rbind(res, tax_oral1000)

tax_gut500 <- data.frame()
tax_gut500 <- rbind(res, tax_gut500)

tax_gut250 <- data.frame()
tax_gut250 <- rbind(res, tax_gut250)

tax_gut1000 <- data.frame()
tax_gut1000 <- rbind(res, tax_gut1000)

library(BiocParallel)
register(SnowParam(4))

dds <- DESeqDataSetFromMatrix(countData = foreheadr_gf, 
                              colData = foreheadr, 
                              design = ~ factor(sex)+factor(age_cat)+factor(Caucasian)+factor(collection_season)+factor(pop_binary)
                              +EVI1000meg)
dds <- DESeq(dds,sfType = "poscounts")
res <- results(dds, contrast=c("EVI1000meg",'high','low'))
summary(res,0.05)
res <- as.data.frame(res) %>% dplyr::filter(., !is.na(.$padj))
names(res) <- paste("EVI1000m","_",names(res),sep = "")
res$Row.names <- rownames(res)

tax_forehead1000 <- merge(tax_forehead1000,res, by="Row.names",
                    all.x=TRUE,all.y=TRUE)

hp_forehead1000 <- tax_forehead1000[,c("Row.names",
                          "NDVI1000m_baseMean","NDVI1000m_log2FoldChange","NDVI1000m_padj",
                          "EVI1000m_baseMean","EVI1000m_log2FoldChange","EVI1000m_padj")]

hp_oral500 <- subset(hp_oral500, NDVI500m_padj<0.05 | EVI500m_padj <0.05)
hp_oral250 <- subset(hp_oral250, NDVI250m_padj<0.05 | EVI250m_padj <0.05)
hp_oral1000 <- subset(hp_oral1000, NDVI1000m_padj<0.05 | EVI1000m_padj <0.05)
write.csv(hp_hand1000, file = "hp_hand1000_0106.csv")


#整理作图用的表格
hp_orals5 <- read.csv(file="hp_orals5原始.csv")
hp_orals5 <- hp_orals5[,c(6:12)] %>% filter(., g !="g__"& g !="__") %>%
  mutate(g = substr(.$g,4,nchar(.$g)))

which(duplicated(hp_foreheads1$g))
#hp_forehead250$NDVI250m_log2FoldChange[is.na(hp_forehead250$NDVI250m_log2FoldChange)] <- 0

#hp_gut1000 <- hp_forehead1000[-which(hp_gut1000$g == "Clostridium"),] #Clostridium分类错误，去除

#hp_oral500$g <- factor(hp_oral500$g, 
#                     levels=hp_oral500$g[order(hp_oral500$NDVI500m_log2FoldChange, decreasing=FALSE)])

chp_orals5 <- hp_orals4[,c(1,3,6)]


genesname500 <- rbind(chp_hand500[,c(1:3)], chp_forehead500,chp_gut500[,c(1:3)], chp_oral500[,c(1:3)]) %>%
  distinct(., g, .keep_all = T)

genesname250 <- rbind(chp_hand250[,c(1:3)], chp_forehead250,chp_gut250[,c(1:3)], chp_oral250[,c(1:3)]) %>%
  distinct(., g, .keep_all = T)

genesname1000 <- rbind(chp_hand1000[,c(1:3)], chp_forehead1000, chp_gut1000[,c(1:3)], chp_oral1000) %>%
  distinct(., g, .keep_all = T)


genesnamehand <- rbind(chp_hands1, chp_hands2, chp_hands3, chp_hands4,chp_hands5,chp_hands6) %>%
  distinct(., g, .keep_all = T)hp_forehe


genesnamegut <- rbind(chp_guts1, chp_guts2, chp_guts3, chp_guts4, chp_guts5, chp_guts6) %>%
  distinct(., g, .keep_all = T)


genesnameforehead <- bind_rows(chp_foreheads1,chp_foreheads2,chp_foreheads3,chp_foreheads4,chp_foreheads5,chp_foreheads6 ) %>%
  distinct(., g, .keep_all = T)


genesnameoral <- bind_rows(chp_orals1,chp_orals2,chp_orals3,chp_orals4) %>%
  distinct(., g, .keep_all = T)

write.csv(genesnameoral, file = "genesnameoral_0121.csv")

genesnameoral <- read.csv(file = "genesnameoral_0121.csv")


chp_orals4 <- merge(chp_orals4, genesnameoral[,1:2], by="g", all.x=T) 

#chp_gut1000$EVI1000m_log2FoldChange[is.na(chp_gut1000$EVI1000m_log2FoldChange)] <- 0

#chp_hand500$name <- factor(chp_hand500$name, 
#                               levels=chp_hand500$name[order(chp_hand500$NDVI500m_log2FoldChange, decreasing=FALSE)])

rownames(chp_orals4) <- chp_orals4$name



####作图
chp_orals3$EVI500m_log2FoldChange[which(is.na(chp_orals3$EVI500m_log2FoldChange))] <- 0

col_fun1 = colorRamp2(c(-1, 0, 1),c("#0099ff", "white", "#ce472e"))
circos.par(gap.after = c(15))
hp <- circos.heatmap(chp_orals4[,c(2:3)], col = col_fun1,track.height= 0.2,
               rownames.side = "outside",rownames.cex = 0.4, cell.border = "black",
               cell.lwd = 0.1)
#circos.text(3, 1, "inside", facing = "inside", cex = 0.8)
lgd = Legend(title = "lfc", col_fun = col_fun1)
grid.draw(lgd)
circos.clear()
#输出4*4


circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = list("EVI500m","NDVI500m")
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(0.1, "mm"), #x坐标
                1:n ,#y坐标
                cn, #标签
                cex = 0.5, adj = c(0, 1))
  }
}, bg.border = NA)
circos.clear()



#合并P值
chp_hands6 <- merge(chp_hands6, hp_hands6[,c(1,4,7)], by = "g", all.x=T)
chp_guts6 <- merge(chp_guts6, hp_guts6[,c(1,4,7)], by = "g", all.x=T)
chp_foreheads6 <- merge(chp_foreheads6, hp_foreheads6[,c(1,4,7)], by = "g", all.x=T)
chp_orals4 <- merge(chp_orals4, hp_orals4[,c(1,4,7)], by = "g", all.x=T)

#标注p值
chp_orals4$pn <- 
  ifelse(chp_orals4$NDVI500m_padj < 0.05, "*", "")
chp_orals4$pe <- 
  ifelse(chp_orals4$EVI500m_padj < 0.05, "*", "")

#循环输出
df <- list(chp_orals1, chp_orals2, chp_orals3, chp_orals4)

for (i in 1:4) {
  x <- df[[i]]
  filename <- paste("chp_orals",i,".csv",sep = "")
  write.csv(x, filename)
}



####物种敏感性分析####

####排除+额外调整

#anti,ibd,diabetes,pet,bmi,edu

taxo_hands1 <- data.frame()

taxo_hands2 <- data.frame()

taxo_hands3 <- data.frame()

taxo_hands4 <- data.frame()

taxo_hands5 <- data.frame()

taxo_foreheads1 <- data.frame()

m <- filter(foreheadr,  antibiotic_month != "Yes")
br <- foreheadr_gf[,m$X.SampleID]

dds <- DESeqDataSetFromMatrix(countData = br, 
                              colData = m, 
                              design = ~ factor(sex)+factor(age_cat)+factor(Caucasian)+factor(collection_season)
                              +factor(bmi_cat)+NDVI500meg)
dds <- DESeq(dds,sfType = "poscounts")
res <- results(dds, contrast=c("EVI500meg",'high','low'))
summary(res,0.05)
res <- as.data.frame(res)
res <- dplyr::filter(res, !is.na(res$padj))
names(res) <- paste("EVI500m","_",names(res),sep = "")
res$Row.names <- rownames(res)

taxo_hands5 <- rbind(res, taxo_hands5)

taxo_hands5 <- merge(taxo_hands5,res, by="Row.names",
                      all.x=TRUE,all.y=TRUE)

hp_hands5 <- taxo_hands5[,c("Row.names",
                              "NDVI500m_baseMean","NDVI500m_log2FoldChange","NDVI500m_padj",
                              "EVI500m_baseMean","EVI500m_log2FoldChange","EVI500m_padj")] %>%
  subset(.,NDVI500m_padj<0.05 |EVI500m_padj<0.05)

write.csv(hp_hands5, file = "hp_hands5原始.csv")



Deseq_s1N <- function(x1,x2){
  m <- filter(x1, x2 != "Yes")
  br <- gutr_gf[,m$X.SampleID]
  
  dds <- DESeqDataSetFromMatrix(countData = br, 
                                colData = m, 
                                design = ~ factor(sex)+factor(age_cat)+factor(Caucasian)+factor(collection_season)+factor(pop_binary)
                               +NDVI500meg)
  dds <- DESeq(dds,sfType = "poscounts")
  res <- results(dds, contrast=c("NDVI500meg",'high','low'))
  summary(res,0.05)
  res <- as.data.frame(res)
  res <- dplyr::filter(res, !is.na(res$padj))
  names(res) <- paste("NDVI500m","_",names(res),sep = "")
  res$Row.names <- rownames(res)
  return(res)
}

Deseq_s1E <- function(x1,x2){
  m <- filter(x1, x2 != "Yes")
  br <- gutr_gf[,m$X.SampleID]
  
  dds <- DESeqDataSetFromMatrix(countData = br, 
                                colData = m, 
                                design = ~ factor(sex)+factor(age_cat)+factor(Caucasian)+factor(collection_season)+factor(pop_binary)
                                +EVI500meg)
  dds <- DESeq(dds,sfType = "poscounts")
  res <- results(dds, contrast=c("EVI500meg",'high','low'))
  summary(res,0.05)
  res <- as.data.frame(res)
  res <- dplyr::filter(res, !is.na(res$padj))
  names(res) <- paste("EVI500m","_",names(res),sep = "")
  res$Row.names <- rownames(res)
  return(res)
}



#额外调整公式
Deseq_s2N <- function(x1,x2){
  m <- filter(x1, diet_cat != "Not provided")
  br <- x2[,m$X.SampleID]
  
  dds <- DESeqDataSetFromMatrix(countData = br, 
                                colData = m, 
                                design = ~ factor(sex)+factor(age_cat)+factor(Caucasian)+factor(collection_season)+factor(pop_binary)
                                +factor(diet_cat)+NDVI500meg)
  dds <- DESeq(dds,sfType = "poscounts")
  res <- results(dds, contrast=c("NDVI500meg",'high','low'))
  summary(res,0.05)
  res <- as.data.frame(res)
  res <- dplyr::filter(res, !is.na(res$padj))
  names(res) <- paste("NDVI500m","_",names(res),sep = "")
  res$Row.names <- rownames(res)
  return(res)
}


Deseq_s2E <- function(x1,x2){
  m <- filter(x1, diet_cat != "Not provided")
  br <- x2[,m$X.SampleID]
  
  dds <- DESeqDataSetFromMatrix(countData = br, 
                                colData = m, 
                                design = ~ factor(sex)+factor(age_cat)+factor(Caucasian)+factor(collection_season)+factor(pop_binary)
                                +factor(diet_cat)+EVI500meg)
  dds <- DESeq(dds,sfType = "poscounts")
  res <- results(dds, contrast=c("EVI500meg",'high','low'))
  summary(res,0.05)
  res <- as.data.frame(res)
  res <- dplyr::filter(res, !is.na(res$padj))
  names(res) <- paste("EVI500m","_",names(res),sep = "")
  res$Row.names <- rownames(res)
  return(res)
}




#排除
res1 <- Deseq_s1N(gutr,gutr$diabetes_cat)
res2 <- Deseq_s1E(gutr,gutr$diabetes_cat)
tax_guts3 <-data.frame()
tax_guts3 <- rbind(res1, tax_guts3) %>% 
  merge(.,res2, by="Row.names",all.x=TRUE,all.y=TRUE)

hp_guts3 <- tax_guts3[,c("Row.names",
                            "NDVI500m_baseMean","NDVI500m_log2FoldChange","NDVI500m_padj",
                            "EVI500m_baseMean","EVI500m_log2FoldChange","EVI500m_padj")] %>%
  subset(.,NDVI500m_padj<0.05 |EVI500m_padj<0.05)

write.csv(hp_guts3, file = "hp_guts3原始.csv")

#额外调整
res1 <- Deseq_s2N(foreheadr,foreheadr_gf)
res2 <- Deseq_s2E(foreheadr,foreheadr_gf)
taxo_foreheads6 <-data.frame()
taxo_foreheads6 <- rbind(res1, taxo_foreheads6) %>% 
  merge(.,res2, by="Row.names",all.x=TRUE,all.y=TRUE)

hp_foreheads6 <- taxo_foreheads6[,c("Row.names",
                                    "NDVI500m_baseMean","NDVI500m_log2FoldChange","NDVI500m_padj",
                                    "EVI500m_baseMean","EVI500m_log2FoldChange","EVI500m_padj")] %>%
  subset(.,NDVI500m_padj<0.05 |EVI500m_padj<0.05)

write.csv(hp_foreheads6, file = "hp_foreheads6原始.csv")




hp_oral2 <- hp_oral2[,c(6:12)] %>% filter(., g !="g__"& g !="__") 


hp_oral2$g <- substr(hp_oral2$g,4,nchar(hp_oral2$g))
hp_oral2$NDVI500m_log2FoldChange[is.na(hp_oral2$NDVI500m_log2FoldChange)] <- 0
hp_oral2$g <- factor(hp_oral2$g, 
                     levels=hp_oral2$g[order(hp_oral2$NDVI500m_log2FoldChange, decreasing=FALSE)])
