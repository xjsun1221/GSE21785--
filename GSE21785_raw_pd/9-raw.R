# BiocManager::install(c( 'oligo' ),ask = F,update = F)
library(oligo) 

#BiocManager::install(c( 'pd.hg.u133a' ),ask = F,update = F)
library(pd.hg.u133a)

dir='F:/14个gse的分析结果/GSE21785_RAW/'
od=getwd()
setwd(dir)
celFiles <- list.celfiles(listGzipped = T)
celFiles
affyRaw <- read.celfiles( celFiles )
setwd(od)
eSet <- rma(affyRaw)
eSet
# http://math.usu.edu/jrstevens/stat5570/1.4.Preprocess_4up.pdf

exp = exprs(eSet)
pd = pData(eSet)
#处理表达矩阵列名和pd行名，去掉.CEL.gz
identical(rownames(pd),colnames(exp))
library(stringr)
rownames(pd) = str_remove(rownames(pd),".CEL.gz")
colnames(exp) = rownames(pd)
boxplot(exp)
pb2 = boxplot(exp)
#pd有问题，从geo复制
tmp ="GSM542661	Tubulus_LD5
GSM542662	Tubulus_LD6
GSM542663	Tubulus_LD7
GSM542664	Tubulus_LD8
GSM542665	Tubulus_LD9
GSM542666	Tubulus_LD10
GSM542667	Glomerulus_LD5
GSM542668	Glomerulus_LD6
GSM542669	Glomerulus_LD7
GSM542670	Glomerulus_LD8
GSM542671	Glomerulus_LD9
GSM542672	Glomerulus_LD10
"
p1 = str_replace_all(tmp,"\t","\n")%>%
  str_replace_all("\n",",")%>%
  str_split(",",simplify = T)%>%
  as.character()%>%
  str_subset("GSM")
identical(rownames(pd),p1)
p2 = str_replace_all(tmp,"\t","\n")%>%
  str_replace_all("\n",",")%>%
  str_split(",",simplify = T)%>%
  as.character()%>%
  str_subset("LD")
p2
group_list = ifelse(str_detect(p2,"Tubulus"),"Tubulus","Glomerulus")
group_list=factor(group_list,levels = c("Tubulus","Glomerulus"))
pairinfo = factor(c(5:10,5:10))
#exp = log2(exp+1)
#PCA
{
  dat=as.data.frame(t(exp))
  library(FactoMineR)#画主成分分析图需要加载这两个包
  library(factoextra) 
  # pca的统一操作走起
  dat.pca <- PCA(dat, graph = FALSE)
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               #palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  dir.create("GSE21785_RAW")
  ggsave(paste0("GSE21785_RAW","/PCA.png"))
}


#差异分析，用limma包----
library(limma)
design=model.matrix(~group_list+pairinfo)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
head(deg)

#为deg数据框添加几列----
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
#tibble::rownames_to_column()
head(deg)
#2.加symbol列，火山图要用
#id转换，查找芯片平台对应的包
#eSet[[1]]@annotation
#http://www.bio-info-trainee.com/1399.html
#hgu133a
if(!require(hgu133a.db))BiocManager::install("hgu133a.db")
library(hgu133a.db)
ls("package:hgu133a.db")
ids <- toTable(hgu133aSYMBOL)
head(ids)
#merge
deg <- inner_join(deg,ids,by="probe_id")
deg <- deg[!duplicated(deg$symbol),]
head(deg)
#3.加change列：上调或下调，火山图要用

logFC_t=2 #不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
change=ifelse(deg$P.Value>0.01,'stable', 
              ifelse( deg$logFC >logFC_t,'up', 
                      ifelse( deg$logFC < -logFC_t,'down','stable') )
)
deg <- mutate(deg,change)
head(deg)
table(deg$change)
deg <- mutate(deg,v = -log10(P.Value))

#4.加ENTREZID列，后面富集分析要用
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(unique(deg$symbol), fromType = "SYMBOL",
            toType = c( "ENTREZID"),
            OrgDb = org.Hs.eg.db)
head(s2e)
head(deg)
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))

head(deg)

library(dplyr)
dat <- mutate(deg,v=-log10(P.Value))
head(dat)

#exp2,获得行名为基因名的表达矩阵----
exp2 = exp[match(deg$probe_id,rownames(exp)),]
dim(exp2)
rownames(exp2) = deg$symbol
#火山图----
p <- ggplot(data = deg, 
            aes(x = logFC, 
                y = v)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw()
for_label <- deg %>% 
  filter(abs(logFC) >5& P.Value< 0.01)
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
ggsave(paste0("GSE21785_RAW/","volcano.png"))


#配对热图，用exp2画----
x=deg$logFC 
names(x)=deg$symbol
cg=c(names(head(sort(x),20)),
     names(tail(sort(x),20)))
library(pheatmap)
test = data.frame(gsm = colnames(exp),group_list,pairinfo)
test
col = (arrange(test,pairinfo,group_list))$gsm
od = match(col,colnames(exp))
n=exp2[cg,od]

annotation_col=data.frame(group= as.character(group_list)[od],
                          pair = as.character(pairinfo)[od])
rownames(annotation_col)=colnames(n) 

png(filename = "GSE21785_RAW/heatmapB.png")
pheatmap(n,show_colnames =F,
         #show_rownames = F,
         scale = "row",
         cluster_cols = F, 
         annotation_col=annotation_col,
         gaps_col = c(2,4,6,8,10)
) 
dev.off()

# 配对样本的箱线图----

dat <- data.frame(pairinfo=pairinfo,group=group_list,t(exp2))
#配对样本箱线图批量绘制,画10张玩玩

library(ggplot2)
x = colnames(dat)[3:12]
pl = list()
for(i in 1:length(x)){
  pl[[i]] = ggplot(dat, aes_string("group",colnames(dat)[i+2],fill="group")) +
    geom_boxplot() +
    geom_point(size=2, alpha=0.5) +
    geom_line(aes(group=pairinfo), colour="black", linetype="11") +
    xlab("") +
    ylab(paste("Expression of",colnames(dat)[i+2]))+
    theme_classic()+
    theme(legend.position = "none")
}
#拼图
library(patchwork)
pb_top10= wrap_plots(pl,nrow = 2)+plot_annotation(tag_levels = 'A')
ggsave(plot = pb_top10,filename = paste0("box.png"),width = 15,height = 6)
