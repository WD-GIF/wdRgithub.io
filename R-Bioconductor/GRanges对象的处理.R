options(stringsAsFactors = F)
suppressPackageStartupMessages({library(BiocStyle)
  library(GenomicRanges)})
df<-data.frame(
  seqnames=rep(c("chr1","chr2","chr1","chr3"),c(1,2,3,4)),
start = c(101, 105, 125, 132, 134, 152, 153, 160, 166, 170),
end = c(104, 120, 133, 132, 155, 154, 159, 166, 171, 190),
strand = rep(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
score = 1:10,
GC = seq(1, 0, length=10),
row.names = head(letters, 10)
)
df
gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
seqinfo(gr) <- Seqinfo(genome="hg38")
seqnames(gr)
ranges(gr)
score(gr)
start(gr)
end(gr)
strand(gr)
seqinfo(gr)
mcols(gr)$score
subset(gr,strand=="+"&score>5,select=GC)#select Onlt metadata
grMod<-gr
grMod[2]<-gr[1]
grMod
gr
rep(gr[2],times=3)
rev(gr)
gr[IRanges(start=c(2,7),end=c(3,9))]#即gr[c(2,3,7,8,9)]
sp<-split(gr,rep(1:2,each=5))#分成两组，每组5条染色体
sp
library(BSgenome.Hsapiens.UCSC.hg19)
getSeq(Hsapiens,sp)
split(gr,~strand)
c(sp[[1]])
stack(sp,index.var="group")#增加一个分组信息
aggregate(gr,score~strand,mean)#按照strand分组，得到一个score的对应值
aggregate(gr,~strand,n_score=lengths(score),mean_score=mean(score))
g<-gr[1:3]
g<-append(g,gr[10])#rbind的替代
flank(g,10)
flank(g,10,start=F)
promoters(g)#一般是上游200bp左右，下游2000bp左右，生成启动子
unstrand(g)
flank(unstrand(g),10)
shift(g,5)
g
resize(g,30)
reduce(gr,ignore.strand=T)#剔除了有重叠的片段和不同类型链的差别，将其合并到一起
gaps(g)
disjoin(g)
cov<-coverage(g)
cov[1:3]
cov_gr<-GRanges(cov)
cov_gr
cov<-coverage(cov_gr,weight="score")
cov
GPos(cov[1:3])
rg <- reduce(gr, with.revmap=TRUE)
rg$score <- mean(extractList(gr$score, rg$revmap))
##################################################
g2<-head(gr,n=2)
g2
union(g,g2)
g2
g
intersect(g,g2)
setdiff(g,g2)
g3<-g[1:2]
ranges(g3[1])<-IRanges(start=105,end=112)
g3
g2
punion(g2,g3)#加p是合并或取序列的差集，并非是直接的看整段序列的始末是否相同
pintersect(g2,g3)
methods(class="GRanges")
set.seed(66+105+111+99+49+56)
pos<-sample(1:200,size=30L)
size<-10L
end<-size+pos-1L
chrom <- sample(paste0("chr", 1:3), size = 30L, replace = TRUE)
query_df <- data.frame(chrom = chrom,
                       start = pos,
                       end = end)
query_dfs <- split(query_df, 1:3)
q1 <- rename(query_dfs[[1L]], start = "pos")
q2 <- rename(query_dfs[[2L]], chrom = "ch", start = "st")
q3 <- rename(query_dfs[[3L]], end = "last")
q1
q1 <- makeGRangesFromDataFrame(q1, start.field = "pos")
q2 <- makeGRangesFromDataFrame(q2, seqnames.field = "ch",
                               start.field = "st")
q3 <- makeGRangesFromDataFrame(q3, end.field = "last")
query <- mstack(q1, q2, q3, .index.var="replicate")
query
subject<-gr
subsetByOverlaps(query,subject,ignore.strand=T)
hits <- findOverlaps(query, subject, ignore.strand=TRUE)
hits
joined <- query[queryHits(hits)]
joined
ranges(joined) <- ranges(pintersect(joined, subject[subjectHits(hits)]))
hitsByQuery <- as(hits, "List")
hitsByQuery
