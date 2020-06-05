#读取Affymatrix芯片
library(affy)
## 把所有的cel文件都放在一个文件夹(比如/home/jmzeng/your/celfiles文件夹)下面，文件夹地址赋给变量 dir_cels
dir_cels='/home/jmzeng/your/celfiles'
#perform mas5 normalization
affy_data = ReadAffy(celfile.path=dir_cels)
eset.mas5 = mas5(affy_data)
exprSet.nologs = exprs(eset.mas5)
exprSet = log(exprSet.nologs, 2)  #transform to Log_2 if needed
### 上面是是mas5)，下面是rma，两种不同的芯片数据处理方法，最后都是要生成表达矩阵。
library(affy)
data <- ReadAffy(celfile.path=dir_cels) 
eset <- rma(data)
write.exprs(eset,file="data.txt")
#用oligo包来读取
library(oligo)
celFiles<-list.celfiles()
affyRaw <- read.celfiles(celFiles)
library(pd.mogene.2.0.st)  
## 根据芯片平台来载入芯片设计包，没办法自动选择
## 不同的芯片探针包不一样：mogene20sttranscriptcluster.db
eset <- rma(affyRaw)
write.exprs(eset,file="data.txt")