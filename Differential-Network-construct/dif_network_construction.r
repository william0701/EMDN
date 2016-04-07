

# load the gene expression data and p-value: 


 
dif_network_construction <- function(geneexpr,pvalue,delta)
{
library(PCIT);
num_gene = dim(geneexpr)[1];
if(num_gene == dim(pvalue)[1])
{



num_gene = dim(geneexpr)[1]
geneexpr[is.na(geneexpr)] = 0;

geneexpr = data.matrix(geneexpr)
#a1 = t(a1)

#network = array(0,dim=c(num_gene,num_gene))
network = cor(geneexpr,method="pearson")
rownames(network) = colnames(geneexpr)
colnames(network) = colnames(geneexpr)
#diag(network) = 0;
#save(network,file="methylation.RData")
print(dim(network))
print(length(which(is.na(network))))
network[network<0]=0;
print("working");

stage1_network = network;
system.time(stage1_result <- pcit(stage1_network));
meaningful.idx.stage1 <- idx(stage1_result);
unmeaningful.idx.stage1 = idxInvert(dim(stage1_network)[1], meaningful.idx.stage1);
stage1_network[unmeaningful.idx.stage1]=0;

network = stage1_network

#************ End of PCIT ******************



#************ Differential network ***************
network[network>delta] = 1;
pvalue = abs(pvalue -1E-50)
t1 = sqrt(pvalue[,1]%*%t(pvalue[,1]))/2

network = network*t1
network = network/max(network)

save(network,file="Dif_network.RData")
 

  }else
  {
  print("genes in expression does not match the gene p-value")

  }



  }


 
