
# this algorithm is the to rank genes for each cancer networks 

rankgene_new <- function(networks,alpha)
{ 
 print(dim(networks));
 print(length(networks[is.na(networks)]));
networks[is.na(networks)] = 0;
connetworks = networks;
 
num_gene = dim(connetworks)[1];

num_netw = dim(connetworks)[3]


Nnetworks = connetworks;
 
for(i in 1:num_netw)
{ 
  
  tempmatrix = c();
  diag(connetworks[,,i])=0;
  tempmatrix = connetworks[,,i];
  tempmatrix[is.na(tempmatrix)] = 0;
  degree = array(0,dim=c(num_gene));
  degree = rowSums(tempmatrix);
  #print(which(degree==0));
 
  degree = 1/sqrt(degree);
   degree[is.infinite(degree)] = 0
  print(max(degree));
  print(min(degree));
  print(length(degree[is.na(degree)])); 
  degree[is.na(degree)] = 0;
   for(j in 1:num_gene){  tempmatrix[j,]=degree[j]*tempmatrix[j,];   }
   print(length(tempmatrix[is.na(tempmatrix)]));
   for(j in 1:num_gene){  tempmatrix[,j]=degree[j]*tempmatrix[,j];   }
   Nnetworks[,,i] = tempmatrix;
    print(length(tempmatrix[is.na(tempmatrix)]));
  
   
}

 

num_gene = dim(Nnetworks)[1];
 
 
#********************************************************************* 

#******************  The network propagation *************************




lone_out_rank = c();

iteration = 10;


for(tempi1 in 1:2)
{
      temp_gene_rank = array(1,dim=c(num_gene,1));
      
        for(j in 1:iteration)
     {
      temp_gene_rank = Nnetworks[,,tempi1]%*%temp_gene_rank;
      # print(j);
     }
  
 
   lone_out_rank = cbind(lone_out_rank,temp_gene_rank);
 
}


 
for(i in 1:2){lone_out_rank[,i] = scale(lone_out_rank[,i],center=TRUE,scale=TRUE);}
 
ggrank = data.matrix(lone_out_rank);
 
#*********************************************************************

 ggrank;
}


