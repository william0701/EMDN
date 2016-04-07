


#***************************** Rand module version ************************************** 
lname = load("test-data/1.RData");

source("./rankgene_new.r")
networks = network;
ggrank = rankgene_new(networks,0.2);

save(ggrank,file="1_rank_result.RData")


#print(hello)
num_netw = 2;
num_node = dim(networks)[1];
num_seed = 1000; 
  num_network <<- 2;
ggrank[ggrank>2] = 2;
ggrank[ggrank< -2] = -2;   
generank = rowSums(ggrank);
#seedgene = order(generank,decreasing=TRUE);
seedgene = c(3573);
degree=array(0,dim=c(num_node,num_netw));
#par(mfrow=c(2,2))
#**** here we should transform the modules ************

#threshold = c(0.3,0.2);#threshold = c(0.3,0.20,0.15); #modification


  for(tempi in 1:num_netw)
 {

  diag(networks[,,tempi])=0; 
  degree[,tempi]=rowSums(networks[,,tempi]); 

 }


#******************************************************
 
seeds = sort(seedgene[1:num_seed]);
 
 

source("./nmodule.R");
n_modules =nmodule(networks,seeds);

#nmodule merges once the overlap is 
if(length(n_modules)>1)
{
source("./overlap.R")
n_modules = overlap(n_modules);
}
source("./entropy.R")
entropy_value = entropy(n_modules);
save(n_modules,entropy_value,file="FEM_1.RData");#66
 


 
