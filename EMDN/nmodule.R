
# the goal of the alogrithm is the n-modules in all these networkr 
nmodule <- function(networks,seed_v){

print("Yes, we are in the n-module extraction precedure");
adjmatrix = networks;
adjmatrix[is.na(adjmatrix)] = 0;
num_vertex = dim(networks)[1];
num_network = dim(networks)[3];

for(i in 1:num_netw){ diag(networks[,,i])=0; }



ori_module <- vector(mode='list', length=length(seed_v)); # to store the n-modules
print("yes, we finish the module construction");



#**************** initial procedure ***************************
initial_index = 2;  # 1: maximal 2: unoverlapping;
ori_module <- vector(mode='list', length=length(seed_v)); # to store the n-modules

if(initial_index==1)
{
for (i in 1:length(seed_v))
{ 
  ori_module[[i]]<- list(name=paste("module", i, sep=" "), entropy=100,members=seed_v[i]);


  
   tempmatrix = networks[ori_module[[i]]$members,,];

   # print(dim(tempmatrix)); 7737* 2
   can_weight = rowSums(tempmatrix);
   can_gene_list = order(can_weight,decreasing=TRUE);

 
   ori_module[[i]]$members=c(ori_module[[i]]$members, can_gene_list[1]);
     
   
    tempmatrix = networks[,ori_module[[i]]$members,];
    can_weight = rowSums(tempmatrix[,,1]);
    for(tempi in 2:num_netw){ can_weight = can_weight+rowSums(tempmatrix[,,tempi]);  }
     can_gene_list = order(can_weight,decreasing=TRUE);
    for(tempi in 1:length(can_gene_list)){ 
                     if (!(can_gene_list[tempi] %in% ori_module[[i]]$members))
                        {  ori_module[[i]]$members=c(ori_module[[i]]$members, can_gene_list[tempi]);  break;}  
                        } 

  #  print(ori_module[[i]]$members);


   ori_module[[i]]$members = sort( ori_module[[i]]$members);
   m_entropy=array(0,dim=c(length( ori_module[[i]]$members),7)); #c1:genes;c2-c3: in degree, c4-c5:  degree, c6-c7:entropy
   m_entropy[,1]=ori_module[[i]]$members;                        #c1: genes
   m_entropy[,c(4:5)]=degree[ori_module[[i]]$members,];          #c2-c3: in degree
   for(tempj in 1:num_netw)
   {
       m_entropy[,tempj+1]=rowSums(networks[ori_module[[i]]$members,ori_module[[i]]$members,tempj]); #c
       
     # compute the entropy for each gene
     for(tempj1 in 1:dim(m_entropy)[1])
       {
           if (m_entropy[tempj1,tempj+1] == 0)
             { m_entropy[tempj1,tempj+5] = 2;  }else if ( 2*m_entropy[tempj1,tempj+1] < m_entropy[tempj1,tempj+3])
             { tempprob = m_entropy[tempj1,tempj+1]/m_entropy[tempj1,tempj+3]; 
               m_entropy[tempj1,tempj+5] = 2 + tempprob*log2(tempprob)+(1-tempprob)*log2(1-tempprob);  
             }else{  tempprob = m_entropy[tempj1,tempj+1]/m_entropy[tempj1,tempj+3]; 
                   m_entropy[tempj1,tempj+5] =  -tempprob*log2(tempprob)-(1-tempprob)*log2(1-tempprob);   }

          

       }
      
   }
   #print(ori_module[[i]]$members);
  ori_module[[i]]$p_entropy = m_entropy;
 # print(ori_module[[i]]$p_entropy);
   ori_module[[i]]$v_entropy = sum(m_entropy[,c(6:7)])/num_netw/length(ori_module[[i]]$members);
  
}
}else
{
init_index_vect = rep(0,num_vertex);

for (i in 1:length(seed_v))
{ 
  ori_module[[i]]<- list(name=paste("module", i, sep=" "), entropy=100,members=seed_v[i]);


  
   tempmatrix = networks[ori_module[[i]]$members,,];

   # print(dim(tempmatrix)); 7737* 2
   can_weight = rowSums(tempmatrix);
   can_gene_list = order(can_weight,decreasing=TRUE);
    temp.index = 1;
  while(length(ori_module[[i]]$members)<3)
  {  
     if(!init_index_vect[can_gene_list[temp.index]] ){  ori_module[[i]]$members = c(ori_module[[i]]$members,can_gene_list[temp.index]); init_index_vect[can_gene_list[temp.index]]=1; }
     temp.index = temp.index +1;
  }     
 



   ori_module[[i]]$members = sort( ori_module[[i]]$members);
   m_entropy=array(0,dim=c(length( ori_module[[i]]$members),7)); #c1:genes;c2-c3: in degree, c4-c5:  degree, c6-c7:entropy
   m_entropy[,1]=ori_module[[i]]$members;                        #c1: genes
   m_entropy[,c(4:5)]=degree[ori_module[[i]]$members,];          #c2-c3: in degree
   for(tempj in 1:num_netw)
   {
       m_entropy[,tempj+1]=rowSums(networks[ori_module[[i]]$members,ori_module[[i]]$members,tempj]); #c
       
     # compute the entropy for each gene
     for(tempj1 in 1:dim(m_entropy)[1])
       {
           if (m_entropy[tempj1,tempj+1] == 0)
             { m_entropy[tempj1,tempj+5] = 2;  }else if ( 2*m_entropy[tempj1,tempj+1] < m_entropy[tempj1,tempj+3])
             { tempprob = m_entropy[tempj1,tempj+1]/m_entropy[tempj1,tempj+3]; 
               m_entropy[tempj1,tempj+5] = 2 + tempprob*log2(tempprob)+(1-tempprob)*log2(1-tempprob);  
             }else{  tempprob = m_entropy[tempj1,tempj+1]/m_entropy[tempj1,tempj+3]; 
                   m_entropy[tempj1,tempj+5] =  -tempprob*log2(tempprob)-(1-tempprob)*log2(1-tempprob);   }

          

       }
      
   }
   #print(ori_module[[i]]$members);
  ori_module[[i]]$p_entropy = m_entropy;
 # print(ori_module[[i]]$p_entropy);
   ori_module[[i]]$v_entropy = sum(m_entropy[,c(6:7)])/num_netw/length(ori_module[[i]]$members);
  
}






}
#*************************************************************
#***** The candidates should meet two requirements:
#*****   1. should be highly co-expressed 
#*****   2. sensitivity to the drug 

beta = 0.01;
alpha = 0.1;         # the threshold for the 
min_num_netw = 1;
size_cand = 100;
beta2 = 0.001;
m_size = 200; # the size of the 

 
#**************** expand procedure ***************************


for (i in 1:length(seed_v)) #for checking
{    
  # tempindex indicates the maximum size of a module 
   
   print("the modules");
   print(i);
   tempindex=0;


  while((length(ori_module[[i]]$members)<m_size)&(tempindex<1))
   {

   source("./Cand_neibor.R");
   neibor = Cand_neibor(ori_module[[i]]$members,beta,min_num_netw,size_cand); 
  


    if (length(neibor)==0)
      { 
       print("no candidates"); tempindex=101;  
      }else 
      {
     
     c_node = c(ori_module[[i]]$members,neibor);
     c_matrix = networks[c_node,ori_module[[i]]$members,];
     source("./multip_entropy.R");
     candid =multip_entropy(c_matrix,c_node,ori_module[[i]]$p_entropy); #c_node is the union of $members and neibors
     
      #print(candid[[1]]);
        if ((ori_module[[i]]$v_entropy-candid[[1]])<beta2|is.na(ori_module[[i]]$v_entropy-candid[[1]]))
          { 
            #print("the entropy is: ");
            #print(ori_module[[i]]$v_entropy);
           #print(candid[[1]]);
           #print("the new node");
           # print(candid[[3]]); 
          #  print("no further improvement");  
            tempindex=101; }else
         {  
            # print("the entropy is: ");
            #print(ori_module[[i]]$v_entropy);
            #print(candid[[1]]);
           # print("the new node");
           #  print(candid[[3]]);
            ori_module[[i]]$members = c(ori_module[[i]]$members,candid[[3]]);
             ori_module[[i]]$members = sort( ori_module[[i]]$members);
            tempmatrix = c();
            tempmatrix = candid[[2]];
            ori_module[[i]]$members=sort(tempmatrix[,1]);
            tempmatrix= tempmatrix[order(tempmatrix[,1]),];
            ori_module[[i]]$p_entropy = c();
            ori_module[[i]]$p_entropy = tempmatrix;
            ori_module[[i]]$v_entropy=candid[[1]];
          }  
      
    
      

    }

    
    }

  ori_module[[i]]$matrix = networks[sort(ori_module[[i]]$members),sort(ori_module[[i]]$members),];
 
  print(ori_module[[i]]$members);

}
#************************************************************* 

 

ori_module;

}
