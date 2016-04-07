
# this is a fast version of N-modules: cmatrix: the adjacent matrix, cnode: the set of nodes; curentropy: the current entropy for members 

multip_entropy <- function(cmatrix, cnode, curentropy){
  
#*************Step 1: compute the individual partial entropy for the candidate genes ********************
#print("hello welcome to the fast one"); 
num_mem = dim(curentropy)[1];                      # number of members
num_can = length(cnode) - num_mem;                 # number of candidates




 
#*******************Step 2: compute the  entropy changes for each input node ***************************************
nsize = num_mem+1;
Min_evalue = 1000000;
Min_ematrix = c();
Add_node = c();
for(tempi in 1:num_can) # for testing
{

 

can_pentropy = array(0,dim=c(nsize,7));

can_pentropy[1:num_mem,] = curentropy;
can_pentropy[nsize,1] = cnode[num_mem+tempi];
can_pentropy[nsize,c(4:5)] = degree[can_pentropy[nsize,1],];
#print(can_pentropy);
rel_matrix = cmatrix[num_mem+tempi,,];
can_pentropy[nsize,c(2:3)] = colSums(rel_matrix);
#print(rel_matrix);
can_pentropy[1:num_mem,c(2:3)] = can_pentropy[1:num_mem,c(2:3)]+rel_matrix;
#print(can_pentropy);

for(tempi1 in 1:nsize)
  {
     for(tempi2 in 1:num_netw)
        {
        if ( can_pentropy[tempi1,tempi2+1]==0)
            {  can_pentropy[tempi1,tempi2+5] = 2;  }else if ( 2* can_pentropy[tempi1,tempi2+1] <  can_pentropy[tempi1,tempi2+3])
           { tempprob =  can_pentropy[tempi1,tempi2+1]/ can_pentropy[tempi1,tempi2+3]; 
                can_pentropy[tempi1,tempi2+5] = 2 + tempprob*log2(tempprob)+(1-tempprob)*log2(1-tempprob);  
            }else{  tempprob =  can_pentropy[tempi1,tempi2+1]/ can_pentropy[tempi1,tempi2+3]; 
                    can_pentropy[tempi1,tempi2+5] =  -tempprob*log2(tempprob)-(1-tempprob)*log2(1-tempprob);  }
         } # tempi2
   }# tempi1 

temp_entropy_value = sum(can_pentropy[,c(6:7)])/num_netw/(num_mem+1);
temp_entropy_value[is.na(temp_entropy_value)] = 100;

ori_entropy_value = colSums(curentropy[,6:7])/num_mem;
 cur_entropy_value = colSums(can_pentropy[,6:7])/(num_mem+1);
threshold = min(ori_entropy_value-cur_entropy_value);
 if(length(cur_entropy_value[is.na(cur_entropy_value)])==0)
{
if ((temp_entropy_value<Min_evalue)&(threshold>0.001)) { Min_evalue = temp_entropy_value; Min_ematrix = can_pentropy;  Add_node = cnode[tempi+num_mem]; }
}
}
 

#print(Add_node); 

list(Min_evalue,Min_ematrix,Add_node);  
   
 

}
