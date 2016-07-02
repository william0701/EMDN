
# select candidates for seeds with various constrains: 
# cur_node: the current nodes of the n-module
# threshold: the values the candidate must satisfying 
# the strategy: 1. get the submatrix; 2. select only gene only one time based on the submatrix;

Cand_neibor <- function(cur_node,threshold,min_num,size_cand)
{
 
 neibor=c();
 #print("yes, in the candidate");

 
rem_node  = setdiff(c(1:num_node),cur_node);
num_rem_node = length(rem_node);
rem_networks = networks[rem_node,cur_node,]; 
#print(dim(rem_networks));
rem_connect = array(0,dim=c(num_rem_node,num_netw+1));

for(i in 1:num_netw){rem_connect[,i]=rowSums(rem_networks[,,i]);   }
rem_connect[,3]=rowSums(rem_connect[,1:2]);

#******************* First level candidate ***********************************************************
for(i in 1:num_rem_node)
{
     num_count = 0;
     index_count = rep(0,2);
     for(tempi in 1:num_netw){ if((max(rem_connect[i,1:3])>threshold)&(min(rem_connect[i,1:3])>0)){ num_count = num_count+1; }  }
   if(num_count>=min_num) 
   {
    neibor=c(neibor,i); 
  }
}

#print(length(neibor));
#******************* Second level candidate *********************************************************** 
if(length(neibor)>size_cand) 
{
rem_node = rem_node[neibor];

rem_networks = rem_networks[neibor,,];
can_weight = rowSums(rem_networks[,,1]);
for(tempi in 2:num_netw){ can_weight = can_weight+rowSums(rem_networks[,,tempi]);  } 


can_weight_list = order(can_weight,decreasing=TRUE);



neibor = rem_node[can_weight_list[1:size_cand]];
}
else
{
neibor = rem_node[neibor];

}

    

neibor;
}
