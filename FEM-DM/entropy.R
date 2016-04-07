
# the entropy value of each module

entropy<-function(module)
{

 num_module = length(module);
 num_network = 2;
  
 mentropy = array(0,dim=c(num_module,num_network+1));
wdegree = degree;
for(i in 1:length(module))
{
  #print(i);
  adjmatrix = networks[module[[i]]$members,module[[i]]$members,];
  modulelength = length(module[[i]]$members);
  for(j in 1:num_network)
  {
    tempmatrix = c();
    tempmatrix = adjmatrix[,,j];
    degree1 = wdegree[module[[i]]$members,j];#diag(tempmatrix); 
    tempmatrix[is.na(tempmatrix)] = 0;
    diag(tempmatrix) = 0;
    indegree = rowSums(tempmatrix);
     #print(degree);
     entropyvalue = 0;
     for(k in 1:modulelength)
        {

          #if(degree[k]==0){prob=0; print(tempmatrix);}else{prob = indegree[k]/degree[k];}
         prob = indegree[k]/degree1[k];
          if(prob==0)
            { 
                entropyvalue=entropyvalue+2; 
            }else if(prob<0.5)
                {   entropyvalue= entropyvalue+2+prob*log2(prob)+(1-prob)*log2(1-prob);  
                }else
                    {    entropyvalue = entropyvalue- prob*log2(prob)-(1-prob)*log2(1-prob);    }
       
        }

      mentropy[i,j]=log2(entropyvalue/modulelength);
 

  }

  
mentropy[i,num_network+1]=sum(mentropy[i,1:num_network])/num_network;



}


mentropy;

} 
