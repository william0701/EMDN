overlap <- function(n_modules)
{

library(Matrix);
#lname=load("combine_nmodule(0.6_0.01).RData");
num_module = 0;
mname = c();

for(i in 1:length(n_modules))
{

temp.matrix = networks[n_modules[[i]]$members,n_modules[[i]]$members,];
temp.unconnect = c();
for(j in 1:num_network)
{
  temp.matrix1 = temp.matrix;
  diag(temp.matrix1[,,j]) = 0;
  temp.unconnect = c(which(rowSums(temp.matrix1)==0));

}
if(length(temp.unconnect)>0)
{
temp.b = n_modules[[i]]$members;

n_modules[[i]]$members = temp.b[-temp.unconnect];
}
}

for(i in 1:length(n_modules))
{
 if (length(n_modules[[i]]$members)>9)
   {
    num_module=num_module+1;
    mname = c(mname,i);
    }


}
#print("stop1 is okay");   # BREAK POINT1
delete = setdiff(c(1:length(n_modules)),mname);
#print(mname);
if (length(delete)>0)
{
for (i in 1:length(delete))
{
   n_modules[[delete[i]+1-i]]<-NULL;
}
}

#print(length(n_modules));
#print("STOP2 is okay");   # BREAK POINT2

#********************* remove the identical modules ************************

print(length(n_modules));

start=1;
index=1;
while(index)
{
  tempcount = c();
  for(i in (start+1):length(n_modules))
   {
    temp1 = intersect(n_modules[[start]]$members,n_modules[[i]]$members);
    temp2 = unique(union(n_modules[[start]]$members,n_modules[[i]]$members));
    if (length(temp1)>0.1*length(temp2))
       {
        tempcount=c(tempcount,i);
       }   
   }

 if(length(tempcount)>0){n_modules[tempcount]<-NULL;}
   

start=start+1;

if(start<length(n_modules)){index=1;}else{index=0;}
} 

print(length(n_modules));

#***************************************************************************


 
n_modules;
}
