# FEM-DM algorithm User Manual
* __Version__: 1.0 <br>
* __Date__: 09/20/2015 <br>
* __Authors__: Xiaoke Ma, Zaiyi Liu, Lin Gao, Peizhuo Wang, Wanxin Tang<br>
* __Maintainer__: Xiaoke Ma (xkma@xidian.edu.cn)<br>
* __Depends__: R (>2.15.1)<br>
* __License__: GPL (>=2) <br>

##Introduction
The FEM-DM algorithm is designed for identifying methylated gene modules in multiple differential networks. It takes as inputs multiple edge-weighted gene networks and a set of prior probabilities representing the importance of a gene for the conditions under study (no priori information is also practicable). Along with network topological features, the prior probabilities are used to rank and select seeds to initialize module search. 
##Software tutorial
* __step 1: construct the differential co-methylation (co-expression) networks__
```R
  source("dif_network_construction.r")
  dif_network = dif_network_construction(genemeth, pvalue, delta)
  #genemeth: the matrix denoting methyaltion data where row corresponds to gene and 
   column to sample
  #pvalue: the pvalue of gene differential methylation between two cohorts
  # delta : the cutoff parameter for co-methylation (co-expression) network
```


