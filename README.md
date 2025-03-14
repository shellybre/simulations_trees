code to create simulated trees with random simulated genomes. The trees given fixed branch lengths. 
The code finds using wasserstein_distance distance between two genomes based on their island distribution. 
The user can choose : 
  lengths it wants to compare to. 
  branch length in tree (left and right). 
  relationship it want to check : between parent and child (always left child). 
  between two brothers. 
  between grandparent and grandchild (always takes opposite branchs = left then right or right then left) 
  number of checks it wants.
it creates bar chart showing the percentage it chose each given length  
