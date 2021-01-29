## -----------------------------------------------------------------------------
library(contingency)
tab <- rprobMat(10, 2, 3)
tab

## -----------------------------------------------------------------------------
tab[c(1,4,5),]

## -----------------------------------------------------------------------------
tab[,1,1,]

## -----------------------------------------------------------------------------
tab[,1,1,,drop=FALSE]

## -----------------------------------------------------------------------------
margin(tab, 2:3)         # margin of second and third dimensions
conditional(tab, 2, 1)  # second dimension conditional on first

## -----------------------------------------------------------------------------
                         # as above but sequence of cells
margin2(tab, 2:3)        # in table is retained
conditional2(tab, 2, 1)  

## -----------------------------------------------------------------------------
tab2 <- rprobMat(10,2,3)
kl(tab, tab2)   # pairwise Kullback-Leibler divergence
                       # mutual information between
mutualInf(tab, 2, 3)   # second and third dimensions
mutualInf(tab, 2, 3, cond=1)   # conditional mutual information

