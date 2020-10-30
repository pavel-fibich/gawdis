gawdis
======

**gawdis** is was designed in the manuscript *Towards a more balanced combination of multiple traits when computing functional differences between species* by Francesco de Bello, Zoltan Botta-Dukat, Jan Leps & Pavel Fibich (under review in MEE).

**gawdis** R package provides 'gawdis' function to produce multi-trait dissimilarity with more uniform contributions of different traits. The approach is based on minimizing the differences in the correlation between the dissimilarity of each trait, or groups of traits, and the multi-trait dissimilarity. This is done using either an analytic or a numerical solution, both available in the function.

**gawdis** will be available on CRAN https://cran.r-project.org/web/packages/gawdis you just need

```
install.packages("gawdis")
library("gawdis")
```

and then run it similarly as gowdis from FD package (https://cran.r-project.org/web/packages/FD).

Examples
======
For nice examples see vignettes on CRAN https://cran.r-project.org/web/packages/gawdis or directly in R :
```
  vignette("gawdis")
```
To check the package follow examples
```
  library(gawdis)
  bodysize<-c(10, 20, 30, 40, 50, NA, 70)
  carnivory<-c(1, 1, 0, 1, 0,1, 0)
  red<-c(1, 0, 0.5, 0, 0.2, 0, 1)
  yellow<-c(0, 1, 0, 0, 0.3, 1, 0)
  blue<-c(0, 0, 0.5,1, 0.5, 0, 0)
  colors.fuzzy<-cbind(red, yellow, blue)
  names(bodysize)<-paste("sp", 1:7, sep="")
  names(carnivory)<-paste("sp", 1:7, sep="")
  rownames(colors.fuzzy)<-paste("sp", 1:7, sep="")
  tall<-as.data.frame(cbind(bodysize, carnivory, colors.fuzzy))
  
  dissim.bodysize<-gowdis(tall[, "bodysize", drop=F])
  dissim.carnivory<-gowdis(tall[, "carnivory", drop=F])
  dissim.colour<-gowdis(tall[, 3:5])/max(gowdis(tall[, 3:5]))
  dissim.bodysize
  dissim.carnivory
  dissim.colour
  
  dall<-list(as.matrix(dissim.bodysize), as.matrix(dissim.carnivory), as.matrix(dissim.colour))
  mean.dissim.all<-as.dist(apply(simplify2array(dall), c(1, 2), mean, na.rm=T), 2)
  mean.dissim.all### this is the correct one
  
  x=tall; w.type="equal"; groups =c(1, 2, 3, 3, 3); fuzzy=TRUE
  (gaw.tallo<-gawdis(tall, w.type="optimized", groups =c(1, 2, 3,3,3), fuzzy=TRUE))
  (gaw.talle<-gawdis(tall, w.type="equal", groups =c(1, 2, 3, 3, 3), fuzzy=TRUE))
  attr(gaw.talle,"weights")
  
  #make groups, when there are traits that are either very much correlated
  #or from the same organs#
  #example from the tussock dataset, with many leaf traits#
  library("FD")
  head(tussock$trait)
  head(tussock$trait[, 3:7])
  cor(tussock$trait[, 3:7], use = "complete")
  #select fewer traits and log-transform when needed#
  tussock.trait<-tussock$trait[, c("height", "LDMC", "leafN","leafS",
  "leafP", "SLA", "seedmass", "raunkiaer", "pollination", "clonality",
  "growthform")]
  tussock.trait.log<-tussock.trait
  #some traits needed log-tranformation, just creating a matrix to
  #store the new data
  tussock.trait.log$height<-log(tussock.trait$height)
  tussock.trait.log$seedmass<-log(tussock.trait$seedmass)
  tussock.trait.log$leafS<-log(tussock.trait$leafS)
  colnames(tussock.trait.log)
  #run the function and test trait contributions#
  #there are NAs so the iteration approach is the only possible#
  #only 20 iterations are used, because of tests#
  #use definitely at least opti.maxiter=100
  gaw.groups<-gawdis(tussock.trait.log, w.type = "optimized",
  opti.maxiter=20,groups.weight=TRUE,groups = c(1,2,2,2,2,2,3,4,5,6,7))
  cors.gaw.gr<-attr(gaw.groups,"correls")
  cors.gaw.gr[12]<-attr(gaw.groups,"group.correls")[2]
  names(cors.gaw.gr)[12]<-"leaves"
  cors.gaw.gr
  #correlation of single traits dissimilarity,
  #including all leaf traits together ("leaves") with
  #multi-trait dissimilarity

  
```  

