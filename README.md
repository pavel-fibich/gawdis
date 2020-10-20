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
```  

