#' @title gawdis function
#'
#' @description \code{gawdis()}, is an extension of the function \link[FD]{gowdis}, in the package FD, for Gower distance (Gower 1971) as fully described in de Bello et al. (2021). It provides a solution to the problem of unequal traits contribution when combining different traits in a multi-trait dissimilarity (without care the contribution of some traits can much stronger than others, i.e. the correlation of the dissimilarity of individual trait, with the multi-trait dissimilarity, will be much stronger for some traits, particularly categorical ones). The solution to this problem is based on minimizing the differences in the correlation between the dissimilarity of each individual trait (or type of traits) and the multi-trait one. Such a task can be resolved analytically or using iterative explorations, depending on the type of data (basically is NA is available only the iterative approach is possible). Both approaches assess ways to provide an equal contribution of traits to the combined trait dissimilarity. Iterative exploration borrows an algorithm from genetic analyses (GA), with the package for genetic algorithms GA, Morrall (2003). This approach is used to minimize standard deviation (SD) of Pearson correlations between the Gower dissimilarity based on single traits and Gower distances combining all traits together, with a proper weight on each variable. GA iteratively explores the space of trait weights by trying several sets of weights (population of candidate solutions), and combines them by processes inspired from the biology (e.g. selection, mutation and crossover) to get new sets of weights (new generation) with better fitness than previous one (Morrall 2003). The best fitness in our case are weights with the minimal SD of correlations. GA is thus doing an optimization, meaning that the more interactions is used the better solution should be found (although still there is a random effect applied), but also greater computing time is necessary. When the \code{groups} are given, first a combined traits distance between species is computed for each group separately as a distance for all traits inside the group together. The computation of the distance depends also on if the traits should be weighted inside the groups. If so, the weights are at the first found by \code{gawdis()} applied on the matrix with just traits inside the group (gawdis() founds the best weights for the group). If traits should not be weighted inside the groups, directly just a standard Gower distances is applied for all traits inside the group.
#'
#' @param x Matrix or data frame containing the variables. Variables can be numeric, ordered, or factor. Symmetric or asymmetric binary variables should be numeric and only contain 0 and 1. Character variables will be converted to factor. NAs are tolerated.
#' @param W Vector listing the weights for the variables in x. W is considered only if \code{w.type} is \code{user}, for \code{w.type="equal"} all weights having the same value and for other w.type’s the weights are computed (see \code{w.type}).
#' @param asym.bin Vector listing the asymmetric binary variables in x.
#' @param ord	Character string specifying the method to be used for ordinal variables (i.e. ordered). "podani" refers to Eqs. 2a-b of Podani (1999), while "metric" refers to his Eq. 3 (see ‘Details’); both options convert ordinal variables to ranks. "classic" simply treats ordinal variables as continuous variables.
#' @param w.type Type of used method. \code{w.type = "analytic"} (default option) – weights optimized by a mathematical algorithm (no NAs are allowed in this option); \code{w.type = "optimized"} – weights optimized by genetic/optimization algorithm based on iteractions; \code{w.type = "equal"} – equal weights, \code{w.type = "user"} – user defined weights are used. Note that is \code{w.type = "analytic"} in case of NAs, the function will apply \code{w.type = "equal"}.
#' @param groups Vector for traits grouping, i.e. defining group of traits that are considered to be reflecting similar biological information (e.g. many leaf traits in plants covering similar information). By default each trait is treated separately (\code{groups = NULL}). In order to define groups use the same values, e.g. \code{groups = c(1,2,2,2,3,3)} in case of 6 variables attributed to 3 groups, with the length of vector that should be the same as \code{ncol(x)}.
#' @param groups.weight Option to weight traits inside the groups. By default it is set to FALSE, all traits inside the groups have the same weights, meaning that some traits will have a greater contribution within the group; TRUE means that gawdis will determine different weights of traits inside the groups, before combining this group with other traits outside the group.
#' @param fuzzy Vector including groups which are defining a single variable, like in the case of fuzzy coding and dummy variables. In this case, use the argument \code{groups} to define which columns belong to the groups. If \code{fuzzy} includes group name (from \code{groups} argument), then the function will transform distances between species within specified group to have maximum value set to 1 (e.g. for \code{groups=c(1,1,2,2,2),fuzzy=c(2)} only distances of group 2 will be transformed). Default is NULL, not to transform distances of any group. Having both \code{groups.weight=TRUE, fuzzy=TRUE} is not possible, therefore \code{!is.null(fuzzy)} leads to overwriting \code{groups.weight} to FALSE.
#' @param silent If to print warnings and detailed information during the computation.
#' @param opti.getSpecDists Allows to use own code that defines the function \code{getSpecDists(tr,gr,gr.weight)} for computing distances between species for each trait (traits are passed as tr argument). It can be given, or pre-defined function doing the same things as gowdis is used (it is not necessary to specify it). If groups and groups.weight arguments are given in gawdis, then they are passed to \code{getSpecDists()} as gr and gr.weight arguments.
#' @param opti.f This is the criteria used to equalize the contribution of traits to the multi-trait dissimilarity. It can be specified. Alternative, by default, the approach is minimizing the differences in the correlations between the dissimilarity on individual trait and the multi-trait approach. Specifically the  1/SD of correlations (SD=standard deviation) is used, i.e. all traits will tend to have a similar correlation with the multi-trait dissimilarity. opti.f is fitness function that is maximalized by genetic algorithm.
#' @param opti.min.weight Set minimum value for weights of traits.
#' @param opti.max.weight Set maximum value for weights of traits.
#' @param opti.maxiter Maximum number of iterations to run before the GA search is halted, see ?ga from GA package. The default is 300 which was found to be quite reliable. The greater numbers increase the computation time.
#'
#' @usage
#' gawdis(x,W = NULL, asym.bin = NULL, ord = c("podani", "metric", "classic"),
#' w.type = c("analytic", "optimized", "equal", "user"), groups = NULL,
#' groups.weight = FALSE, fuzzy = NULL, opti.getSpecDists = NULL,
#' opti.f = NULL,opti.min.weight = 0.01, opti.max.weight = 1,
#' opti.maxiter = 300, silent = FALSE)
#'
#' @keywords gawdis gowdis
#' @return  An object of class dist with the following attributes: Labels, Types (the variable types, where 'C' is continuous/numeric, 'O' is ordinal, 'B' is symmetric binary, 'A' is asymmetric binary, and 'N' is nominal), Size, Metric. Including attributes 1) “correls” with the correlations of each trait with the multi-trait dissimilarity, 2) “weights” for the weights of traits, 3) “group.correls” with weights of groups, 4)”components” with between species transformed distances, and 5) “cor.mat” with correlations between traits.
#' @references  de Bello, F. et al. (2021) Towards a more balanced combination of multiple traits when computing functional differences between species. Methods in Ecology and Evolution, doi: https://doi.org/10.1111/2041-210X.13537.
#'
#' Gower, J. C. (1971) A general coefficient of similarity and some of its properties. Biometrics 27: 857-871.
#'
#' Podani, J. (1999) Extending Gower's general coefficient of similarity to ordinal characters. Taxon 48:331-340.
#'
#' Morrall, D. (2003) Ecological Applications of Genetic Algorithms. Springer, Berlin, Heidelberg.
#'
#' Laliberté, E., and Legendre, P. (2010) A distance-based framework for measuring functional diversity from multiple traits. Ecology 91:299-305.
#'
#' Laliberté, E., Legendre, P., and Shipley, B. (2014). FD: measuring functional diversity from multiple traits, and other tools for functional ecology. R package version 1.0-12. https://cran.r-project.org/package=FD.
#'
#' @seealso \link[FD]{gowdis} from FD package.
#' @examples
#' library(FD) # input data
#' #the gowdis and gawdis functions provide the same results#
#' ex1 <- gowdis(dummy$trait)
#' #using gawdis in the same way as gowdis
#' ex1.gaw1 <- gawdis(dummy$trait, w.type ="equal")
#' \donttest{plot(ex1, ex1.gaw1); abline(0, 1)}
#' #but when doing so, some traits have stronger contribution on the
#' #multi-trait dissimilarity particularly factorial and binary traits#
#' attr(ex1.gaw1, "correls")
#' #correlation of single-trait dissimilarity with multi-trait one#
#'
#' #the gawdis function finds the best weights to equalize trait
#' #contributions this can be done in two ways: analytic=using formulas;
#' #optimized=using iterations both approaches give very similar results
#' #but only the latter can work with NAs#
#' #for the sake of comparisons here NAs are removed#
#' analytical<-gawdis(dummy$trait[,c(2,4,6,8)], w.type ="analytic")
#' #it is not needed to add the argument w.type, this is the approach
#' #used by default if not defined#
#' attr(analytical, "correls")
#' attr(analytical, "weights") #weights finally given to traits
#' iters<-gawdis(dummy$trait[,c(2,4,6,8)], w.type ="optimized", opti.maxiter=2)
#' #here we used 'only' 2 iterations, to speed up the process of tests and
#' #because it better to use at least opti.maxiter=100#
#' attr(iters, "correls")
#' #correlations are not equal, but enough close to each other
#' attr(iters, "weights")
#' \donttest{plot(analytical, iters); abline(0, 1)}
#'
#' #the function can be used also for fuzzy coded/dummy variables traits#
#' #let's create some data#
#' bodysize<-c(10, 20, 30, 40, 50, NA, 70)
#' carnivory<-c(1, 1, 0, 1, 0,1, 0)
#' red<-c(1, 0, 0.5, 0, 0.2, 0, 1)
#' yellow<-c(0, 1, 0, 0, 0.3, 1, 0)
#' blue<-c(0, 0, 0.5,1, 0.5, 0, 0)
#' colors.fuzzy<-cbind(red, yellow, blue)
#' names(bodysize)<-paste("sp", 1:7, sep="")
#' names(carnivory)<-paste("sp", 1:7, sep="")
#' rownames(colors.fuzzy)<-paste("sp", 1:7, sep="")
#' tall<-as.data.frame(cbind(bodysize, carnivory, colors.fuzzy))
#' tall
#' #use groups and fuzzy to treat the 3 columns related to traits
#' #as one traits#
#' gaw.tall<-gawdis(tall, w.type="equal", groups =c(1, 2, 3,3,3),fuzzy=c(3))
#' attr(gaw.tall,"weights")
#' #to get optimized results just change w.type="optimized"
#'
gawdis <- function (x, W=NULL, asym.bin = NULL, ord = c("podani", "metric","classic"),
                     w.type = c("analytic","optimized","equal","user"),
                     groups=NULL, groups.weight=FALSE, fuzzy=NULL,
                     opti.getSpecDists=NULL, opti.f=NULL, opti.min.weight=0.01,
                     opti.max.weight=1, opti.maxiter=300, silent=FALSE)

  {
    if (length(dx <- dim(x)) != 2 || !(is.data.frame(x) || is.numeric(x)))
        stop("x is not a dataframe or a numeric matrix\n")
    if( (!is.null(fuzzy)) & (groups.weight) ){
      warning("(fuzzy!=NULL) and groups.weight=TRUE is not possible, groups.weight is set to FALSE")
      groups.weight <- FALSE
    }
    if( ( !is.null(fuzzy)) & ( length(groups)<1 ) ){
      warning("fuzzy!=NULL requires having groups argument set, distances will not be transformed")
    }

    w.type <- match.arg(w.type)

    if( ( length(groups)>=1 ) & (length(unique(groups)) == 1) & (w.type !="equal") ){
      warning("just one group (~one trait) is automatically setting w.type to equal")
      w.type <- "equal"
    }
    if (!silent) {
      if( length(groups)>=1  ){
        print(paste("Running w.type=",w.type," on groups=c(",paste(groups,collapse = ","),")",sep=""))
      } else {
        print(paste("Running w.type=",w.type, sep=""))
      }
    }
    n <- dx[1]
    p <- dx[2]
    ord <- match.arg(ord)

    varnames <- dimnames(x)[[2]]
    xorig = x

    if (is.data.frame(x)) {
        type <- sapply(x, data.class)
    } else {
        type <- rep("numeric", p)
        names(type) <- colnames(x)
    }

    if (any(type == "character"))
        for (i in 1:p) if (type[i] == "character")
            x[, i] <- as.factor(x[, i])
    is.bin <- function(k) all(k[!is.na(k)] %in% c(0, 1))
    bin.var <- rep(NA, p)
    names(bin.var) <- varnames
    for (i in 1:p) bin.var[i] <- is.bin(x[, i])
    if (any(type[bin.var] != "numeric"))
        stop("Binary variables should be of class 'numeric'\n")
    type[type %in% c("numeric", "integer")] <- 1
    type[type == "ordered"] <- 2
    type[type %in% c("factor", "character")] <- 3
    type[bin.var] <- 4
    if (!is.null(asym.bin)) {
        if (!all(bin.var[asym.bin]))
            stop("Asymetric binary variables must only contain 0 or 1\n")
        else type[asym.bin] <- 5
    }
    type <- as.numeric(type)

    x <- data.matrix(x)
    if (any(type == 2)) {
        if (ord != "classic")
            for (i in 1:p) if (type[i] == 2)
                x[, i] <- rank(x[, i], na.last = "keep")
            else for (i in 1:p) if (type[i] == 2)
                x[, i] <- as.numeric(x[, i])
    }
    range.Data <- function(v) {
        r.Data <- range(v, na.rm = T)
        res <- r.Data[2] - r.Data[1]
        return(res)
    }
    range2 <- apply(x, 2, range.Data)
    comp.Timax <- function(v) {
        Ti.max <- max(v, na.rm = T)
        no.na <- v[!is.na(v)]
        res <- length(no.na[no.na == Ti.max])
        return(res)
    }
    Timax <- apply(x, 2, comp.Timax)
    comp.Timin <- function(v) {
        Ti.min <- min(v, na.rm = T)
        no.na <- v[!is.na(v)]
        res <- length(no.na[no.na == Ti.min])
        return(res)
    }
    Timin <- apply(x, 2, comp.Timin)
    if (ord == "podani") {
        pod <- 1
    } else {
      pod <- 2
    }

# Begining of new part
    if (w.type == "analytic") if (any(is.na(x))) {
        warning("Analytic solution cannot be calculated due to missing trait values. Equal weighting will be applied.")
        w.type<-"equal"
    }

    d.raw<-matrix(NA,nrow=n*(n-1)/2,ncol=p)
    # w <- rep(1, p)/sum(rep(1, p))
    # res <- .C("gowdis", as.double(x), as.double(w), as.integer(type),
    #           as.integer(n), as.integer(p), as.double(range2), as.integer(pod),
    #           as.double(Timax), as.double(Timin), res = double(n *
    #                                                              (n - 1)/2), NAOK = T, PACKAGE = "FD")$res
    for (i in 1:p) {
      if (type[i]==1) d.raw[,i]<-as.numeric(dist(x[,i]/range2[i]))
      if (type[i]==2)
          {
          if (ord != 'podani')  d.raw[,i]<-as.numeric(dist(x[,i]/range2[i]))
          else
             {
             tie<-rep(NA,n)
             for (j in 1:n) tie[j]<-sum(x[,i]==x[j,i],na.rm=T)
             w<-matrix(NA,nrow=n,ncol=n)
             for (j in 1:n) w[j,]<-abs(x[,i]-x[j,i])-((tie[j]-1)/2)-((tie-1)/2)
             w<-(w/(range2[i]-((Timax[i]-1)/2)-((Timin[i]-1)/2)))
             for (j in 1:n) for (k in 1:n) if ( is.na(x[j,i]) | is.na(x[k,i]) ) {
               w[j,k] <-NA
               } else { if ((x[j,i] == x[k,i])) w[j,k] <- 0  }
             d.raw[,i]<-1-as.numeric(as.dist(1-w))
             }
          }
      if (type[i]==3)
          {
          w<-matrix(NA,nrow=n,ncol=n)
          for (j in 1:n) w[j,]<-as.numeric(x[,i]==x[j,i])
          d.raw[,i]<-as.numeric(as.dist(1-w))
          }
      if (type[i]==4) d.raw[,i]<-as.numeric(dist(x[,i]))
      if (type[i]==5)
          {
          w1<-x[,i]%*%t(x[,i])
          w1<-1-w1
          w2<-(1-x[,i])%*%t(1-x[,i])
          is.na(w1)<-(w2==1)
          d.raw[,i]<-as.numeric(as.dist(w1))
          }
    }

    tnames = dimnames(x)[[2]]

    cor.mat<-cor(d.raw,use="pairwise.complete.obs")

#### analytic
    if ( w.type == "analytic" ) {
       if ( !is.null(groups) ) { # have groups
         w<-c()
         d.raw2 <-NULL
         k <- 0
         tnames <- paste("gr",unique(groups),sep="")
         for (i in unique(groups) ) {
         k=k+1
           ii = (1:length(groups))[i == groups]
           if ( length(ii) > 1 ){ # groups of traits
             if ( groups.weight ) { # weight traits inside the groups?
               print ("Traits inside the group were weighted - analytic.")
               group.gaw = gawdis(as.data.frame(xorig[, ii]), w.type="analytic", groups=NULL,groups.weight = F, fuzzy=fuzzy, silent = T )
             } else { # no weighting inside the groups
               print ("Traits inside the group were not weighted - analytic.")
               group.gaw = gawdis(as.data.frame(xorig[, ii]), w.type="equal", groups=rep(i,length(ii)),groups.weight = F, fuzzy=fuzzy, silent = T )
             }
             if ( is.null(d.raw2) ) {
               d.raw2<-matrix(NA,nrow=nrow(d.raw),ncol=length(unique(groups)))
               d.raw2[,k] <- c(group.gaw)
             } else {
               d.raw2[,k] <- c(group.gaw)
             }
             w<-c(w, attr(group.gaw,"weights") )
           }
         }
         w<-w/sum(w)

         cor.mat2<-cor(d.raw2,use="pairwise.complete.obs")
         sigma2<-apply(d.raw2,2,sd,na.rm=T)*(nrow(d.raw2)-1)
         A<-matrix(NA,nrow=length(sigma2),ncol=length(sigma2))
         A[1,]<-1
         for (i in 2:length(sigma2))
           for (j in 1:length(sigma2)) A[i,j]<-sigma2[j]*(cor.mat2[i,j]-cor.mat2[1,j])
         w2<-solve(A,c(1,rep(0,length(sigma2)-1)))

         w.mat<-t(matrix(rep(w,nrow(d.raw)),nrow=p,ncol=nrow(d.raw)))
         d.rawna<-p-rowSums(is.na(d.raw))
         for (dn in which(d.rawna!=p)) w.mat[dn, !is.na(d.raw[dn,]) ] <- w[!is.na(d.raw[dn,])] / sum(w[!is.na(d.raw[dn,])])
         res2<-rowSums(d.raw*w.mat,na.rm=T)
         correls2<-cor(cbind(res2,d.raw),use="pairwise.complete.obs")[-1,1]

         d.raw<-d.raw2
         p<-length(sigma2)
         w3<-w
         w<-w2
         cor.mat<-cor.mat2

       } else { # no groups given
        sigma<-apply(d.raw,2,sd,na.rm=T)*(nrow(d.raw)-1)
        A<-matrix(NA,nrow=p,ncol=p)
        A[1,]<-1
        for (i in 2:p)
           for (j in 1:p) A[i,j]<-sigma[j]*(cor.mat[i,j]-cor.mat[1,j])
        w<-solve(A,c(1,rep(0,p-1)))
        }
    } # end of analytic
#### optimized
    if ( w.type=="optimized") {
      # run GAgawdis
      opti.res <- GAgawdis(tr=xorig, asym.bin=asym.bin, ord=ord, gr=groups, gr.weight = groups.weight,
                           fuzzy=fuzzy, getSpecDists = opti.getSpecDists,
                           f = opti.f,min.weight = opti.min.weight, max.weight = opti.max.weight,
                           maxiter = opti.maxiter)
      w <- c(opti.res$weights)/sum(c(opti.res$weights))
    } # end of w.type=optimized

#### equal
    if (w.type=="equal") {
      if ( !is.null(groups) ) { # have groups
        ugr <- unique(groups)
        w <-c()
        for(oneg in groups) w<-c(w, sum(groups==oneg)/ ( (length(ugr))*sum(groups==oneg)*sum(groups==oneg)) )
      } else {
        w <- rep(1, p)/sum(rep(1, p))
      }
    }
#### user
    if (w.type=="user") {
      if (!missing(W)) {
         if (length(W) != p | !is.numeric(W))
             stop("w needs to be a numeric vector of length = number of variables in x\n")
         if (all(W == 0))
             stop("Cannot have only 0's in 'w'\n")
         w <- W/sum(W)
         }
      else w <- rep(1, p)/sum(rep(1, p))
    }
##### end of w.type specialization
    #####
    if ( ( length(groups)<1 ) | (length(unique(groups))==1) | (w.type %in% c("analytic","optimized") ) ) {
      # groups are solved inside the methods
      #warning("not separating groups")
      w.mat<-t(matrix(rep(w,nrow(d.raw)),nrow=p,ncol=nrow(d.raw)))
      d.rawna<-p-rowSums(is.na(d.raw))
      for (dn in which(d.rawna!=p)) w.mat[dn, !is.na(d.raw[dn,]) ] <- w[!is.na(d.raw[dn,])] / sum(w[!is.na(d.raw[dn,])])
      res<-rowSums(d.raw*w.mat,na.rm=T)
      if ( (!is.null(fuzzy)) && (length(unique(groups))==1) && (unique(groups) %in% fuzzy) ) res <- res/max(res,na.rm=T)
    } else {
      # groups for equal and user must be solved separately
      #warning("separating groups")
      ugr <- unique(groups)
      pgaw <- NULL
      pgaw.na <-NULL # divide distance by just non NA distances
      for (i in ugr ) {
        ii = (1:length(groups))[i == groups]
        if (is.null(pgaw)){
          pgaw <- gawdis(as.data.frame(x[,ii]),W=w[ii],w.type=w.type,groups=rep(i,length(ii)), fuzzy=fuzzy)
          pgaw.na <- !is.na(gowdis(as.data.frame(x[,ii])))
          if ( (!is.null(fuzzy)) && ( i %in% fuzzy ) ) pgaw <- pgaw / max(pgaw,na.rm=T)
        } else {
          onegaw <- gawdis(as.data.frame(x[,ii]),W=w[ii],w.type=w.type,groups=rep(i,length(ii)), fuzzy=fuzzy)
          pgaw.na1 <- !is.na(gowdis(as.data.frame(x[,ii])))
          pgaw.na <- pgaw.na + pgaw.na1
          if ( (!is.null(fuzzy)) && ( i %in% fuzzy ) ) {
            pgaw <- pgaw + (onegaw/ max(onegaw,na.rm=T))
          } else {
            pgaw <- pgaw + onegaw
          }
        }
      }
      res <- pgaw/pgaw.na
    }



    correls<-cor(cbind(res,d.raw),use="pairwise.complete.obs")[-1,1]

    if ( any(w<0) & (!silent) ) warning("Some weights are negative! Try to remove the corresponding trait from the analysis or re-run the function with w.type='optimized'")
# End of new part

    type[type == 1] <- "N"
    type[type == 2] <- "O"
    type[type == 3] <- "C"
    type[type == 4] <- "B"
    type[type == 5] <- "A"
    if (any(is.na(res))) attr(res, "NA.message") <- "NA's in the dissimilarity matrix!"
    attr(res, "Labels") <- dimnames(x)[[1]]
    attr(res, "Size") <- n
    attr(res, "Metric") <- "Gower"
    attr(res, "Types") <- type

# Four new attributes
    names(correls)<- tnames
    attr(res, "correls") <- correls
    names(w)<- tnames
    attr(res, "weights") <- w
    if (!is.null(groups) & (w.type == "analytic") ){
      names(w3)<- dimnames(x)[[2]]
      attr(res, "indiv.weights") <- w3
      attr(res, "indiv.correls") <- correls2
    }
    if (!is.null(groups) & (w.type == "optimized") ){
      d.rawg<-data.frame(row.names = 1:length(c(opti.res$spedis[[1]])) )
      for (sdi in 1:length(opti.res$spedis)) d.rawg[[names(opti.res$spedis)[sdi]]] <- c(opti.res$spedis[[names(opti.res$spedis)[sdi]]])
      correlsg<-cor(cbind(res,d.rawg),use="pairwise.complete.obs")[-1,1]
      attr(res, "group.correls") <- correlsg
    }
    colnames(d.raw)<- tnames
    attr(res, "components") <- d.raw
    colnames(cor.mat)<-tnames
    rownames(cor.mat)<-tnames
    attr(res, "cor.mat") <- cor.mat

    d.raw.f <- apply(d.raw,2,function(x){ table(x)/sum(table(x)) })
    d.raw.a <- unlist(lapply(d.raw.f,function(x) {any(x>0.66)}))
    d.raw.n <- names(d.raw.a)[d.raw.a]

    if ( (length(d.raw.n)>1) & (!silent) ) warning(paste("Consider removing traits: ",paste(d.raw.n,collapse = " "),", because of unbalanced distribution.",sep=""))
    if ( (sd(correls) > 0.05) & (is.null(groups)) & (!silent) ) {
      warning("The effects of invididual traits are higly variable, consider excluding some less variable or unbalanced traits.")
      if (w.type =="optimized") warning("Or consider increase of opti.maxiter.")
    }
    class(res) <- "dist"
    return(res)
}
