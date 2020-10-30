#' @title Internal Genetic Algorithm gawdis function
#'
#' @description Internal part of \code{gawdis()} function for running genetic algorithm
#'
#' @param tr Matrix or data frame containing the variables. Variables can be numeric, ordered, or factor. Symmetric or asymmetric binary variables should be numeric and only contain 0 and 1. Character variables will be converted to factor. NAs are tolerated.
#' @param gr Vector for traits grouping, i.e. defining group of traits that are considered to be reflecting similar biological information (e.g. many leaf traits in plants covering similar information). By default each trait is treated separately (\code{groups=NULL}). In order to define groups use the same values, e.g. \code{groups = c(1,2,2,2,3,3)} in case of 6 variables attributed to 3 groups, with the length of vector that should be the same as \code{ncol(x)}.
#' @param asym.bin Vector listing the asymmetric binary variables in x.
#' @param ord	Character string specifying the method to be used for ordinal variables (i.e. ordered). \code{podani} refers to Eqs. 2a-b of Podani (1999), while "metric" refers to his Eq. 3 (see ‘Details’); both options convert ordinal variables to ranks. "classic" simply treats ordinal variables as continuous variables.
#' @param gr.weight Option to weight traits inside the groups. By default it is set to FALSE, all traits inside the groups have the same weights, meaning that some traits will have a greater contribution within the group; TRUE means that \code{gawdis()} will determine different weights of traits inside the groups, before combining this group with other traits outside the group.
#' @param fuzzy Set to TRUE in case there is a group of columns, in x, which is defining a single variable, like in the case of fuzzy coding and dummy variables. In this case, use the argument \code{groups} to define which columns belong to this group. If set to TRUE the function will make sure distances between species within groups to have maximum value set to 1. Default is FALSE, not to transform between species distances. Having groups.weight and fuzzy set both to TRUE is not possible, therefore fuzzy=TRUE leads to overwriting groups.weight to FALSE.
#' @param getSpecDists Allows to use own code that defines the function \code{getSpecDists(tr,gr,gr.weight)} for computing distances between species for each trait (traits are passed as tr argument). It can be given, or pre-defined function doing the same things as \code{gowdis()} is used (it is not necessary to specify it). If groups and groups.weight arguments are given in gawdis, then they are passed to \code{getSpecDists()} as gr and gr.weight arguments.
#' @param f This is the criteria used to equalize the contribution of traits to the multi-trait dissimilarity. It can be specified. Alternative, by default, the approach is minimizing the differences in the correlations between the dissimilarity on individual trait and the multi-trait approach. Specifically the  1/SD of correlations (SD=standard deviation) is used, i.e. all traits will tend to have a similar correlation with the multi-trait dissimilarity. opti.f is fitness function that is maximalized by genetic algorithm.
#' @param min.weight Set minimum value for weights of traits.
#' @param max.weight Set maximum value for weights of traits.
#' @param maxiter Maximum number of iterations to run before the GA search is halted, see \code{?ga} from GA package. The default is 300 which was found to be quite reliable. The greater numbers increase the computation time.
#' @param monitor If to monit progress of genetic algorithm.
#' @param ... Arguments passed to GA
#'
#' @usage
#' GAgawdis( tr = NULL, asym.bin = NULL, ord = "podani",gr = NULL,
#' gr.weight = FALSE, fuzzy = FALSE, getSpecDists = NULL,
#' f = NULL, min.weight = 0.001, max.weight = 1, maxiter = 300,
#' monitor = FALSE, ... )
#'
#' @keywords gawdis gowdis
#' @return  Returns 'diss' as dissimilarity, weights as solution of GA, ga as GA, spedis as species distance.
#' @examples
#' #GAgawdis() is not exptected to be run directly, but you can try it by
#' \donttest{
#'  library(FD)
#'  GAgawdis(dummy$trait,maxiter=100)
#'  }
#'

GAgawdis <- function(tr=NULL, asym.bin = NULL, ord = "podani",
                  gr=NULL, gr.weight=FALSE, fuzzy=FALSE, getSpecDists=NULL, f=NULL, min.weight=0.001,max.weight=1,
                  maxiter=300, monitor=FALSE, ...) {


  ############### INTERNAL FUNCTION DEFINITIONS if they are not given
  # get list of distances between species for each trait or groups of traits
  if (is.null (getSpecDists) ) {
    getSpecDists <- function (tr=NULL, gr=NULL, gr.weight=FALSE) {

      getIndTraitDist<-function(onetr){
        n<-length(onetr)
        if ( any("ordered" %in% class(onetr)) ) { # ordered factor
          if (ord != "classic"){
            x <- rank(onetr, na.last = "keep")
          } else {
            x <- as.numeric(onetr)
          }
          if (ord != 'podani')  {
            x<- dist(x)/ max(dist(x),na.rm = T)
          } else {
            tie<-rep(NA,n); for (j in 1:n) tie[j]<-sum(x==x[j],na.rm=T)
            w<-matrix(NA,nrow=n,ncol=n)
            for (j in 1:n) w[j,]<-abs(x-x[j])-((tie[j]-1)/2)-((tie-1)/2)
            Timax2 <- sum( x == max(x,na.rm=T),na.rm=T)
            Timin2 <- sum( x == min(x,na.rm=T),na.rm=T)
            w <- w/(max(dist(x),na.rm = T) -((Timax2-1)/2)-((Timin2-1)/2))
            for (j in 1:n) for (k in 1:n) {
              if ( is.na(x[j]) | is.na(x[k]) ) {
                w[j,k] <-NA
              } else { if (x[j] == x[k]) w[j,k] <- 0 }
            }
            x <- 1-as.dist(1-w)
          }
        } else if ( any("factor" %in% class(onetr)) ) { # factor
          x = gowdis(as.data.frame(onetr),asym.bin = asym.bin, ord=ord)
        } else if (all(onetr %in% 0:1)) { # binary
          x = dist(onetr)
        } else {# integer, float
          x= dist(onetr) / max(dist(onetr), na.rm =T)
        }
        return(x)
      }
      speciesdists = list()
      if ( (! is.null(tr))  && ( is.data.frame(tr) ) ){ # tr is data frame
        if (is.null(gr)) { # group not given
          for (i in 1:ncol(tr)) {
            speciesdists[[names(tr)[i]]] = getIndTraitDist(tr[,i])
          }
        } else { # group given
          for (i in unique(gr) ) {
            ii = (1:length(gr))[i == gr]
            if ( length(ii) > 1 ){ # groups of traits
              if ( gr.weight ) { # weight traits inside the groups?
                print ("Traits inside the group were weighted - optimized.")
                group.gaw = GAgawdis( as.data.frame(tr[, ii]) )
                speciesdists[[ paste(names(tr)[ii], collapse = ".gr.") ]] = gowdis( as.data.frame(tr[, ii]),w = group.gaw$weights ,
                                                                                    asym.bin = asym.bin, ord=ord)
              } else { # no weighting inside the groups
                print ("Traits inside the group were not weighted - optimized.")
                ggow = gowdis( as.data.frame(tr[, ii]) , asym.bin = asym.bin, ord=ord)
                if (fuzzy) ggow <- ggow/max(ggow,na.rm=T)
                speciesdists[[ paste(names(tr)[ii], collapse = ".gr.") ]] = ggow
              }
            } else { # individual trait
              speciesdists[[names(tr)[ii]]] = getIndTraitDist(tr[,ii])
            }
          }
        }
      } else {
        print("Argument should exist and be data.frame!")
      }
      return(speciesdists)
    }
  }
  # fitness function 1/std of cors
  if (is.null(f)){
    f <-function (traits, speciesdists, ...){
      # distance between species for all traits together
      dis.gow.w1<-gowdis(traits, w=..., asym.bin = asym.bin, ord=ord)
      # dis.gow.w1 <- mat.multi(as.matrix(spdis),w=...)
      corlist= lapply(speciesdists, FUN=cor, y=dis.gow.w1, use = "pairwise.complete.obs")
      tomin= sd( unlist(corlist) )
      return(1/tomin )
    }
  }
  # groups checking
  if (!is.null(gr) & (length(gr) != ncol(tr)) ) {
    print ("If groups is given, that is should have the same length as ncol(tr)!")
  } else {
    if ( (!is.null(tr))  && is.data.frame(tr) ) {
      speciesdists = getSpecDists(tr,gr,gr.weight)

      GA <- ga(type = "real-valued", fitness =  f, traits=tr, speciesdists=speciesdists,
               lower = c( rep(min.weight, ncol(tr)) ), upper = c( rep(max.weight,ncol(tr)) ),
               maxiter = maxiter, monitor = monitor, ... )
      dis.sol=gowdis(tr, w = GA@solution/sum(GA@solution,na.rm=T),asym.bin = asym.bin, ord=ord )
      #dis.sol <- mat.multi(as.matrix(spdi), w=c(GA@solution))
      return(list (diss=dis.sol, weights=GA@solution, fmax=GA@fitnessValue, ga=GA, spedis=speciesdists))
    } else {
      print("Argument should exist and be data.frame with columns corresponding to traits!")
    }
  }
}
