library(fitdistrplus)
library(numbers)
library(distributions3)
library(VineCopula)
library (TSP) 
library(rsample)
library(Dict)
library(TSdist)
library(ggpubr)
library(tidyverse)  
library(cluster)    
library(factoextra) 
library(dplyr)
library(ggplot2)
library(xlsx)


findDistribution <- function(dataset, distribution){
  chooseParameters<-list()
  chooseDistributions<-list()
  k=0
  for (column in colnames(dataset)){
    k=k+1
    AIC<-list()
    param<-list()
    for (i in 1:length(distribution)){
      p <- as.vector(dataset[column])[,1]
      fit <- fitdist(p, distribution[i])
      AIC[i]<-fit$aic
      param[i]<-list(addtoList(fit$estimate))
    }
    index<-which.min(AIC)
    chooseDistributions[k]<-distribution[index]
    chooseParameters[k]<-param[index]
  }
  return(list(chooseDistributions,chooseParameters))
}


multiplicationVectorMatric<-function(matrix,wagi){
  u<-as.matrix(matrix)
  x<-u
  for(i in 1:ncol(matrix)){
    x[,i]<-wagi[i]*u[,i]
  }
  return(x)
}

getMatixforlines<-function(lines){
  R<-matrix(c(1, 0.5, 0.5 ,0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25, 0.25, 0.25,
              0.5, 1, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25,
              0.5, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.25, 0.5, 0.25,
              0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.5, 0.5, 0.25, 0.5, 0.5,
              0.5, 0.25, 0.25, 0.25, 1, 0.5 ,0.5, 0.25, 0.5 ,0.5, 0.25, 0.25,
              0.25, 0.25, 0.25, 0.25, 0.5,1, 0.5, 0.25, 0.5 ,0.5 ,0.25, 0.25,
              0.5 ,0.5 ,0.25, 0.25, 0.5 ,0.5 ,1, 0.25, 0.5 ,0.5 ,0.25, 0.25,
              0.25, 0.5, 0.5 ,0.5 ,0.25, 0.25, 0.25, 1, 0.5, 0.25, 0.25, 0.5,
              0.5, 0.5 ,0.5 ,0.5, 0.5, 0.5, 0.5 ,0.5, 1, 0.25, 0.5, 0.25,
              0.25, 0.25, 0.25, 0.25, 0.5, 0.5 ,0.5 ,0.25, 0.25, 1, 0.25, 0.25,
              0.25, 0.25, 0.5, 0.5 ,0.25, 0.25, 0.25, 0.25, 0.5 ,0.25, 1, 0.25,
              0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5 ,0.25, 0.25, 0.25, 1)
            ,nrow = 12,ncol = 12)
  return(R[lines,lines])
}


asCall <- function(fun, param){
  cc <-
    if (length(param) == 0)
      quote(FUN(x))
  else if(is.list(param)) {
    as.call(c(quote(FUN), c(quote(x), as.expression(param))))
  } else {
    as.call(c(quote(FUN), c(quote(x), substitute(param))))
  }
  cc[[1]] <- as.name(fun)
  cc
}

quantileDistributions <- function(data,margins,paramMargins,wwagi) {
  dim <- length(margins)
  x<- as.matrix(data)
  u<-x
  #x<-multiplicationVectorMatric(x,wagi)
  for (i in 1:dim) {
    qdf.expr <- asCall(paste0("q", margins[[i]]), paramMargins[[i]])
    x[,i] <- eval(qdf.expr, list(x = u[,i]))
  }
  ww<-multiplicationVectorMatric(x,wwagi)
  Q<-quantile(rowSums(ww), 0.995) - mean(rowSums(ww))
  return(list(x,Q))
}


Qmoduls <- function(margins,paramMargins,wwagi,means_distributions) {
  Qs<-c()
  dim <- length(margins)
  Q = rep(0,dim)
  for (i in 1:dim) {
    qdf.expr <- asCall(paste0("q", margins[[i]]), paramMargins[[i]])
    Qs <- c(Qs,eval(qdf.expr, list(x = 0.995)))
    #Qs <- eval(qdf.expr, list(x = 0.995))
  }
  Qs_wagi<-(Qs-means_distributions)*wwagi
  return(list(Qs_wagi,Qs-means_distributions))
}

QcopulaCvine<-function(data,margins,paramMargins,wagi){
  dane_kopulowe<-pobs (data)
  dane_kop<-data.frame(dane_kopulowe)
  copula<-as.copuladata(dane_kop) 
  familyset=c(1,2,3,4,5,6,7,8,9,10)
  C_vine<- RVineStructureSelect(data =copula ,
                                type = 1,
                                familyset  = familyset, 
                                treecrit ='tau',
                                selectioncrit = 'logLik',
                                rotations = TRUE, method = 'mle' ) 
  familyC <- C_vine$family
  parC<-C_vine$par
  parC2<-C_vine$par2
  MatrixC<-C_vine$Matrix
  set.seed(10)
  
  RVMC = RVineMatrix(Matrix=MatrixC,family=familyC,par=parC,
                     par2=parC2)
  RVMC$type
  set.seed(10)
  Data_C_vine = as.copuladata(RVineSim(10000,RVMC))
  #Data_C_vine<-multiplicationVectorMatric(Data_C_vine,wagi)
  Q<-quantileDistributions(Data_C_vine,margins,paramMargins,wagi)
  return(Q)
}

QcopulaDvine<-function(data,margins,paramMargins,wagi){
  dane_kopulowe<-pobs(data)
  dane_kop<-data.frame(dane_kopulowe)
  copula<-as.copuladata(dane_kop) 
  familyset=c(1,2,3,4,5,6,7,8,9,10)
  d <-dim(copula ) [ 2 ] 
  M<-1-abs(TauMatrix( copula ) )
  hamilton <-TSP::insert_dummy(TSP(M) , label = 'cut' )
  sol <- solve_TSP(hamilton, 
                   method='repetitive_nn' ) 
  order <- cut_tour( sol , 'cut' )
  DVM<- D2RVine(order ,
                family = rep(0 , d*(d-1)/2) ,
                par = rep(0 , d*(d-1)/2)) 
  D_vine<- RVineCopSelect ( data =copula , family =familyset,
                            Matrix = DVM$Matrix , selectioncrit = 'logLik',
                            method = 'mle' , rotations = TRUE)
  
  familyD <- D_vine$family
  parD<-D_vine$par
  parD2<-D_vine$par2
  MatrixD<-D_vine$Matrix
  
  RVMD = RVineMatrix(Matrix=MatrixD,family=familyD,par=parD,par2=parD2)
  RVMD$type
  set.seed(10)
  Data_D_vine = as.copuladata(RVineSim(10000,RVMD))
  #Data_D_vine<-multiplicationVectorMatric(Data_D_vine,wagi)
  
  Q<-quantileDistributions(Data_D_vine,margins,paramMargins,wagi)
  return(Q)
}

#Optymalizuje wagi z jakimi ubezpieczyciel powinien zainwestowac w dana linie biznesowa
#w trojaki sposob

#1 wagi wybieram losowo
solve_random<-function(datainput,Qmod,iter){
  w_best<-runif(7)^4
  best<-w_best/sum(w_best)
  warst<-c(1,1,1,1,1,1,1)
  for(i in 1:iter){
    w_cand<-runif(7)^4
    cand<-w_cand/sum(w_cand)
    x_best<-multiplicationVectorMatric(datainput,best)
    x_warst<-multiplicationVectorMatric(datainput,warst)
    x_cand<-multiplicationVectorMatric(datainput,cand)
    Q_best<-1-(quantile(rowSums(x_best), 0.995) - mean(rowSums(x_best)))/sum(Qmod*best)
    Q_cand<-1-(quantile(rowSums(x_cand), 0.995) - mean(rowSums(x_cand)))/sum(Qmod*cand)
    Q_warst<-1-(quantile(rowSums(x_warst), 0.995) - mean(rowSums(x_warst)))/sum(Qmod*warst)
    if (Q_best>Q_cand){best<-cand}
    else if (Q_warst>Q_cand){warst<-cand}
  }
  return(best)
}

#2 Wybieram wszystkie kombinacje wag

solve_propor<-function(datainput,Qmod,v){
  a<-v[1]
  b<-v[2]
  c<-v[3]
  x <- expand.grid(seq(a,b,c),
                   seq(a,b,c),
                   seq(a,b,c),
                   seq(a,b,c),
                   seq(a,b,c),
                   seq(a,b,c),
                   seq(a,b,c))
  x<-data.frame(x)
  x_combinations <- x[rowSums(x) == 1, ]
  x_combinations<-data.frame(x_combinations)
  tbest<-c(1,0,0,0,0,0,0)
  x_combinations<-x_combinations
  for (i in 1:nrow(x_combinations)){
    candidat<-as.matrix(x_combinations[i,])
    x_best<-multiplicationVectorMatric(datainput,tbest)
    x_cand<-multiplicationVectorMatric(datainput,candidat)
    Q_best<-1-(quantile(rowSums(x_best), 0.995) - mean(rowSums(x_best)))/sum(Qmod*tbest)
    Q_cand<-1-(quantile(rowSums(x_cand), 0.995) - mean(rowSums(x_cand)))/sum(Qmod*candidat)
    if (Q_best<Q_cand){
      tbest<-candidat
    }
  }
  return(tbest)
}

#3 Algorytm optymalizacji napisany przezemnie. Daje najlepsze wyniki spośród 
#stosowanych trzech algorytmow.

solve_smart<-function(datainput,Qmod,smart_inter,smart_start,smart_eps,smart_step){
  n=ncol(datainput)
  current_weight<-c(rep(1/n,n))
  step<-smart_start
  twarst<-c(1,1,1,1,1,1,1)
  while(step>smart_eps){
    for (i in 1:smart_inter){
      ind1 <- sample(1:n, 1)
      ind2 <- sample(1:n, 1)
      if(ind1==ind2 || current_weight[ind2]<step){next}
      proposition<-current_weight
      proposition[ind1] <- proposition[ind1] + step
      proposition[ind2] <- proposition[ind2] - step
      x_current<-multiplicationVectorMatric(datainput,current_weight)
      x_twarst<-multiplicationVectorMatric(datainput,twarst)
      x_proposition<-multiplicationVectorMatric(datainput,proposition)
      Q_current<-1-(quantile(rowSums(x_current), 0.995) - mean(rowSums(x_current)))/sum(Qmod*current_weight)
      Q_twarst<-1-(quantile(rowSums(x_twarst), 0.995) - mean(rowSums(x_twarst)))/sum(Qmod*twarst)
      Q_proposition<-1-(quantile(rowSums(x_proposition), 0.995) - mean(rowSums(x_proposition)))/sum(Qmod*proposition)
      if(Q_proposition<Q_current){current_weight<-proposition}
    }
    step<-step*smart_step
  }
  return(current_weight)
}








