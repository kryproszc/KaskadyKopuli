library('funkcje')
#DOpasowanie rozkladow do danych
rozklady <- c('norm', 'lnorm', 'gamma', 'weibull', 'exp','logis')

o<-findDistribution(daneAnaliza,rozklady)
margins<-o[[1]]
paramMargins<-o[[2]]

fit_rzecz.1 <- fitdist(daneAnaliza[,1], "norm")
plot(fit_rzecz.1)
hist(daneAnaliza[,2])
fit_rzecz.2 <- fitdist(daneAnaliza[,2], "logis")
plot(fit_rzecz.2)
fit_rzecz.3 <- fitdist(daneAnaliza[,3], "logis")
plot(fit_rzecz.3)
fit_rzecz.4 <- fitdist(daneAnaliza[,4], "norm")
plot(fit_rzecz.4)
fit_rzecz.5 <- fitdist(daneAnaliza[,5], "logis")
plot(fit_rzecz.5)
fit_rzecz.6 <- fitdist(daneAnaliza[,6], "logis")
plot(fit_rzecz.6)
fit_rzecz.7 <- fitdist(daneAnaliza[,7], "weibull")
plot(fit_rzecz.7)


#Wyliczenie wartosci oczekiwanych dla rozkladow
(E.1<-as.numeric(fit_rzecz.1$estimate[1]))
(E.2<-as.numeric(fit_rzecz.2$estimate[1]))
(E.3<-as.numeric(fit_rzecz.3$estimate[1]))
(E.4<-as.numeric(fit_rzecz.4$estimate[1]))
(E.5<-as.numeric(fit_rzecz.5$estimate[1]))
(E.6<-as.numeric(fit_rzecz.6$estimate[1]))
(E.7<-as.numeric(gamma(1+1/fit_rzecz.7$estimate[1])/fit_rzecz.7$estimate[2]))
(E<-c(E.1,E.2,E.3, E.4, E.5, E.6, E.7))
QModuls<-Qmoduls(margins,paramMargins,c(rep(1/7,7)),E)


#Zakladam, ze w kazda linie biznesowa ubezpieczyciel inwestuje tyle samo kapitalu
wagi <-c(rep(1,7)/7)
QDvine<-QcopulaDvine(daneAnaliza,margins,paramMargins,wagi)
QDvine<-QcopulaCvine(daneAnaliza,margins,paramMargins,wagi)


#Wyznaczam SCR oraz efekty dywersyfikacji 
Q <- data.frame(
  Metoda = c('Cvine','Dvine')
  ,Kwantyl = c(QCvine[[2]],QDvine[[2]])
  ,ED = 1-c(QCvine[[2]],QDvine[[2]])/sum(QModuls))
Q

#Wyznaczam wagi, ktore maksymalizuja efekt dywersyfikacji dla kaskady C
wagi_losowe<-solve_random(QCvine[[1]],margins,paramMargins,QModuls[[2]],10000)
wagi_algorytm<-solve_smart(QCvine[[1]],margins,paramMargins,QModuls[[2]],10000,1,0.0001,0.5)
wagi_kombinacje<-solve_propor(QCvine[[1]],c(0,1,0.1))