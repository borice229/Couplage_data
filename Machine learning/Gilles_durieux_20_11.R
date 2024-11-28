library(tidyverse)
library(doParallel)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
 
seed=123

# Simulation
n<-1000
sigma<-sqrt(0.1)

xi<-sort(runif(n))

si<- rnorm(n,mean=0, sd=sigma)

m_x <- function(x){
  return (sin(2*pi*x^3))^3
}

yi<-m_x(xi) + si


plot(xi,yi)

# Fonction qui determine l'estimateur

# Noyaux utilisés 
rectangle <- function(x) 
{ 
  0.5*(abs(x)<=1) 
} 

triangle <- function(x) 
{ 
  (1-abs(x))*(abs(x)<=1) 
} 

gauss <- function(x) 
{ 
  (1/(2*pi)^0.5)*exp(-x^2/2) 
} 

epanechnikov <- function(x)
{
  ((3/4)*(1-x^2))*(abs(x)<=1)
}

data <- data.frame(cbind(xi,yi,mx=m_x(xi),si)) 

# Estimateur de Nadaraya-Watson 
Nada_wa <- function(x,Xi,Yi,h,noy="G") 
{ 
  if (noy=="G") 
    noyau=gauss 
  if (noy=="R") 
    noyau=rectangle 
  if (noy=="T") 
    noyau=triangle 
  if (noy=="E") 
    noyau=epanechnikov 
  dXh <- (x-Xi)/h 
  num <- sum( noyau(dXh)*Yi ) 
  denom <- sum( noyau(dXh) ) 
  num/denom 
  return(num/denom) 
} 


# Estimateur de Validation

CV <- function(x,Xi,Yi,h,noy="G"){
  Nada_Wat<-data(cbind(h,mapply(Nada_wa, x,Xi,Yi,h,noy="G")))
  N<-length(h)
  cv<-rep(0,N)
  for (i in N ){
    cv <-sum((Yi -Nada_Wat[i,2])^2)
  }
  return (stock<-which.min(cbind(h,cv)))
}

 

# Calcul pour différents h (h=0.5, h=0.1 et h=0.01)
h <- 0.05
resu1 <- NULL 
for(i in 1: n) 
{ 
  resu1[i] = Nada_wa(data$xi[i],data$xi,data$yi,h) 
} 
lines(data$xi,resu1,col="yellow",lty=1) 

h<-rnorm(10,mean=0,sd=0.01)

CV(data$xi[1],data$xi,data$yi,h)

CV <- function(X,Y,h,noy="G") 
{ 
  E <- NULL 
  for (i in 1:n) 
  { 
    E[i] <- (Y[i]- Nada_wa(X[i],X[-i],Y[-i],h,noy))^2 
  } 
  return(mean(E)) 
} 

hCV <- function(X,Y,noy="G") 
{ 
  CV1 <- function(h) 
  {
    CV(X,Y,h,noy="G") 
    #print(CV(X,Y,h,noy="G"))
  } 
  grilleh=1:200/1000 
  E <- sapply(grilleh,CV1) 
   #E <- mclapply(grilleh,CV1,mc.cores=no_cores)
  return(list(h=grilleh[which.min(E)],erreur=E)) 
}

# V1 ce sont les indivus ici les hurtres
# V2 le temps
# V3 l'ecart en mm 
new_data_V6<- read.table("DBH.txt") %>% 
  filter(V1==6) %>% 
  rename(ind=V1,temp=V2,Ecart=V3)



