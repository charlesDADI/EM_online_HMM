
Nsimulations = 100
K = 30 # on fera K itérations de l’algorithme

### Pour stocker les valeurs

Mu1Array = array(0,dim = c(Nsimulations, K))
Mu2Array = array(0,dim = c(Nsimulations, K))
Sigma1Array = array(0,dim = c(Nsimulations, K))
Sigma2Array = array(0,dim = c(Nsimulations, K))

for(simu in 1:Nsimulations){
  
  ################################################################
  # Données issues d'un mélange gaussien 
  ################################################################
  GaussianSmoothing=function(N,p,mu1,mu2,sigma1,sigma2){
    X=NULL
    X<-matrix(data=0,nrow=N,ncol=1)
    ### Boucle permettant de tirer N valeurs issues
    ### d’un  m’elange gaussien :
    for (i in 2:N) {
      if(runif(1)<p){X[i]=rnorm(1,mu1,sigma1)}
      else{X[i]=rnorm(1,mu2,sigma2)}
    }
    X
  }
  X=GaussianSmoothing(N=100,p=0.4,mu1=-2,mu2=2,sigma1=0.5,sigma2=0.5)
  #hist(X,col="lightblue",breaks=100)
  
  #############
  # Algo EM : #
  #############
  #Définition et initialisation des parametres :
  t1=proc.time()# horloge au début de l'algo.
  K = 30 # on fera K itérations de l’algorithme
  
  hatP<-matrix(data=0,nrow=K,ncol=2) #matrice K,2 ..1ere colonne : p, 2nd:1-p
  mu<-matrix(data=0,nrow=K,ncol=2) 
  sigma<-matrix(data=0,nrow=K,ncol=2) 
  
  hatP[1,1]=runif(1) #Initialisation de p uniforme
  hatP[1,2]=runif(1)
  mu[1,1] = runif(1,min=-3,max=0) 
  mu[1,2]= runif(1,min=0,max=2)
  sigma[1,1] = 1
  sigma[1,2] = 1  
  
  for (i in 2:K) {
    
    vrais1 = hatP[i-1,1]*dnorm(X,mean=mu[i-1,1],sd=sigma[i-1,1])
    vrais2 = (1-hatP[i-1,1])*dnorm(X,mean=mu[i-1,2],sd=sigma[i-1,2])
    
    vrais12 = vrais1 / (vrais1 + vrais2) # probas a posteriori p_{i,1}
    vrais22 = vrais2 / (vrais1 + vrais2) # probas a posteriori p_{i,2}
    ## Mise a jour de lambda1 = P(Z=1 | X,Thteta) :
    hatP[i,1]=mean(vrais12)   #p estimé
    hatP[i,2]= 1-hatP[i,1] 
    
    ## Mise  à jour de mu1 et mu2 :
    mu[i,1] = sum(vrais12*X)/sum(vrais12)
    mu[i,2] = sum(vrais22*X)/sum(vrais22)
    ## Mise a jour de sigma1 et sigma2 :
    sigma[i,1] = sqrt(sum(vrais12*(X-mu[i,1])^2)/(sum(vrais12)))
    sigma[i,2] = sqrt(sum(vrais22*(X-mu[i,2])^2)/(sum(vrais22)))
    
    ## Stockage des valeurs
    Mu1Array[simu, i] = mu[i,1]
    Mu2Array[simu, i] = mu[i,2]
    Sigma1Array[simu, i] = sigma[i,1]
    Sigma2Array[simu, i] = sigma[i,2]
  }
}

Mu1Median = rep(0,K)
Mu1Sup = rep(0,K)
Mu1Inf = rep(0,K)
for(i in 1:K){
  s = summary(Mu1Array[,i])
  Mu1Inf[i] = s[2]
  Mu1Median[i] = s[3]
  Mu1Sup[i] = s[5]
}

Mu2Median = rep(0,K)
Mu2Sup = rep(0,K)
Mu2Inf = rep(0,K)
for(i in 1:K){
  s = summary(Mu2Array[,i])
  Mu2Inf[i] = s[2]
  Mu2Median[i] = s[3]
  Mu2Sup[i] = s[5]
}

Sigma1Median = rep(0,K)
Sigma1Sup = rep(0,K)
Sigma1Inf = rep(0,K)
for(i in 1:K){
  s = summary(Sigma1Array[,i])
  Sigma1Inf[i] = s[2]
  Sigma1Median[i] = s[3]
  Sigma1Sup[i] = s[5]
}

Sigma2Median = rep(0,K)
Sigma2Sup = rep(0,K)
Sigma2Inf = rep(0,K)
for(i in 1:K){
  s = summary(Sigma2Array[,i])
  Sigma2Inf[i] = s[2]
  Sigma2Median[i] = s[3]
  Sigma2Sup[i] = s[5]
}

############################################
### Graphiques
############################################

par(mfrow=c(2,2))

plot(1:K, Sigma1Median, main = "Sigma1", ylim = c(0,1), pch = 16, type='o', col = "red", cex = 0.6)
grid(col = "grey", lty = "dotted", lwd = 1) 
points(1:K, Sigma1Inf, pch = 20, col = "green", type='l', lwd = 1, cex = 1) 
points(1:K, Sigma1Sup, pch = 20, col = "green", type='l', lwd = 1, cex = 1) 

plot(1:K, Sigma2Median, main = "Sigma2", ylim = c(0,1), pch = 16, type='o', col = "red", cex = 0.6)
grid(col = "grey", lty = "dotted", lwd = 1) 
points(1:K, Sigma2Inf, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 
points(1:K, Sigma2Sup, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 

plot(1:K, Mu1Median, main = "Mu_1", ylim = c(-2,0), pch = 16, type='o', col = "red", cex = 0.6)
grid(col = "grey", lty = "dotted", lwd = 1) 
points(1:K, Mu1Inf, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 
points(1:K, Mu1Sup, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 

plot(1:K, Mu2Median, main = "Mu_2", ylim = c(0,2), pch = 16, type='o', col = "red", cex = 0.6)
grid(col = "grey", lty = "dotted", lwd = 1) 
points(1:K, Mu2Inf, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 
points(1:K, Mu2Sup, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 

t2=proc.time()-t1 # calcul le temps pour effecuter algo EM 
t2
