

source("Em_online_HMM_fonction.R")

Nsimulations = 100
NbEMiterations = 1000

### Pour stocker les valeurs

Q11Array = array(0,dim = c(Nsimulations, NbEMiterations))
Q22Array = array(0,dim = c(Nsimulations, NbEMiterations))
Mu1Array = array(0,dim = c(Nsimulations, NbEMiterations))
Mu2Array = array(0,dim = c(Nsimulations, NbEMiterations))
VArray = array(0,dim = c(Nsimulations, NbEMiterations))

#######################################
### Simulations
#######################################

for(simulation in 1:Nsimulations){
  
  #######################################
  ### Génération des observations
  #######################################
  
  NbObs = 1000
  
  nmin = 20 ## Modifier pour voir ce que ça change
  
  Q<-matrix(data=c(0.95,0.3,0.05,0.7),nrow=2,ncol=2)# to define matrix Q
  mu=matrix(data=c(0,1),nrow=1,ncol=2) # to define vector mu 
  
  Y=MarkovChain(N=NbObs,q=Q,mu1=mu[1],mu2=mu[2],sigma=.5)
  
  #######################################
  ### Algorithme
  #######################################
  
  Mu=c(-0.5,0.5)
  q = matrix(c(0.7,0.5,0.3,0.5), nrow = 2)
  v = 0.5
  nu = c(-0.5, 0.5)  ### nu est la densité de X0 donc soit +- 0.5            
  gamma = 1
  
  ### Initialisation des paramètres 
  phi = c( nu[1]*g(y,v)[1] , nu[2]*g(y,v)[2] ) / (g(y,v)[1] + g(y,v)[2]) 
  
  # Donne Rho[i,j,k] :
  Rho = array(0, dim = c(2,2,2)) 
  
  #Donne RhoD[i,k,d+1]
  RhoD = array(0,dim = c(2,2,3))
  for(i in 1:2){
    for(k in 1:2){
      for(d in 1:3){
        RhoD[i,k,d] = delta(i-k)*(y^(d-1))
      }
    }
  }
  
  #simulation = 1
  
  Q11Array[simulation, 1] = q[1,1]
  Q22Array[simulation, 1] = q[2,2]
  Mu1Array[simulation, 1] = Mu[1]
  Mu2Array[simulation, 1] = Mu[2]
  VArray[simulation, 1] = v
  
  ### C'est parti pour la boucle !
  for(n in 2:NbEMiterations){
    y=Y[n]
    gamma = n^(-0.6)
    
    ### Approx. Filter Update
    newPhi = actualiserPhi() 
    
    ### E-Step
    newR = actualiserR() 
    newRho = actualiserRho()
    newRhoD = actualiserRhoD() 
    
    ### M-Step
    if(n>=nmin){
      newS = actualiserS()
      q = actualiserQ() 
      newSD = actualiserSD() 
      Mu = actualiserMu() 
      v = actualiserV()
    }
    
    phi = newPhi
    Rho = newRho
    RhoD = newRhoD
    
    Q11Array[simulation, n] = q[1,1]
    Q22Array[simulation, n] = q[2,2]
    Mu1Array[simulation, n] = Mu[1]
    Mu2Array[simulation, n] = Mu[2]
    VArray[simulation, n] = v
  }
}

Q11Median = rep(0,NbEMiterations)
Q11Sup = rep(0,NbEMiterations)
Q11Inf = rep(0,NbEMiterations)
for(i in 1:NbEMiterations){
  s = summary(Q11Array[,i])
  Q11Inf[i] = s[2]
  Q11Median[i] = s[4]
  Q11Sup[i] = s[5]
}

Q22Median = rep(0,NbEMiterations)
Q22Sup = rep(0,NbEMiterations)
Q22Inf = rep(0,NbEMiterations)
for(i in 1:NbEMiterations){
  s = summary(Q22Array[,i])
  Q22Inf[i] = s[2]
  Q22Median[i] = s[4]
  Q22Sup[i] = s[5]
}

Mu1Median = rep(0,NbEMiterations)
Mu1Sup = rep(0,NbEMiterations)
Mu1Inf = rep(0,NbEMiterations)
for(i in 1:NbEMiterations){
  s = summary(Mu1Array[,i])
  Mu1Inf[i] = s[2]
  Mu1Median[i] = s[4]
  Mu1Sup[i] = s[5]
}

Mu2Median = rep(0,NbEMiterations)
Mu2Sup = rep(0,NbEMiterations)
Mu2Inf = rep(0,NbEMiterations)
for(i in 1:NbEMiterations){
  s = summary(Mu2Array[,i])
  Mu2Inf[i] = s[2]
  Mu2Median[i] = s[4]
  Mu2Sup[i] = s[5]
}

par(mfrow=c(2,2))

plot(1:NbEMiterations, Q11Median, main = "Q_11", ylim = c(0,1), pch = 16, type='o', col = "red", cex = 0.6)
grid(col = "grey", lty = "dotted", lwd = 1) 
points(1:NbEMiterations, Q11Inf, pch = 20, col = "green", type='l', lwd = 1, cex = 1) 
points(1:NbEMiterations, Q11Sup, pch = 20, col = "green", type='l', lwd = 1, cex = 1) 

plot(1:NbEMiterations, Q22Median, main = "Q_22", ylim = c(0,1), pch = 16, type='o', col = "red", cex = 0.6)
grid(col = "grey", lty = "dotted", lwd = 1) 
points(1:NbEMiterations, Q22Inf, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 
points(1:NbEMiterations, Q22Sup, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 

plot(1:NbEMiterations, Mu1Median, main = "Mu_1", ylim = c(-0.6,0.2), pch = 16, type='o', col = "red", cex = 0.6)
grid(col = "grey", lty = "dotted", lwd = 1) 
points(1:NbEMiterations, Mu1Inf, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 
points(1:NbEMiterations, Mu1Sup, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 

plot(1:NbEMiterations, Mu2Median, main = "Mu_2", ylim = c(0,1.2), pch = 16, type='o', col = "red", cex = 0.6)
grid(col = "grey", lty = "dotted", lwd = 1) 
points(1:NbEMiterations, Mu2Inf, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 
points(1:NbEMiterations, Mu2Sup, pch = 20, col = "green", type='l', lwd = 1, cex = 0.6) 



