
#########################################
### Génération des données
#########################################

MarkovChain=function(N,q,mu1,mu2,sigma){
  X=NULL
  X<-matrix(nrow=N,ncol=1)
  X[1,1]=mu2
  V=rnorm(N,0,sigma)

  for(i in 2:N){
    u=runif(1)
    if(X[i-1,1]==mu1){
      if(u<q[1,1]){X[i,1]=mu1}else X[i,1]=mu2}
    if(X[i-1,1]==mu2){
      if(u<q[2,2]){X[i,1]=mu2}else X[i,1]=mu1}
  }
  Y=X+rnorm(N,0,sigma)
  Y
}

## Fonction g et delta

g = function(y,v){
  c(exp(-((y - Mu[1])^2) / (2*v)), exp(-((y - Mu[2])^2) / (2*v)))
}

delta = function(x){
  1*(x==0)
}

#########################################
### Approx. Filter Update
#########################################

actualiserPhi = function(){
  newPhi = c(0,0)
  for(i in 1:2){
    #numerateur = (phi[1]*q[1,i]*g(y,v)[i]+phi[2]*q[2,i]*g(y,v)[i]) 
    #denominateur = (phi[1]*q[1,1]*g(y,v)[1]+phi[2]*q[2,1]*g(y,v)[1]+phi[1]*q[1,2]*g(y,v)[2]+phi[2]*q[2,2]*g(y,v)[2])
    #newPhi[i] = numerateur / denominateur
    newPhi[i] = (phi[1]*q[1,i]*g(y,v)[i]+phi[2]*q[2,i]*g(y,v)[i]) 
  }
  newPhi/(sum(newPhi))
}


#########################################
### E-Step
#########################################

actualiserR = function(){
  newR = array(0, dim = c(2,2))
  for(i in 1:2){
    for(j in 1:2){
      newR[i,j] = phi[i]*q[i,j]/(phi[1]*q[1,j] + phi[2]*q[2,j])
    }
  }
  newR
}

actualiserRho = function(){
  newRho = array(0, dim = c(2,2,2))
  for(i in 1:2){
    for(j in 1:2){
      for(k in 1:2){
        newRho[i,j,k] = gamma*delta(j-k)*newR[i,j] + (1-gamma)*(Rho[i,j,1]*newR[1,k] + Rho[i,j,2]*newR[2,k])
      }
    }
  }
  newRho
}

actualiserRhoD = function(){
  newRhoD = array(0, dim = c(2,2,3))
  for(i in 1:2){
    for(k in 1:2){
      for(d in 1:3){
        newRhoD[i,k,d] = gamma*delta(i-k)*(y^(d-1)) + (1-gamma)*(RhoD[i,1,d]*newR[1,k] + RhoD[i,2,d]*newR[2,k])
      }
    }
  }
  newRhoD
}


#########################################
### M-Step
#########################################

actualiserS = function(){
  newS = array(0, dim = c(2,2))
  for(i in 1:2){
    for(j in 1:2){
      newS[i,j] = newRho[i,j,1]*newPhi[1] + newRho[i,j,2]*newPhi[2]
    }
  }
  newS
}

actualiserQ = function(){
  q = array(0,dim = c(2,2))
  for(i in 1:2){
    for(j in 1:2){
      q[i,j] = newS[i,j] / (newS[i,1] + newS[i,2])
    }
  }
  q
}

actualiserSD = function(){
  newSD = array(0, dim = c(2,3))
  for(i in 1:2){
    for(d in 1:3){
      newSD[i,d] = newRhoD[i,1,d]*newPhi[1] + newRhoD[i,2,d]*newPhi[2]
    }
  }
  newSD
}

actualiserMu = function(){
  c(newSD[1,2]/ newSD[1,1], newSD[2,2] / newSD[2,1])
}

actualiserV = function(){
  numerateur = newSD[1,3] - (Mu[1]^2)*newSD[1,1] + newSD[2,3] - (Mu[2]^2)*newSD[2,1]
  numerateur /(newSD[1,1] + newSD[2,1])
}








