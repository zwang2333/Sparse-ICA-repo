##########################################################################################
########################### Functions supporting Sparse ICA ##############################
############################### Last Updates: 1/5/2023 ###################################
######################## Zihang Wang ##### zihang.wang@emory.edu  ########################
##########################################################################################


##########################################################################################
#Match based on L2 distances
matchICA=function(S,template,M=NULL) {
  require(clue)
  n.comp=ncol(S)
  n.obs=nrow(S)
  if(n.comp>n.obs) warning('Input should be n x d')
  if(!all(dim(template)==dim(S))) warning('Template should be n x d')
  S = t(S)
  template = t(template)
  l2.mat1=matrix(NA,nrow=n.comp,ncol=n.comp)
  l2.mat2=l2.mat1
  for (j in 1:n.comp) {
    for (i in 1:n.comp) {
      l2.mat1[i,j]=sum((template[i,]-S[j,])^2)/n.obs
      l2.mat2[i,j]=sum((template[i,]+S[j,])^2)/n.obs
    }
  }
  l2.mat1=sqrt(l2.mat1)
  l2.mat2=sqrt(l2.mat2)
  l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
  map=as.vector(solve_LSAP(l2.mat))
  l2.1=diag(l2.mat1[,map])
  l2.2=diag(l2.mat2[,map])
  sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2)
  perm=diag(n.comp)[,map]%*%diag(sign.change)

  s.perm=t(perm)%*%S
  if(!is.null(M)) {
    M.perm=t(M)%*%perm
    return(list(S=t(s.perm),M=t(M.perm)))
  }  else {
    t(s.perm)
  }
}

##########################################################################################
angleMatchICA=function(Mx,My,Sx=NULL,Sy=NULL) {
  # match the colums of Mx and Mx, using the 
  # n x p parameterization of the JIN decomposition
  require(clue)
  rx=ncol(Mx)
  ry=ncol(My)
  n.comp = max(rx,ry)
  
  n.obs=nrow(Mx)
  if(n.comp>n.obs) warning('Input should be n x r')
  if(n.obs!=nrow(My)) warning('Mx and My need to have the same number of rows')
  Mx = t(Mx)
  Mx = Mx / sqrt(apply(Mx^2,1,sum))
  
  
  My = t(My)
  My = My / sqrt(apply(My^2,1,sum))
  
  if(rx<ry) {
    Mx = rbind(Mx,matrix(0,ry-rx,n.obs))
  }
  if(ry<rx) {
    My = rbind(My,matrix(0,rx-ry,n.obs))
  }
  angle.mat1=acos(Mx%*%t(My))
  angle.mat2=acos(-1*Mx%*%t(My))
  angle.mat=angle.mat1*(angle.mat1<=angle.mat2)+angle.mat2*(angle.mat2<angle.mat1)
  map=as.vector(solve_LSAP(angle.mat))
  angle.mat1.perm = angle.mat1[,map]
  angle.mat2.perm = angle.mat2[,map]
  angle1=diag(angle.mat1.perm)
  angle2=diag(angle.mat2.perm)
  matchedangles = apply(cbind(angle1,angle2),1,min)
  allangles = angle.mat1.perm*(angle.mat1.perm<=angle.mat2.perm)+angle.mat2.perm*(angle.mat2.perm<angle.mat1.perm)
  sign.change=-1*(angle2<angle1)+1*(angle1<=angle2)
  perm=diag(n.comp)[,map]%*%diag(sign.change)
  
  My.perm=t(perm)%*%My

  # reorder components by their matched angles 
  smatchedangles = sort(matchedangles)
  omangles = order(matchedangles)
  sallangles = allangles[omangles[1:rx],omangles[1:ry]]
  sMx = Mx[omangles[1:rx],]
  sMy.perm = My.perm[omangles[1:ry],]

  if(!is.null(Sy)) {
  Sy.perm=t(perm)%*%Sy
  sSy.perm = Sy.perm[omangles[1:ry],]
  sSx = Sx[omangles[1:rx],]
  return(list(Mx=t(sMx),My = t(sMy.perm),matchedangles = smatchedangles,allangles = sallangles,Sx = sSx, Sy = sSy.perm))
  }
  else {
    return(list(Mx=t(sMx),My = t(sMy.perm),matchedangles = smatchedangles,allangles = sallangles))
  }
}


##########################################################################################
# Match mixing matrices:
# This function does not require M to be square:
frobICA<-function(M1=NULL,M2=NULL,S1=NULL,S2=NULL,standardize=FALSE) {
  #MODEL: X = S M + E, so M is d x p
  #standardize: if standardize==TRUE, then standardizes rows of M1 and M2
  #to have unit norm; if using S1 and S2, standardizes columns to have unit variance.
  #standardize=TRUE makes the measure scale invariant.

  require(clue)
  tfun = function(x) all(x==0)
  if(is.null(M1) && is.null(M2) && is.null(S1) && is.null(S2)) stop("need to supply either M1 and M2 or S1 and S2")
  if(!is.null(M1) && !is.null(M2) && !is.null(S1) && !is.null(S2)) {
    stop("provide either (M1 and M2) or (S1 and S2) but not both (M1,M2) and (S1,S2)")
  }
  if(!is.null(M1) && nrow(M1) > ncol(M1)) stop("The input appears to be S1 and S2, but the arguments were not specified; re-run with S1=<object> and S2=<object>")

  if(is.null(M1)) {
    nS = nrow(S1)
    if(nS!=nrow(S2)) stop('S1 and S2 must have the same number of rows')
    if(sum(apply(S1,2,tfun)) + sum(apply(S2,2,tfun))) stop('frobICA not defined when S1 or S2 has a column of all zeros')
    if(standardize) {
      S1 = scale(S1)
      S2 = scale(S2)
    }
    p = ncol(S1)
    q = ncol(S2)
    if(p < q) {
      S1 = cbind(S1,matrix(0,nS,(q-p)))
    }
    if(q < p) {
      S2 = cbind(S2,matrix(0,nS,(p-q)))
    }
    Stemp = matchICA(S=S1,template=S2)
    n.comp = max(q,p)
    indices = c(1:n.comp)[!(apply(Stemp,2,tfun) | apply(S2,2,tfun))]
    return(sqrt(sum((Stemp[,indices] - S2[,indices])^2))/sqrt(nS*min(p,q)))
  }

  else {
    if(sum(apply(M1,1,tfun)) + sum(apply(M2,1,tfun))) stop('frobICA not defined when M1 or M2 has a row of all zeros')
    if(standardize) {
      temp = diag((diag(M1%*%t(M1)))^(-1/2))
      M1 = temp%*%M1
      temp = diag((diag(M2%*%t(M2)))^(-1/2))
      M2 = temp%*%M2
    }
    p = ncol(M1)
    if(p!=ncol(M2)) stop("M1 and M2 must have the same number of columns")
    d = nrow(M1)
    q = nrow(M2)
    n.comp=max(d,q)
    if(n.comp > p) warning("M should be d x p")
    if(d<q) {
      M1 = rbind(M1,matrix(0,(q-d),p))
    }
    if(q<d) {
      M2 = rbind(M2,matrix(0,(d-q),p))
    }
    l2.mat1=l2.mat2=matrix(NA,nrow=n.comp,ncol=n.comp)
    for (j in 1:n.comp) {
      for (i in 1:n.comp) {
        #since signs are arbitrary, take min of plus and minus:
        l2.mat1[i,j]=sum((M2[i,]-M1[j,])^2)
        l2.mat2[i,j]=sum((M2[i,]+M1[j,])^2)
      }
    }
    l2.mat1=sqrt(l2.mat1)
    l2.mat2=sqrt(l2.mat2)
    #take the min of plus/min l2 distances. This is okay because solve_LSAP is one to one
    l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
    map=as.vector(solve_LSAP(l2.mat))
    #retain relevant l2 distances:
    l2.1=diag(l2.mat1[,map])
    l2.2=diag(l2.mat2[,map])
    #sign.change is for re-ordered matrix 2
    sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2)
    perm=diag(n.comp)[,map]%*%diag(sign.change)
    M.perm=t(perm)%*%M1
    indices = c(1:n.comp)[!(apply(M.perm,1,tfun) | apply(M2,1,tfun))]
    return(sqrt(sum((M.perm[indices,]-M2[indices,])^2))/sqrt(p*min(d,q)))
  }
}



##########################################################################################
# Returns square root of the precision matrix:
covwhitener <- function(X,n.comp=ncol(X),center.row=FALSE) {
  require(MASS)
  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Creates model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  covmat = cov(x.center)
  evdcov = eigen(covmat,symmetric = TRUE)
  whitener = evdcov$vectors%*%diag(1/sqrt(evdcov$values))%*%t(evdcov$vectors)
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=whitener,Z=x.center%*%whitener,mean=apply(X,2,mean)))
}


##########################################################################################
SimFMRI123 = function(snr = 1, noisyICA=TRUE, nTR=50, nImages=1, phi=0.5, dim.data=c(33,33), var.inactive=0) {
  require(neuRosim)
  require(steadyICA)
  m = nImages
  #Latent components are fixed for each simulation:
  x1 = rep(3,5)
  y1 = c(3:7)
  s1.coords = cbind(x1,y1)
  s1 = specifyregion(dim = dim.data, coord = s1.coords, form = "manual")
  s1[s1!=0] = seq(0.5,1,length=length(x1))
  x2 = c(8,8,8,9,10,9,10,10,10,9,8)
  y2 = c(15,14,13,13,13,15,15,16,17,17,17)
  
  s2.coords = cbind(c(x2,x2+7),c(y2,y2))
  s2 = specifyregion(dim=dim.data, coord = s2.coords, form = 'manual')
  s2[s2!=0] = seq(0.5,1,length=2*length(x2))
  
  x3 = c(13,14,15,15,15,14,13,15,15,14,13)
  y3 = c(19,19,19,20,21,21,21,22,23,23,23)
  
  s3.coords = cbind(c(x3,x3+7,x3+14),c(y3,y3,y3))
  s3 = specifyregion(dim=dim.data, coord = s3.coords, form = 'manual')
  s3[s3!=0] = seq(0.5,1,length=3*length(x3))
  
  sim.S = cbind(as.vector(s1),as.vector(s2),as.vector(s3))
  
  if(m>1) {
    t.sim.S = sim.S
    for(i in 1:(m-1)) t.sim.S = rbind(t.sim.S,sim.S)
    sim.S = t.sim.S
    rm(t.sim.S)
  }
  
  ## Add small amount of Gaussian noise to inactive voxels
  # nInactive = sum(sim.S == 0)
  # baseline = rnorm(nInactive,mean=0,sd=sqrt(var.inactive))
  # sim.S[sim.S==0] = baseline
  ## Add Gaussian random field noise to inactive voxels
  if (var.inactive==0){
    nInactive = sum(sim.S == 0)
    baseline = rnorm(nInactive,mean=0,sd=sqrt(var.inactive))
    sim.S[sim.S==0] = baseline
  }else{
    require(scales)
    baseline = spatialnoise(dim = dim.data, sigma=sqrt(var.inactive), nscan = 3, method = "gaussRF", FWHM = 6)
    my_baseline = cbind(rescale(as.vector(baseline[,,1]),c(-0.5,0.5)),rescale(as.vector(baseline[,,2]),c(-0.5,0.5)),rescale(as.vector(baseline[,,3]),c(-0.5,0.5)))
    # sim.S = my_baseline
    sim.S[sim.S==0] = my_baseline[sim.S==0]
  }
  
  
  ##For noise, simulate Gaussian random field. Unique for each simulation:
  if(noisyICA)  nscan = nTR else nscan = nTR-3
  sim.GRF = NULL
  for(k in 1:m) {
    t.sim.GRF <- spatialnoise(dim = dim.data, sigma=1, nscan = nscan, method = "gaussRF", FWHM = 6)
    dim(t.sim.GRF) <- c(prod(dim.data),nscan)
    sim.GRF = rbind(sim.GRF,t.sim.GRF)
  }
  
  ##Mixmat:
  #create timecourses for latent components:
  totaltime <- nTR
  nOnsets = 5+1
  onsets <- seq(from=1, to=totaltime, length=nOnsets)
  dur <- totaltime/10
  #s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, accuracy = 1)
  # use double-gamma ?
  row1 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(1,3)]), durations = list(dur), effectsize = 1, TR = 1, conv = "gamma")
  row2 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,5)]), durations = list(dur), effectsize = 1, TR=1, conv='gamma')
  #NOTE: Time courses can not be identical.
  row3 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,4)]), durations=list(dur), effectsize=1, TR=1, conv='gamma')
  
  sim.Ms = matrix(c(row1,row2,row3),nrow=3,byrow=TRUE)
  sim.Xs = sim.S%*%sim.Ms
  
  if(noisyICA)  {
    sim.Mn = NULL
    sim.Xn = sim.GRF
    for(t in 2:nTR) sim.Xn[,t] = phi*sim.Xn[,t-1]+sim.Xn[,t]
  }  else {
    sim.Mn = matrix(rnorm(nscan*nTR,0,1),nrow=nscan,ncol=nTR)
    for(t in 2:nTR) sim.Mn[,t] = phi*sim.Mn[,t-1] + sim.Mn[,t]
    sim.Xn = sim.GRF%*%sim.Mn
  }
  #sim.Xs = sim.Xs/sqrt(mean(sim.Xs^2)) 
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2)) 
  sim.Xs = sim.Xs/sd(as.vector(sim.Xs)) #standardize so we can control SNR
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn)) 
  sim.Xs = sqrt(snr)*sim.Xs
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)
  
  if(noisyICA) { 
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.Xn, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  } else {
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.GRF, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  }
}


##########################################################################################
# realistic simulation
SimFMRIreal = function(snr = 1, noisyICA=TRUE, nTR=100, nImages=1, phi=0.5, dim.data=59412, var.inactive=0, components=c(6,17,20),true_data,template_xifti,template_time) {
  require(neuRosim)
  #require(steadyICA)
  
  # load true components from group ICA of ABIDE data
  s1 = true_data$estS_sign[,components[1]]
  s2 = true_data$estS_sign[,components[2]]
  s3 = true_data$estS_sign[,components[3]]
  #Latent components are fixed for each simulation:
  sim.S = cbind(as.vector(s1),as.vector(s2),as.vector(s3))
  
  if(nImages>1) {
    t.sim.S = sim.S
    for(i in 1:(m-1)) t.sim.S = rbind(t.sim.S,sim.S)
    sim.S = t.sim.S
    rm(t.sim.S)
  }
  
  ## Add small amount of Gaussian noise to inactive voxels
  nInactive = sum(sim.S == 0)
  baseline = rnorm(nInactive,mean=0,sd=sqrt(var.inactive))
  sim.S[sim.S==0] = baseline
  
  ##For noise, simulate Gaussian random field. Unique for each simulation:
  if(noisyICA)  nscan = nTR else nscan = nTR-3
  sim.GRF = NULL
  for(k in 1:nImages) {
    #t.sim.GRF = spatialnoise(dim = dim.data, sigma=1, nscan = nscan, method = "gaussRF", FWHM = 6)
    #dim(t.sim.GRF) = c(prod(dim.data),nscan)
    #sim.GRF = rbind(sim.GRF,t.sim.GRF)
    sim.GRF = matrix(rnorm(dim.data*nTR),nrow = dim.data)
    template_image$data$cortex_left=sim.GRF[1:29696,]
    template_image$data$cortex_right=sim.GRF[29697:59412,]
    template_image$data$subcort=matrix(0,nrow = 31870,ncol = nTR)
    
    template_smoothed = smooth_xifti(template_image,surf_FWHM = 6)
    sim.GRF[1:29696,] = template_smoothed$data$cortex_left
    sim.GRF[29697:59412,] = template_smoothed$data$cortex_right
  }
  
  # #Mixmat:
  # #create timecourses for latent components:
  # totaltime = nTR
  # nOnsets = 12+1
  # onsets = seq(from=1, to=totaltime, length=nOnsets)
  # dur = totaltime/10
  # #s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, accuracy = 1)
  # row1 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(1,4,6,9)]), durations = list(dur), effectsize = 1, TR = 1, conv = "gamma")
  # row2 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,5,7,11)]), durations = list(dur), effectsize = 1, TR=1, conv='gamma')
  # #NOTE: Time courses can not be identical.
  # row3 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(3,8,10,12)]), durations=list(dur), effectsize=1, TR=1, conv='gamma')
  
  # ##Mixmat:
  # #create timecourses for latent components:
  # totaltime <- nTR
  # nOnsets = 5+1
  # onsets <- seq(from=1, to=totaltime, length=nOnsets)
  # dur <- totaltime/10
  # #s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, accuracy = 1)
  # # use double-gamma ?
  # row1 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(1,3)]), durations = list(dur), effectsize = 1, TR = 1, conv = "gamma")
  # row2 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,5)]), durations = list(dur), effectsize = 1, TR=1, conv='gamma')
  # #NOTE: Time courses can not be identical.
  # row3 <- specifydesign(totaltime = totaltime, onsets = list(onsets[c(2,4)]), durations=list(dur), effectsize=1, TR=1, conv='gamma')

  sim.Ms = template_time[,1:nTR]
  #sim.Ms = matrix(c(row1,row2,row3),nrow=3,byrow=TRUE)
  sim.Xs = sim.S%*%sim.Ms
  
  if(noisyICA)  {
    sim.Mn = NULL
    sim.Xn = sim.GRF
    for(t in 2:nTR) sim.Xn[,t] = phi*sim.Xn[,t-1]+sim.Xn[,t]
  }  else {
    sim.Mn = matrix(rnorm(nscan*nTR,0,1),nrow=nscan,ncol=nTR)
    for(t in 2:nTR) sim.Mn[,t] = phi*sim.Mn[,t-1] + sim.Mn[,t]
    sim.Xn = sim.GRF%*%sim.Mn
  }
  #sim.Xs = sim.Xs/sqrt(mean(sim.Xs^2)) 
  #sim.Xn = sim.Xn/sqrt(mean(sim.Xn^2)) 
  sim.Xs = sim.Xs/sd(as.vector(sim.Xs)) #standardize so we can control SNR
  sim.Xn = sim.Xn/sd(as.vector(sim.Xn)) 
  sim.Xs = sqrt(snr)*sim.Xs
  sim.X = sim.Xs + sim.Xn
  sim.X.whitened = whitener(X=sim.X)
  
  if(noisyICA) { 
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.Xn, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  } else {
    return(list(S = sim.S, Ms = sim.Ms, X=sim.X, Mn = sim.Mn, N = sim.GRF, scaled.S = scale(sim.S),scaled.X = scale(sim.X), whitened.X = sim.X.whitened$Z, whitener = sim.X.whitened$whitener))
  }
}

##########################################################################################
# Sign change for S matrix to image
signchange = function(S,M=NULL) {
  signskew = sign(apply(S,1,function(x) mean(x^3)))
  newS = signskew*S        # t(signskew*t(S))
  if(!is.null(M)) {
    #newM = t(signskew*t(M))
    newM = t(signskew*t(M))
  }else {newM=M}
  return(list(S=newS,M=newM))
}

##########################################################################################
mlcaFP <- function(xData, n.comp = ncol(xData), W.list = NULL, whiten = c('eigenvec','sqrtprec','none'), maxit = 1000, eps = 1e-06, verbose = FALSE, restarts.pbyd = 0, restarts.dbyd = 0, distribution=c('tiltedgaussian','logistic','JB','tanh'), density=FALSE, out.all=FALSE, orth.method=c('svd','givens'), max.comp = FALSE, reinit.max.comp=FALSE, df=0, irlba=FALSE,...) {
  
  #former option:
  #alg.typ = c('parallel','deflation'),
  #alg.typ = match.arg(alg.typ)
  alg.typ = 'parallel'
  
  
  distribution = match.arg(distribution)
  whiten=match.arg(whiten)
  
  if(restarts.dbyd>0 && whiten!='eigenvec') stop('Use whiten=eigenvec with restarts.dbyd')
  ## whiten:
  
  if(irlba) require(irlba)
  if(max.comp) { #if statement evaluates to true for all max.comp!=0
    s.comp = n.comp
    n.comp = max.comp
  }
  if(max.comp=='TRUE') stop('max.comp should be an integer or FALSE')
  if(reinit.max.comp && max.comp==FALSE) stop('Can not reinitialize from max.comp solution if max.comp==FALSE')
  if(reinit.max.comp && alg.typ=='deflation') stop('reinit.max.comp not yet written for deflation algorithm')
  require(ProDenICA)
  #require(multidcov)
  if(distribution=='tiltedgaussian') Gfunc = tiltedgaussian
  if(distribution=='tiltedgaussian' && df==0) stop('df must be greater than 0 for tiltedgaussian')
  if(distribution=='logistic'  && df>0) stop('df should be set to zero when using logistic')
  if(distribution=='logistic') Gfunc = logistic
  if(distribution=='JB'  && df>0) stop('df should be set to zero when using JB')
  if(distribution=='tanh'  && df>0) stop('df should be set to zero when using tanh')
  if(distribution=='tanh') Gfunc = gtanh
  
  if(!is.null(W.list) & class(W.list)!='list') stop('W.list must be a list')
  if(length(W.list) && (restarts.pbyd || restarts.dbyd)) stop('restarts.pbyd and restarts.dbyd must be equal to zero when supplying W.list')
  
  orth.method= match.arg(orth.method)
  p = ncol(xData)
  nRow = nrow(xData)
  d = n.comp
  xData <- scale(xData, center=TRUE, scale=FALSE)
  if (d > p) stop('d must be less than or equal to p')
  if (whiten=='eigenvec') {
    # Use whitener=='eigenvec' so that restarts.dbyd initiates from the
    # span of the first d eigenvectors.
    temp = whitener(X = xData,n.comp = p,irlba=irlba)
    xData = temp$Z
    whitener = temp$whitener
    rm(temp)
  }  else if (whiten=='sqrtprec') {
    est.sigma = cov(xData)  ## Use eigenvalue decomposition rather than SVD.
    evd.sigma = svd(est.sigma)
    whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
    xData = xData%*%whitener
  }
  else {
    whitener = diag(p)
  }
  # warning('TO DO: Look at whitening methods and check for inconsistent options')
  if (is.null(W.list)) {
    if(restarts.pbyd) W.list = gen.inits(p=p,d=d,runs=restarts.pbyd,orth.method=orth.method)
    if(restarts.dbyd) {
      W.temp = gen.inits(p=d,d=d,runs=restarts.dbyd,orth.method=orth.method)
      #pad with zeros:
      zeros = matrix(0,p-d,d)
      W.temp = lapply(W.temp,FUN = function(x) rbind(x,zeros))
      W.list = c(W.list,W.temp)
    }
    
  }
  ## If restarts.pbyd and restarts.dbyd both equal zero:
  if (is.null(W.list)) W.list = gen.inits(p=p,d=d,runs=1,orth.method=orth.method)
  runs = length(W.list)
  out.list = NULL
  loglik.v = numeric(runs)
  for(k in 1:runs) {
    W0 = as.matrix(W.list[[k]])
    if(alg.typ == 'parallel') {
      out.list[[k]] = lca.par(xData=xData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df, ...)
      out.list[[k]]$df = df
    }
    if(alg.typ == 'deflation') {
      out.list[[k]] = lca.def(xData=xData,W0=W0,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=n.comp,df=df,...)
      out.list[[k]]$df = df
    }
    if(max.comp) {
      flist0 = list()
      for (j in 1:d) flist0[[j]] <- Gfunc(out.list[[k]]$S[, j], ...)
      loglik.S <- apply(sapply(flist0, "[[", "Gs"),2,sum)
      orderedLL = order(loglik.S,decreasing=TRUE)
      out.list[[k]]$S = out.list[[k]]$S[,orderedLL[1:s.comp]]
      out.list[[k]]$Ws = out.list[[k]]$Ws[,orderedLL[1:s.comp]]
      out.list[[k]]$loglik = sum(loglik.S[orderedLL[1:s.comp]])
      loglik.v[k] = out.list[[k]]$loglik
    } else {
      loglik.v[k] = out.list[[k]]$loglik
    }
  }
  for(i in 1:runs){
    out.list[[i]]$distribution=distribution
    out.list[[i]]$whitener = whitener
  }
  out = out.list[[which.max(loglik.v)]]
  if(reinit.max.comp) {
    out = lca.par(xData=xData,W0=out$Ws,Gfunc=Gfunc,maxit=maxit,verbose=verbose,density=density,eps=eps,n.comp=s.comp,df=df, ...)
    out$df = df
    out.list[[k+1]] = out
  }
  if(out.all==TRUE) out.list else out
}


##########################################################################################
logistic <- function(xData, scale=sqrt(3)/pi, df=0) {
  #maximizes likelihood given s then calculates gradient w.r.t. w.hat
  #df is not used
  xData = as.vector(xData)
  list(Gs = -xData/scale - log(scale) - 2*log(1+exp(-xData/scale)), gs = -1/scale + 2*exp(-xData/scale)/(scale*(1+exp(-xData/scale))), gps = - 2*exp(-xData/scale) / (scale^2*(1+exp(-xData/scale))^2))
}

##########################################################################################
# jb.stat written by Ze Jin:
jb.stat <- function(x, df=0) {
  n <- length(x)
  s <- sum(x^3)
  k <- sum(x^4)
  Gs <- x^3 * s / n^2 + (x^4 * k / n^2 + 9 / n - 6 * x^4 / n) / 4
  gs <- 6 * x^2 * s / n^2 + (8 * x^3 * (k / n - 3) / n) / 4
  gps <- 6 * (3 * x^4 + 2 * x * s) / n^2 + (24 * x^2 * (k / n - 3) / n + 32 * x^6 / n^2) / 4
  list(Gs = Gs, gs = gs, gps = gps)
}

##########################################################################################
gtanh <- function(s, df=0) {
  a=-1
  #warning('issue with signs -- not equivalent to fastICA, uses negative of G1 from ProDenICA')
  list(Gs = logb(cosh(a * s))/a, gs = tanh(a * s), gps = a * (1 - 
                                                                tanh(a * s)^2))
}


##########################################################################################
##################### Function for generating random starting points #####################
##########################################################################################
gen.inits <- function(p,d,runs,orth.method=c('svd','givens')) {
  orth.method=match.arg(orth.method)
  W.list = list()
  for(i in 1:runs) {
    if(orth.method=='givens') {
      W.list[[i]] <- as.matrix(theta2W(runif(n=choose(p,2),min=0,max=2*pi)))[,1:d]
    } else {
      temp = matrix(rnorm(p*d),p,d)
      W.list[[i]] <- svd(temp)$u
    }
  }
  W.list
}

##########################################################################################
########### Modified Functions for Fast ICA method # Multiple Initialization #############
############################### Last Updates: 3/20/2023 ##################################
##########################################################################################

require(fastICA)

# CALCULATE NEGENTROPY: from Hyvarinen and Oja #
##########################################################################################
calc.negent.hyvarinen = function(s,Gfunc="logcosh",ev.gauss=0.375937268794) {
  if(Gfunc!="logcosh") stop("Add this form of Gfunc to the function")
  if(is.null(ev.gauss)) ev.gauss=mean(log(cosh(rnorm(10000))))
  a=function(x) mean(log(cosh(x)))
  if(!is.null(dim(s))) sum((apply(s,2,a)-ev.gauss)^2) else (mean(log(cosh(s))) - ev.gauss)^2
}
#ENTROPY: normal distribution has a closed form entropy: 0.5*log(abs(det(sigma)))+(d/2)*(1+ log(2*pi)) 
#TO DO: entropy a function of dimension in Hyvarinen et al., p.113

# Functions for implementing multi-start Fast ICA
##########################################################################################
fastICA_restarts = function(X, n.comp, alg.typ = c("parallel","deflation"),
                            fun = c("logcosh","exp"), alpha = 1.0, method = c("R","C"),
                            row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE,
                            restarts = 1,W.list=NULL){
  
  if (restarts==1){
    out = fastICA(X, n.comp, alg.typ,fun, alpha, method,row.norm, maxit, tol,verbose,W.list)
  }else{
    # Randomly make an input for U(0), by default here used orth.method = "svd"
    W.list = gen.inits(p=n.comp, d=n.comp, runs = restarts, orth.method="svd")
    runs = length(W.list)
    # Create a NULL list for storing outcomes
    out.list = NULL
    # Store neg-entropy of Fast ICA
    negentropy = c()
    
    for (k in 1:runs) {
      est.fastICA = fastICA(X, n.comp, alg.typ,fun, alpha, method,row.norm, maxit, tol,verbose,w.init=W.list[[k]])
      # Store the results
      #if (k != 1) {out.list[[k]] = out.list[[1]]}
      out.list[[k]] = est.fastICA
      out.list[[k]]$negentropy = calc.negent.hyvarinen(est.fastICA$S)
      negentropy[k] = out.list[[k]]$negentropy
    }
    out = out.list[[which.max(negentropy)]]
    out$negentropy_restarts = negentropy
  }
  return(out)
}

##########################################################################################
#################################### Whitening Function ##################################
##########################################################################################
whitener <- function(X,n.comp=ncol(X),center.row=FALSE,irlba=FALSE) {
  require(MASS)
  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Corresponds to model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  if(irlba==FALSE) svd.x=svd(x.center,nu=n.comp,nv=n.comp)
  if(irlba==TRUE) {
    require(irlba)
    svd.x=irlba(x.center,nu=n.comp,nv=n.comp)
  }
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=t(ginv(svd.x$v%*%diag(svd.x$d[1:n.comp])/sqrt(n.rep-1))),Z=sqrt(n.rep-1)*svd.x$u,mean=apply(X,2,mean)))
}


##########################################################################################
################## estimate mixing matrix from estimates of components: ##################
##########################################################################################
est.M.ols <- function(sData,xData,intercept=TRUE) {
  if(intercept) coef(lm(xData~sData))[-1,] else coef(lm(xData~sData-1))
}

##########################################################################################
################ soft_thresh to update matrix in relax and split method ##################
##########################################################################################
# Laplace distribution is used, the default nu is 1
soft_thresh_R = function(x, nu = 1, lambda = sqrt(2)/2) {
  xmin = pmax(abs(x)-nu/lambda,0)*sign(x)
  return(xmin)
}

