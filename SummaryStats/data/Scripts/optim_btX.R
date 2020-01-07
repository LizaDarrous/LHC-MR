run_optim = function(par){
  
  LHC_SEM <- function(theta,bX,bY,nX,nY,M,ld,wd){
    
    #print(theta)
    ld2 = ld^2

    iX = theta[1]
    iY = theta[2]
    pX = theta[3]
    pU = theta[4]
    pY = theta[5]
    h2X = theta[6]
    h2Y = theta[7]
    tX = 0 
    tY = theta[8]
    a = theta[9]
    b = 0 
    rho = theta[10]
    
    L = NA
    
    T0 = matrix(c(iX/nX, rho, iY/nY), ncol = 3, nrow = 1)
    T1 = (h2X / (M*pX)) * matrix(c(1, a, a^2), ncol = 3, nrow = 1)
    T2 = (1 / (M*pU)) * matrix(c((tX + (b*tY))^2, (tX + (b*tY)) * (tY + (a*tX)), (tY + (a*tX))^2), ncol = 3, nrow = 1)
    T3 = (h2Y / (M*pY)) * matrix(c(b^2, b, 1), ncol = 3, nrow = 1)
    
    detS0 = (T0[1]*T0[3])-(T0[2]^2)
    
    #print(T0);print(T1);print(T2);print(T3);print(detS0);
    
    U = T1
    detA = (U[1]*U[3])-(U[2]^2)
    detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
    det0 = (ld2*detA)+(ld*detB)+detS0
    dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
    l01 = pX*(1-pU)*(1-pY)*dmv_sub
    
    U = T2;
    detA = (U[1]*U[3])-(U[2]^2)
    detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
    det0 = (ld2*detA)+(ld*detB)+detS0
    dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
    l02 = pU*(1-pX)*(1-pY)*dmv_sub
    
    U = T3;
    detA = (U[1]*U[3])-(U[2]^2)
    detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
    det0 = (ld2*detA)+(ld*detB)+detS0
    dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
    l03 = pY*(1-pX)*(1-pU)*dmv_sub
    
    U = T1+T2;
    detA = (U[1]*U[3])-(U[2]^2)
    detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
    det0 = (ld2*detA)+(ld*detB)+detS0
    dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
    l012 = pX*pU*(1-pY)*dmv_sub
    
    U = T1+T3;
    detA = (U[1]*U[3])-(U[2]^2)
    detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
    det0 = (ld2*detA)+(ld*detB)+detS0
    dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
    l013 = pX*pY*(1-pU)*dmv_sub
    
    U = T2+T3;
    detA = (U[1]*U[3])-(U[2]^2)
    detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
    det0 = (ld2*detA)+(ld*detB)+detS0
    dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
    l023 = pU*pY*(1-pX)*dmv_sub
    
    U = T1+T2+T3;
    detA = (U[1]*U[3])-(U[2]^2)
    detB = (U[1]*T0[3]) + (U[3]*T0[1]) - (2*U[2]*T0[2])
    det0 = (ld2*detA)+(ld*detB)+detS0
    dmv_sub = (det0^(-1/2))*exp(  (-1/2)*(1/det0)*( ((ld*U[3]) + T0[3])*(bX^2) + ((ld*U[1])+T0[1])*(bY^2) - (2*((ld*U[2])+T0[2])*(bX*bY))))
    l0123 = pX*pU*pY*dmv_sub
    
    l0 = (1-pX)*(1-pU)*(1-pY)*(detS0^(-1/2))*(exp( (-1/2)*(1/detS0) * ((T0[3]*(bX^2)) + (T0[1]*(bY^2)) - (2*T0[2]*(bX*bY)))));
    
    L = log(l01+l02+l03+l012+l013+l023+l0123+l0)-log(2*pi)
    L = -sum(wd*L)
    if (is.finite(sum(L))){
      return(L)
    }else{
      return(1e6)
    }
    
  }
  theta=unlist(par)
  
  test = optim(theta, LHC_SEM, 
               bX=bX,bY=bY,M=M, nX=nX, nY=nY, ld=ld, wd=wd, 
               method = "L-BFGS-B",
               lower = args.df1[1,], 
               upper = args.df1[2,],
               control = list(maxit = 5e3, 
                              parscale=args.df1[3,]))
  
  test.res=c(test$value,test$par)
  names(test.res)=c("mLL","iX", "iY", "pX","pU","pY","h2X","h2Y","tY","a","rho")
  return(test.res)
}