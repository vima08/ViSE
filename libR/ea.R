f <- function(x) {
  return(dnorm(x))
}

F <- function(x) {
  return(pnorm(x))
}

f_ms <- function(x, mu, sigma) {
  return(f((x-mu)/sigma))
}

F_ms <- function(x, mu, sigma)  {
  return(F((x-mu)/sigma))
}

b_xnp <- function(x, n, p)  {
  return(dbinom(x,n,p))
}

Py_x <-function(x,y, prob = p) {
  result = 0;
  if (x < y) {
    for(i in (x+1):y) {
      result = result +  b_xnp(i,y,prob)
    }
  }
  return(result)
}

F_gamma <-function(gamma, p, q, l) {
  res = trunc(gamma) + 0.5 - p * l;
  res = res / sqrt(p * q * l);
  return(F(-res))
}

muplus <-function(mu, sigma, l, l0) {
  p = F(mu/sigma);
  q = F(-mu/sigma);
  
  l1 = trunc(l0) + 0.5 - p * l;
  l1 = l1 / sqrt(p * q * l);
  
  fm = f(mu/sigma);
  sfs = sigma * fm / sqrt(p*q*l);
  
  result = mu * F(-l1) + sfs * f(l1);
  
  return(result);
}

muplus_old <-function(l, l0) {
  mu = -0.1;
  sigma = 1;
  result = 0;
  p = F(mu/sigma);
  q = F(-mu/sigma);
  for(x in (trunc(l0) + 1):l) {
    result = result + (mu + sigma*f(mu/sigma)*(x/(p*l) - 1)/q)*b_xnp(x,l,p);
  }
  return(result);
}

alpha0 <- function(mu, sigma) {
  p = F(mu/sigma);
  q = F(-mu/sigma);
  return(p*(1 - mu*q/(sigma*f(mu/sigma))))
}

####################################
res = NULL;
for (i in 1:100) {
  res <- c(res, muplus(-0.1, 1, 100, i));
}
plot(res);

write(x, file = "sad", ncolumns = 1)
plot( -30:50, Md_epsilon(-0.1, 1, 100, 50, 0.5, -30:50))

####################################

#-----------------------
x <- seq(-10, 10, length= 30)
y <- x
func1 <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
z <- outer(x, y, func1)
persp(x, y, z, phi = 30, theta = 30, col = "lightblue",
  xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
  main = "Surface elevation data"
)

#-----------------------


#---------------------------------
x = seq(0.01, 1, length= 100)
y = seq(-1, 1, length= 101)
z = outer(x,y, surf_delta_t)
persp(x, y, z, phi = 45, theta = 45, col = "lightblue",
  xlab = "delta 0.01 -> 1", ylab = "t -1 -> 1",
  main = "Md_G"
)


x = seq(0.01, 1, length= 100)
y = seq(0.01, 1, length= 100)
z = outer(x,y, Gsurf_delta_alpha)
persp(x, y, z, phi = 45, theta = 75, col = "lightblue",
      xlab = "delta 0.01 -> 1", ylab = "alpha 0.01 -> 1",
      main = "Md_G"
)

write.csv2(z, file = "C:/Users/Vitaly/Desktop/1task/file1.csv")

#---------------------------------



myOuter <- function(x,y, func) {
  res = matrix(0, nrow = length(x), ncol = length(y))
  X = 1:length(x);
  Y = 1:length(y);
  for(i in X) {
    for(j in Y) {
      res[i,j] = func(x[i],y[j]);
    }
  }
  return(res);
}

myOuter1 <- function(x,y, func, p) {
  res = matrix(0, nrow = length(x), ncol = length(y))
  X = 1:length(x);
  Y = 1:length(y);
  for(i in X) {
    for(j in Y) {
      res[i,j] = func(x[i],y[j],p);
    }
  }
  return(res);
}

Md_epsA <- function(mu, sigma, l, n, alpha, delta) {
  a = n -l;
  p = F(mu/sigma);
  res = mu * t(myOuter1(delta*n, a, Py_x, p)); # Py_x(alpha*n,a,p)
  for(y in 0:(delta*n)) {
    res = res + muplus(mu, sigma, l, pmin(alpha*n, delta*n - y)) * b_xnp(y, a, p);
  }
  return (res)
}

Md_altA <- function(mu, sigma, l, n, alpha, delta) {
  a = n -l;
  gammaN = alpha*n - a;
  p = F(mu/sigma);
  res = mu * t(myOuter1(alpha*n, l, Py_x, p)); # Py_x(alpha*n,a,p)
  for(x in (gammaN+1):(alpha*n)) {
    res = res + muplus(mu, sigma, a, delta*n - x) * b_xnp(x, l, p);
  }
  return (res)
}

Md_altA_cut <- function(l, delta) {
  mu = -0.02;
  sigma = 1;
  alpha = 0.5;
  n = 100;
  #delta = 0.6;
  a = n -l;
  gammaN = alpha*n - a;
  p = F(mu/sigma);
  res = mu * t(myOuter1(alpha*n, l, Py_x, p)); # Py_x(alpha*n,a,p)
  for(x in (gammaN+1):(alpha*n)) {
    res = res + muplus(mu, sigma, a, delta*n - x) * b_xnp(x, l, p);
  }
  return (res)
}


Md_altA_cut1 <- function(l, delta) {
  mu = -0.02;
  sigma = 1;
  alpha = 0.5;
  n = 100;
  #delta = 0.6;
  a = n -l;
  gammaN = alpha*n - a;
  p = F(mu/sigma);
  res = mu * Py_x(alpha*n,a,p)
  for(x in (gammaN+1):(alpha*n)) {
    res = res + muplus(mu, sigma, a, delta*n - x) * b_xnp(x, l, p);
  }
  return (res)
}

Md_altA_cut_Mu_Delta <- function(mu, delta) {
  l = 80;
  sigma = 1;
  alpha = 0.5; # !!! warning
  n = 100;
  a = n -l;
  gammaN = alpha*n - a;
  p = F(mu/sigma);
  res = mu * Py_x(alpha*n,l,p) # a -> l 
  #print(res)
  for(x in (gammaN+1):(alpha*n)) {
    res = res + muplus(mu, sigma, a, delta*n - x) * b_xnp(x, l, p);
  }
  return (res)
}

Md_altA_cut_Mu_DeltaErr <- function(mu, delta) { 
  l = 80;
  sigma = 1;
  alpha = 0.5;
  n = 100;
  a = n -l;
  gammaN = alpha*n - a;
  p = F(mu/sigma);
  res = mu * Py_x(alpha*n,a,p)
  for(x in (gammaN+1):(alpha*n)) {
    res = res + muplus(mu, sigma, a, delta*n - x) * b_xnp(x, l, p);
  }
  return (res)
}



Md_epsA_cut <- function(l, delta) {
  mu = -0.1;
  sigma = 1;
  alpha = 0.5;
  n = 100;
  a = n -l;
  p = F(mu/sigma);
  res = mu * t(myOuter1(trunc(delta*n), a, Py_x, p)); # Py_x(alpha*n,a,p)
  for(y in 0:(delta*n)) {
    res = res + muplus(mu, sigma, l, pmin(alpha*n, delta*n - y)) * b_xnp(y, a, p);
  }
  return (res)
}

Md_epsA_cut_Mu_Delta <- function(mu, delta) {
  l = 80;
  sigma = 1;
  alpha = 0.5; # !!! warning
  n = 100;
  a = n -l;
  p = F(mu/sigma);
  res = mu * t(myOuter1(trunc(delta*n) + 1, a, Py_x, p)); # Py_x(alpha*n,a,p)
  #print(res)
  for(y in 0:(delta*n)) {
    res = res + muplus(mu, sigma, l, pmin(alpha*n, delta*n - y)) * b_xnp(y, a, p);
  }
  return (res)
}

Md_epsA_cut_Mu_OptDelta <- function(mu) {
  res = NULL;
  X = 1:length(mu)
  for(i in X) {
    delta = max(min((alpha0(mu[i], 1) - 0.5) * 10 + 0.5, 1), 0);
    #res <- c(res, delta);
    res <- c(res, Md_epsA_cut_Mu_Delta(mu[i], delta));
  }
  return (res)
}


plot(myOuter(1:100,100, Md_altA_cut))
plot(Md_epsA(-0.05, 1, 1:100, 100, 0.5, 0.6))

plot(myOuter(1:100, 0.6, Md_altA_cut1))

xMu = seq(-1, 1, length= 101)
yDelta = seq(0.4, 1, length= 100)
z = myOuter(xMu,yDelta, Md_altA_cut_Mu_Delta)
persp(xMu, yDelta, z, phi = 45, theta = 75, col = "lightblue",
      xlab = "xMu -1 -> 1", ylab = "yDelta 0.4 -> 1",
      main = "Md_G"
)

write.csv2(z, file = "C:/Users/Vitaly/Desktop/1task/file1.csv")

z1 = myOuter(xMu,yDelta, Md_altA_cut_Mu_Delta)
z = myOuter(xMu,yDelta, Md_epsA_cut_Mu_Delta)
write.csv2(z, file = "C:/Users/Vitaly/Desktop/1task/egoA.csv")
write.csv2(z1, file = "C:/Users/Vitaly/Desktop/1task/altr.csv")
z2 = z1*0.2 + z*0.8
plot(Md_G(seq(-1, 1, length= 101), 1, 100, 80, 0.5, 0))
plot(Md_epsilon(seq(-1, 1, length= 101), 1, 100, 80, 0.5, 0))
plot(seq(-1, 1, length= 101), Md_epsA_cut_Mu_Delta(seq(-1, 1, length= 101), 0.4))
plot(seq(-1, 1, length= 101), Md_altA_cut_Mu_Delta(seq(-1, 1, length= 101), 0.4))



plot(seq(-1, 1, length= 101), Md_epsilon(seq(-1, 1, length= 101), 1, 100, 80, 0.5, 0))
plot(seq(-1, 1, length= 101), Md_epsA_cut_Mu_Delta(seq(-1, 1, length= 101), 0.5))

q = list(1:10,11:20)
write.csv2(q, file = "C:/Users/Vitaly/Desktop/1task/test.csv")


compareEgo <- function(m1,m2, delta) {
  x = seq(m1, m2, length= 101);
  e1 = Md_epsilon(x, 1, 100, 80, 0.5, 0);
  e2 = Md_epsA_cut_Mu_Delta(x, delta);
  q = list(x,e1,e2);
  write.csv2(q, file = "C:/Users/Vitaly/Desktop/1task/test.csv");
}


Phi <- function(z) {
  if (z < -8) return (0);
  if (z >  8) return (1);
  sum = 0;
  term = z;
  i = 3;
  while (sum + term != sum) {
    sum  = sum + term;
    term = term * z * z / i;
    i = i + 2;
  }
  return (0.5 + sum * f(z));
}

Phi1 <- function(zAr) {
  res = NULL;
  for (i in 1:length(zAr)) {
    res <-c(res, Phi(zAr[i]));
  }
  return (res);
}



xMu = seq(-1, 1, length= 101)
yDelta = seq(0.4, 1, length= 100)
z = myOuter(xMu,yDelta, Md_epsA_cut_Mu_Delta)
z1 = myOuter(xMu,yDelta, Md_altA_cut_Mu_Delta)
z2 = z1*0.2 + z*0.8
persp(xMu, yDelta, z2, phi = 25, theta = -55, col = "lightblue",
      xlab = "xMu -1 -> 1", ylab = "yDelta 0.4 -> 1",
      main = "Md_G"
)

write.csv2(z, file = "C:/Users/Vitaly/Desktop/диплом/1task/egoA_alpha05.csv")
write.csv2(z1, file = "C:/Users/Vitaly/Desktop/диплом/1task/altrA_alpha05.csv")
write.csv2(z2, file = "C:/Users/Vitaly/Desktop/диплом/1task/avgA_alpha05.csv")


Sg <- function(mu, sigma, g, t) {
  return (F((mu - t)*sqrt(g)/sigma))
}

Eg <- function(mu, sigma, g, n, alpha) {
  tGammaN = trunc(alpha * n - g);
  p = F(mu/sigma);
  l = n - g;
  return (pbeta(p, tGammaN + 1, l - tGammaN));
}

Rg <- function(mu, sigma, g, n, alpha) {
  tAlphaN = trunc(alpha * n);
  q = F(-mu/sigma);
  l = n - g;
  return (pbeta(q, l - tAlphaN, tAlphaN + 1));
}

Se <- function(mu, sigma) {
  p = F(mu/sigma);
  return (p)
}

Ee <- function(mu, sigma, g, n, alpha, t) {
  Fg = F((mu - t)*sqrt(g)/sigma)
  tGammaN = trunc(alpha * n - g);
  tAlphaN = trunc(alpha * n);
  p = F(mu/sigma);
  q = F(-mu/sigma);
  l = n - g;
  sum = 0;
  for(i in (tAlphaN+1):l) {
    sum = sum + i * choose(l,i) * p^i * q^(l-i)
  }
  
  for(i in (tGammaN+1):(tAlphaN)) {
    sum = sum + i * choose(l,i) * p^i * q^(l-i) * Fg
  }
  
  
  return (sum/(p*l));
}

Re <- function(mu, sigma, g, n, alpha, t) {
  Fng = F((t-mu)*sqrt(g)/sigma)
  tGammaN = trunc(alpha * n - g);
  tAlphaN = trunc(alpha * n);
  p = F(mu/sigma);
  q = F(-mu/sigma);
  l = n - g;
  sum = 0;
  for(i in 1:tGammaN) {
    sum = sum + (l-i) * choose(l,i) * p^i * q^(l-i)
  }
  
  for(i in (tGammaN+1):(tAlphaN)) {
    sum = sum + (l-i) * choose(l,i) * p^i * q^(l-i) * Fng
  }
  
  
  return (sum/(q*l));
}


plot(xMu, Sg(xMu,1, 10, 0))
plot(xMu, Eg(xMu, 1, 10, 50, 0.5))
plot(xMu, Rg(xMu, 1, 10, 50, 0.5))
plot(xMu, Se(xMu,1))
plot(xMu, Ee(xMu, 1, 10, 50, 0.5, 0))
plot(xMu, Re(xMu, 1, 10, 50, 0.5, 0))

x = list(xMu, 
         Sg(xMu,1, 10, 0),
         Eg(xMu, 1, 10, 50, 0.5),
         Rg(xMu, 1, 10, 50, 0.5),
         Se(xMu,1),
         Ee(xMu, 1, 10, 50, 0.5, 0),
         Re(xMu, 1, 10, 50, 0.5, 0)
         )

write.csv2(x, file = "C:/Users/Vitaly/Desktop/диплом/1task/coefficients.csv")

#---------------------------m groups

PG_m_m0_p <- function(m, m0, p)  {
  return(dbinom(m0 - 1, m - 1, p))
}

PGU_m_m0_p <- function(m, m0, p)  {
  result = 0
  for( x in m0:(m-1)) {
    result = result + dbinom(x, m - 1, p)
  }
  return(result)
} 

Md_G_mg <- function(m, g, m0, mu, sigma)  {
  n = m * g;
  P = F(mu*sqrt(g) / sigma)
  f = f(mu*sqrt(g) / sigma)
  result = mu * PGU_m_m0_p(m, m0, P) + (mu * P + sigma * f / sqrt(g)) * PG_m_m0_p(m, m0, P)
  return(result)
}

xMu = seq(-0.7, 0.7, length= 101)
x = list(xMu, 
         Md_G_mg(8, 7, 4, xMu, 1),
         Md_G_mg(8, 7, 5, xMu, 1),
         muplus(xMu, 1, 56, 28),
         muplus(xMu, 1, 56, 29),
         muplus(xMu, 1, 56, 35),
         muplus(xMu, 1, 56, 36)
)
write.csv2(x, file = "C:/Users/Vitaly/Desktop/диплом/1task/clusterSoc.csv")