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
  for(i in x:y) {
    result = result +  b_xnp(i,y,prob)
  }
  return(result)
}

F_gamma <-function(gamma, p, q, l) {
  res = trunc(gamma) + 0.5 - p * l;
  res = res / sqrt(p * q * l);
  return(F(-res))
}

F_gamma_1 <-function(gamma, p, q, l) {
  result = NULL;
  X = 1:length(l)
  #print(X)
  for(i in X) {
    if (gamma > l[i]) {
      result <- c(result, 0)
    } else if (gamma <= 0) {
      result <- c(result, 1)
    } else {
      res = trunc(gamma) + 0.5 - p * l[i];
      res = res / sqrt(p * q * l[i]);
      result <- c(result, F(-res))
    }
  }
  return(result)
}

F_gamma_2 <-function(gamma, p, q, l) {
  result = NULL;
  X = 1:length(l)
  #print(X)
  for(i in X) {
    if (gamma[i] > l[i]) {
      result <- c(result, 0)
    } else if (gamma[i] <= 0) {
      result <- c(result, 1)
    } else {
      res = trunc(gamma[i]) + 0.5 - p * l[i];
      res = res / sqrt(p * q * l[i]);
      result <- c(result, F(-res))
    }
  }
  return(result)
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

muplus_1 <-function(mu, sigma, l, l0) {
  res = NULL;
  X = 1:length(l0)
  for(i in X) {
    if (l0[i] > l) {
      res <- c(res, 0);
    } else if (l0[i] <= 0) {
      res <- c(res, mu);
    } else {
      p = F(mu/sigma);
      q = F(-mu/sigma);
      
      l1 = trunc(l0[i]) + 0.5 - p * l;
      l1 = l1 / sqrt(p * q * l);
      
      fm = f(mu/sigma);
      sfs = sigma * fm / sqrt(p*q*l);
      
      result = mu * F(-l1) + sfs * f(l1);
      res <- c(res, result)
    }
  }
  return(res);
}

muplus_11 <-function(mu, sigma, l, l0) {
  res = NULL;
  X = 1:length(l)
  for(i in X) {
    if (l0 > l[i]) {
      res <- c(res, 0);
    } else if (l0 <= 0) {
      res <- c(res, mu);
    } else {
      p = F(mu/sigma);
      q = F(-mu/sigma);
      
      l1 = trunc(l0) + 0.5 - p * l[i];
      l1 = l1 / sqrt(p * q * l[i]);
      
      fm = f(mu/sigma);
      sfs = sigma * fm / sqrt(p*q*l[i]);
      
      result = mu * F(-l1) + sfs * f(l1);
      res <- c(res, result)
    }
  }
  return(res);
}

muplus_2 <-function(mu, sigma, l, l0) {
  res = NULL;
  X = 1:length(l0)
  for(i in X) {
    if (l0[i] > l[i]) {
      res <- c(res, 0);
    } else if (l0[i] <= 0) {
      res <- c(res, mu);
    } else {
      p = F(mu/sigma);
      q = F(-mu/sigma);
      
      l1 = trunc(l0[i]) + 0.5 - p * l[i];
      l1 = l1 / sqrt(p * q * l[i]);
      
      fm = f(mu/sigma);
      sfs = sigma * fm / sqrt(p*q*l[i]);
      
      result = mu * F(-l1) + sfs * f(l1);
      res <- c(res, result)
    }
  }
  return(res);
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

t0 <- function(mu, sigma, l, n, alpha) {
  delta = l / n;
  gamma = alpha + delta - 1;
  p = F(mu/sigma)
  q = F(-mu/sigma)
  #result = 1
  result = delta / (1- delta);
  result = result * (muplus_11(mu, sigma, l, alpha*n) - muplus_2(mu, sigma, l, gamma*n))
  #result = muplus_11(mu, sigma, l, alpha*n)
  #result = muplus_2(mu, sigma, l, gamma*n)
  result = result / (F_gamma_2(gamma*n, p, q, l) - F_gamma_1(alpha*n, p, q, l))
  #result = (F_gamma_2(gamma*n, p, q, l) - F_gamma_1(alpha*n, p, q, l))
  return (result)
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

Md_epsilon <- function(mu, sigma, n, l, alpha, t) {
  delta = l / n;
  gamma = alpha + delta - 1;
  g = n - l;
  P = F((mu - t)*sqrt(g)/sigma);
  Q = F((t - mu)*sqrt(g)/sigma);
  res = muplus(mu, sigma, l, gamma * n) * P + muplus(mu, sigma, l, alpha * n) * Q;
  return (res);
}

Md_G <- function(mu, sigma, n, l, alpha, t) {
  delta = l / n;
  gamma = alpha + delta - 1;
  g = n - l;
  P = F((mu - t)*sqrt(g)/sigma);
  Q = F((t - mu)*sqrt(g)/sigma);
  p = F(mu/sigma);
  q = F(-mu/sigma);
  f_t = f((mu - t)*sqrt(g)/sigma);
  
  F_gamma_n = F_gamma(gamma * n, p, q, l);
  F_alpha_n = F_gamma(alpha * n, p, q, l);
  
  sfg = sigma * f_t / sqrt(g);
  
  res = F_gamma_n * (mu * P + sfg) + F_alpha_n * (mu * Q - sfg);
  return (res);
}

Md_avg <- function(mu, sigma, n, l, alpha, t) {
  g = n - l;
  res = (Md_G(mu, sigma, n, l, alpha, t) * g + Md_epsilon(mu, sigma, n, l, alpha, t) * l) / n;
  return (res);
}

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


Gsurf_delta_t <- function(delta, t) {
  n = 100;
  mu = -0.1;
  sigma = 1;
  alpha = 0.5;
  l = delta * n;
  return (Md_G(mu, sigma, n, l, alpha, t))
}

Esurf_delta_t <- function(delta, t) {
  n = 100;
  mu = -0.1;
  sigma = 1;
  alpha = 0.5;
  l = delta * n;
  return (Md_epsilon(mu, sigma, n, l, alpha, t))
}

AVGsurf_delta_t <- function(delta, t) {
  n = 100;
  mu = -0.1;
  sigma = 1;
  alpha = 0.5;
  l = delta * n;
  return (Md_avg(mu, sigma, n, l, alpha, t))
}


Gsurf_delta_alpha <- function(delta, alpha) {
  n = 100;
  mu = -0.1;
  sigma = 1;
  t = 0.1;
  l = delta * n;
  return (Md_G(mu, sigma, n, l, alpha, t))
}

Esurf_delta_alpha <- function(delta, alpha) {
  n = 100;
  mu = -0.1;
  sigma = 1;
  t = 0.1;
  l = delta * n;
  return (Md_epsilon(mu, sigma, n, l, alpha, t))
}

AVGsurf_delta_alpha <- function(delta, alpha) {
  n = 100;
  mu = -0.1;
  sigma = 1;
  t = 0.1;
  l = delta * n;
  return (Md_avg(mu, sigma, n, l, alpha, t))
}

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

Md_epsA <- function(mu, sigma, l, n, alpha) {
  a = n -l;
  p = F(mu/sigma);
  res = mu * t(myOuter1(alpha*n, a, Py_x, p)); # Py_x(alpha*n,a,p)
  for(y in 0:alpha*n) {
    res = res + muplus(mu, sigma, l, alpha*n - y) * b_xnp(y, a, p);
  }
  return (res)
}

Md_altA <- function(mu, sigma, l, n, alpha) {
  a = n -l;
  p = F(mu/sigma);
  res = mu * t(myOuter1(alpha*n, a, Py_x, p)); # Py_x(alpha*n,a,p)
  for(y in 0:alpha*n) {
    res = res + muplus(mu, sigma, l, alpha*n - y) * b_xnp(y, a, p);
  }
  return (res)
}

#--------------------------------------------
xMu = seq(-6, 4, length= 101)
x = list(xMu, 
         Md_G(-1, 10, 100, 50, 0.5, xMu),
         Md_epsilon(-1, 10, 100, 50, 0.5, xMu),
         (Md_G(-1, 10, 100, 50, 0.5, xMu) + Md_epsilon(-1, 10, 100, 50, 0.5, xMu)) / 2
)
write.csv2(x, file = "C:/Users/Vitaly/Desktop/диплом/1task/thresholdGr.csv")

#------------------------------------------------------------

x = list((1:99)/100, 
          t0(0.1, 1, 1:99, 100, 0.15),
          t0(0.1, 1, 1:99, 100, 0.46),
          t0(0.1, 1, 1:99, 100, 0.6),
          t0(0.1, 1, 1:99, 100, 0.9)
    )

write.csv2(x, file = "C:/Users/Vitaly/Desktop/диплом/1task/t0Gr.csv")

#-------------------------------------------------------------
t0Surf_delta_alpha <-function(delta, alpha) {
  mu = 0.1
  sigma = 1
  n = 100
  l = delta*n
  return (t0(mu, sigma, l, n, alpha))
}


x = seq(0.15, 0.95, length= 80)
y = seq(0.15, 0.95, length= 80)
z = myOuter(x,y, t0Surf_delta_alpha)
persp(x, y, z, phi = 45, theta = -45, col = "lightblue",
      xlab = "delta 0.01 -> 0.99", ylab = "alpha 0.01 -> 0.99",
      main = "t0"
)

write.csv2(z, file = "C:/Users/Vitaly/Desktop/диплом/1task/t0Surf.csv")
