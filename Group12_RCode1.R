##Signal in different bases
##Loading libraries
library(R1magic)
library(mrbsizeR)
library(wavelets)

##Only 1 plot
par(mfrow = c(1, 1))

##Dimension of Signal
N <- 1024

##Time
t <- 1:N

##Our signal, not sparse
x <-  function(t) {10*sin(20*pi*t/N) - 5*sin(60*pi*t/N) + 4*sin(100*pi*t/N)}
f_x <- x(t)

##Plotting the Signal
plot(x(t), type = "l")

##discrete fourier transform basis
dftbasis <- DFTMatrix0(N)

##Sowing how Fourier transform basis looks like, it's complex
print(DFTMatrix0(4))

##Getting coefficients of the signal in the discrete fourier bases
coef_dft <- t(dftbasis)%*%f_x

##Plot real part of the coefficients, super sparse, only 6 significant coefficients
plot(Re(coef_dft))

##We use identify(Re(real_coef)) to identify which coefficients are significant
##Define Coefficient vector with only 6 nonzero entries
sparse_coef <- rep(0, N)
sparse_coef[c(11, 31, 51, 975, 995, 1015)] <- coef_dft[c(11, 31, 51, 975, 995, 1015)]
plot(x(t), type = "l")

##Perfect reconstruction with only 6 nonzero coefficients
lines(Re(dftbasis%*%-sparse_coef),col="red")
##We do the Same with discrete Cosine Basis
dct_basis <- dctMatrix(N)

## Coef vector
coef_dct <- t(dct_basis)%*%f_x

##not as sparse in dct basis
plot(coef_dct)

##we take only around 5% of the most significant coefficients, use the quantile function to chose this value
q <- quantile(abs(coef_dct), 0.975)
coef_dct[(abs(coef_dct) <= q)] <- 0
##998 out of 1024 coefficients set to 0
sum(coef_dct == 0)

##Reconstruction with 26 highest coefficients
plot(coef_dct)
plot(f_x, type = "l")
##good reconstruction
lines(dct_basis%*%coef_dct, col = "red")



##Now we do everything with Haar bases, this time we dont expect a sparse representation

##Definining Haar matrix

Haar_matrix   <- matrix(c(1, 1, 1, -1), ncol = 2,byrow = "FALSE")
for(i in 1:9){
  Haar_matrix <- cbind(Haar_matrix%x%c(1, 1), diag(2^i)%x%c(1, -1))
}

##Normalize it
for (i in 1:1024){
  Haar_matrix[ ,i ] <- Haar_matrix[ ,i ] / norm(Haar_matrix[ , i], type = "2")
}

##Our Haar Basis
haar_basis <- Haar_matrix

##Coef vector
coef_haar <- t(haar_basis)%*%f_x

##Not as sparse in haar basis
plot(coef_haar)

##Reconstruction with about 10% of the coefficients
qhaar <- quantile(abs(coef_haar), 0.95)
coef_haar[(abs(coef_haar)) <= q] <- 0 

##Around 100 nonzero coefficients (10% of coefficients nonzero)
sum(coef_haar == 0)

##Plot coefficients after setting most to 0
plot(coef_haar)

plot(f_x,type = "l")
##Reconstruction not perfect, not continuous
lines(haar_basis%*%coef_haar, col = "red")

