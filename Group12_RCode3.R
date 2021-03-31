##Trying to reconstruct a 32x32 image with compressed sensing using discrete cosine and haar basis
##Load packages
library(linprog)
library(imager)
library(R1magic) 
library(mrbsizeR)
library(wavelets)

##Haarmatrix with 2^10 columns & rows (using a property of the Haar Matrix including kronecker product)
Haar_matrix   <- matrix(c(1, 1, 1, -1), ncol = 2, byrow = "FALSE")
for(i in 1:9){
  Haar_matrix <- cbind(Haar_matrix%x%c(1, 1), diag(2^i)%x%c(1, -1))
}

##Normalize the Haar Matrix
for (i in 1:1024){
  Haar_matrix[ , i] <- Haar_matrix[ , i] / norm(Haar_matrix[ , i], type = "2")
}

##Haar Basis and Discrete Cosine Basis
haar_basis <- Haar_matrix
dct_basis  <- dctMatrix(1024)

##32x32 Image
im                <- imfill(32, 32)
##creating image that I want a sparse recovery for
for(i in 1:32){
  im[ , i]        <- cos(i) / 8
}
im[3:4, 8:24]     <- 1
im[3:10, 8:9]     <- 1
im[3:10, 16:17]   <- 1
im[3:10, 24:25]   <- 1
im[12:20, 8:9]    <- 1
im[15:16, 8:25]   <- 1
im[22:23, 8:25]   <- 1
im[22:29, 16:17]  <- 1
im[29:30, 8:25]   <- 1
##Plotting image
par(mfrow = c(1, 1))
plot(im)

##Converting 32x32 Image into a 1024 dimensional vecot
vector_image <- as.vector(im)
str(vector_image)

##Set seed for reconstruction
set.seed(1)

##Defining dimensions, want to reduce 1024 dimensional signal to M dimensional vector 
##And then recover it in a sparse way
##N: dimension of signal
##Different M values: Dimensions, we project our data to
##First we look at the influence of randomness by chosing the same M twice
M                   <- c(250, 250, 750, 750)
##Then we take different M values
M                   <- c(250, 550, 750, 1000)
N                   <- 1024

##Defining Coefficient vectors, before for loop, for discrete cosine transform and haar transform
coef_haar_norm      <- matrix(NA, nrow = N, ncol = length(M))
coef_dct_norm       <- matrix(NA, nrow = N, ncol = length(M))
coef_haar_bin       <- matrix(NA, nrow = N, ncol = length(M))
coef_dct_bin        <- matrix(NA, nrow = N, ncol = length(M))


##Defining Recovered Signal as NA matrix, before for loop
newsignal_haar_norm <- matrix(NA, nrow = N, ncol = length(M))
newsignal_dct_norm  <- matrix(NA, nrow = N, ncol = length(M))
newsignal_haar_bin  <- matrix(NA, nrow = N, ncol = length(M))
newsignal_dct_bin   <- matrix(NA, nrow = N, ncol = length(M))

##Doing everything (Reconstruction of image) for all values of M, 
##For Gaussian and Binomial random matrix and for Haar and Cosine basis
for(i in 1:length(M)){

  ##Random matrix to construct projection y , both chosen to have column norm 1
  phi_norm <- matrix(data = rnorm(N * M[i], sd = 1 / sqrt(M[i])), nrow = M[i], ncol = N)
  phi_bin  <- matrix((2 * rbinom(M[i] * N, 1, 0.5) - 1) / sqrt(M[i]), nrow = M[i], ncol = N)
  
  ##Random Projection vector y, y has  dimension M[i] (250,550,750,1000), as compared to 1024 of original image
  y_norm   <- phi_norm%*%vector_image
  y_bin    <- phi_bin%*%vector_image
  
  ##Minimize b such that y=Phi%*%dftbasis%*%b   (solve linear program)
  ##Using the function lp
  ##Important for LP function: Note that every variable is assumed to be >= 0!
  ## -> We have to split in in positive and negative part

  program_dct_norm  <- lp(direction = "min", objective.in = rep(1,2*N), 
                          const.mat = cbind(phi_norm%*%dct_basis, -phi_norm%*%dct_basis), 
                          const.dir = rep("==", M[i]), const.rhs = y_norm)
  program_haar_norm <- lp(direction = "min", objective.in = rep(1, 2*N), 
                          const.mat = cbind(phi_norm%*%haar_basis, -phi_norm%*%haar_basis),
                          const.dir = rep("==", M[i]), const.rhs = y_norm)
  program_dct_bin   <- lp(direction = "min", objective.in = rep(1, 2*N), 
                          const.mat = cbind(phi_bin%*%dct_basis,  -phi_bin%*%dct_basis), 
                          const.dir = rep("==", M[i]), const.rhs = y_bin)
  program_haar_bin  <- lp(direction = "min", objective.in = rep(1, 2*N), 
                          const.mat = cbind(phi_bin%*%haar_basis, -phi_bin%*%haar_basis),
                          const.dir = rep("==", M[i]), const.rhs = y_bin)

  
  ##Coefficients have Constraint to be > 0 so we get 2*N dimensional vectors here we construct the
  ##Coefficient vector from positive and negative part
  coef_dct_norm[ , i]  <- program_dct_norm$solution[1:N]  - program_dct_norm$solution[(N+1):(2*N)]
  coef_haar_norm[ , i] <- program_haar_norm$solution[1:N] - program_haar_norm$solution[(N+1):(2*N)]
  coef_dct_bin[ , i]   <- program_dct_bin$solution[1:N]   - program_dct_bin$solution[(N+1):(2*N)]
  coef_haar_bin[ , i]  <- program_haar_bin$solution[1:N]  - program_haar_bin$solution[(N+1):(2*N)]

  ##We reconstruct the signal
  newsignal_haar_norm[ , i] <- haar_basis%*%coef_haar_norm[ , i]
  newsignal_dct_norm[ , i]  <- dct_basis%*%coef_dct_norm[ , i]
  newsignal_haar_bin[ , i]  <- haar_basis%*%coef_haar_bin[ , i]
  newsignal_dct_bin[ , i]   <- dct_basis%*%coef_dct_bin[ , i]

}

##Now we do some plots to compare different bases/random matrices
par(mfrow = c(4, 3))

##Compare Haar and DCT Basis for Gaussian Matrix
for (i in 1:length(M)){
  plot(im,main = "Original image")
  plot(as.cimg(as.vector(newsignal_haar_norm[ , i])), main = (paste("Haar & Normal, M=", M[i])))
  plot(as.cimg(as.vector(newsignal_dct_norm[ , i])),  main = (paste("Cosine & Normal, M=", M[i])))
}

&##Compare Haar and DCT Basis for Binomial Matrix
for (i in 1:length(M)){
  plot(im,main = "Original Image")
  plot(as.cimg(as.vector(newsignal_haar_bin[ , i])),  main = (paste("Haar & Binomial, M=", M[i])))
  plot(as.cimg(as.vector(newsignal_dct_bin[ , i])),   main = (paste("Cosine & Binomial, M=", M[i])))
}

##Compare Normal and Binomial Matrix for Haar Basis
for (i in 1:length(M)){
  plot(im, main = "Original Image")
  plot(as.cimg(as.vector(newsignal_haar_norm[ , i])), main = (paste("Haar & Normal, M=", M[i])))
  plot(as.cimg(as.vector(newsignal_haar_bin[ , i])),  main = (paste("Haar & Binomial, M=", M[i])))
}

##Compare Normal and Binomial Matrix for DCT Basis
for (i in 1:length(M)){
  plot(im,main = "Original Image")
  plot(as.cimg(as.vector(newsignal_dct_norm[ , i])),  main = (paste("Cosine & Normal, M=", M[i])))
  plot(as.cimg(as.vector(newsignal_dct_bin[ , i])),   main = (paste("Cosine & Binomial, M=", M[i])))
}

##Now we take a look at the coefficient vectors for the Gaussian random matrix
par(mfrow = c(4, 2))
for (i in 1:length(M)){
  plot(coef_haar_norm[ , i], ylab = "Haar Coefficient",  main = paste("Haar and Normal, M=", M[i]))
  plot(coef_dct_norm[ , i], ylab = "Cosine Coefficient", main = paste("Cosine and Normal, M=", M[i]))
}


##Now we take the transposed image, and we will see that the image is not sparse anymore for the haar basis
##That's because R goes over the image column wise and now our columns have this cosine function in it
##And are not constant anymore
im_transposed <- imfill(32,32)
im_transposed[1:32, 1:32, 1, 1] <- t(im[1:32, 1:32, 1, 1])
par(mfrow=c(1, 1))
plot(im_transposed)
vector_image_t <- as.vector(im_transposed)

##We define coefficient vectors, this time we do it only for normal random matrix
coef_haar_norm_t      <- matrix(NA, nrow = N, ncol = length(M))
coef_dct_norm_t       <- matrix(NA, nrow = N, ncol = length(M))

##Defining Recovered Signal as NA matrix
newsignal_haar_norm_t <- matrix(NA, nrow = N, ncol = length(M))
newsignal_dct_norm_t  <- matrix(NA, nrow = N, ncol = length(M))

##Here we have almost the same for loop as before
for(i in 1:length(M)){
  
  ##Random matrix to construct projection y 
  phi_norm <- matrix(data = rnorm(N * M[i], sd = 1 / sqrt(M[i])), nrow = M[i], ncol = N)
  
  ##Random Projection vector y, y has now dimension M[i] (250,550,750,1000), as compared to 1024 from original image
  y_norm   <- phi_norm%*%vector_image_t
  
  ##Solving the linear program with lp 
  program_dct_norm_t  <- lp(direction = "min", objective.in = rep(1, 2*N), 
                          const.mat = cbind(phi_norm%*%dct_basis, -phi_norm%*%dct_basis), 
                          const.dir = rep("==", M[i]), const.rhs = y_norm)
  program_haar_norm_t <- lp(direction = "min", objective.in = rep(1, 2*N), 
                          const.mat = cbind(phi_norm%*%haar_basis, -phi_norm%*%haar_basis),
                          const.dir = rep("==", M[i]), const.rhs = y_norm)
  
  ##Constructing coefficient vector
  coef_dct_norm_t[ , i]  <- program_dct_norm_t$solution[1:N]  - program_dct_norm_t$solution[(N+1):(2*N)]
  coef_haar_norm_t[ , i] <- program_haar_norm_t$solution[1:N] - program_haar_norm_t$solution[(N+1):(2*N)]
  
  ##We reconstruct the signal
  newsignal_haar_norm_t[ , i] <- haar_basis%*%coef_haar_norm_t[ , i]
  newsignal_dct_norm_t[ , i]  <- dct_basis%*%coef_dct_norm_t[ , i]
}
par(mfrow=c(4,3))

##That was the original image, we see we get a sparse solution for the haar basis
for (i in 1:length(M)){
  plot(im, main = "Original Image")
  plot(as.cimg(as.vector(newsignal_haar_norm[ , i])), main = paste("Haar & Normal, M=", M[i]))
  plot(as.cimg(as.vector(newsignal_dct_norm[ , i])),  main = paste("Cosine & Normal,M=", M[i]))
}

##Now we do the same for the transposed image and we see that we do not get a sparse solution anymore
for (i in 1:length(M)){
  plot(im_transposed, main = "Transposed Image")
  plot(as.cimg(as.vector(newsignal_haar_norm_t[ , i])), main = paste("Haar & Normal, M=", M[i]))
  plot(as.cimg(as.vector(newsignal_dct_norm_t[ , i])),  main = paste("Cosine & Normal, M=", M[i]))
}


##Now we do a big simulation for a lot of M values and trying everything also for a uniform random matrix

##Haarmatrix with 2^10 columns & rows
Haar_matrix_long   <- matrix(c(1, 1, 1, -1), ncol = 2, byrow = "FALSE")
for(i in 1:9){
  Haar_matrix_long <- cbind(Haar_matrix_long%x%c(1 , 1), diag(2^i)%x%c(1, -1))
}

##Normalize the Haar Matrix
for (i in 1:1024){
  Haar_matrix_long[ , i] <- Haar_matrix_long[ , i] / norm(Haar_matrix_long[ , i],type = "2")
}

##Haar Basis and Discrete Cosine Basis
haar_basis_long <- Haar_matrix_long
dct_basis_long  <- dctMatrix(1024)


##32x32 Image
im_long <- imfill(32, 32)
##Creating image that i want a sparse recovery for
for(i in 1:32){
  im_long[ , i]       <- cos(i) / 8
}
im_long[3:4, 8:24]    <- 1
im_long[3:10, 8:9]    <- 1
im_long[3:10, 16:17]  <- 1
im_long[3:10, 24:25]  <- 1
im_long[12:20, 8:9]   <- 1
im_long[15:16, 8:25]  <- 1
im_long[22:23, 8:25]  <- 1
im_long[22:29, 16:17] <- 1
im_long[29:30, 8:25]  <- 1
par(mfrow=c(1, 1))
plot(im_long)

##Converting 32x32 Image into a 1024 dimensional vecot
vector_image_long <- as.vector(im_long)
str(vector_image_long)

##Set seed for reconstruction
set.seed(7)

##Defining dimensions, want to reduce 1024 dimensional signal to M dimensional vector and then recover it 
##in a different basis
##N: dimension of signal
##Different M values: Dimensions, we project or data to, we chose each value twice to see effect of randomness
M_long <- c(150, 150 ,200 ,200 ,250 ,250 ,300 ,300, 350, 350, 400, 400, 450, 450, 500, 500, 550, 550, 600, 600,
            650, 650, 700, 700, 750, 750, 800, 800, 850, 850, 900, 900, 950, 950, 1000, 1000, 1024, 1024)
N      <- 1024

##Defining Coefficient vectors, before for loop, for discrete cosine transform and haar transform
coef_haar_norm_long      <- matrix(NA, nrow = N, ncol = length(M_long))
coef_dct_norm_long       <- matrix(NA, nrow = N, ncol = length(M_long))
coef_haar_bin_long       <- matrix(NA, nrow = N, ncol = length(M_long))
coef_dct_bin_long        <- matrix(NA, nrow = N, ncol = length(M_long))
coef_haar_uni_long       <- matrix(NA, nrow = N, ncol = length(M_long))
coef_dct_uni_long        <- matrix(NA, nrow = N, ncol = length(M_long))

##Defining Recovered Signal as NA matrix, before for loop
newsignal_haar_norm_long <- matrix(NA, nrow = N, ncol = length(M_long))
newsignal_dct_norm_long  <- matrix(NA, nrow = N, ncol = length(M_long))
newsignal_haar_bin_long  <- matrix(NA, nrow = N, ncol = length(M_long))
newsignal_dct_bin_long   <- matrix(NA, nrow = N, ncol = length(M_long))
newsignal_haar_uni_long  <- matrix(NA, nrow = N, ncol = length(M_long))
newsignal_dct_uni_long   <- matrix(NA, nrow = N, ncol = length(M_long))

##Doing everything for all values of M, for Gaussian, Binomial and this time also uniform random matrix
##This can take 1-2 hours!
for(i in 1:length(M_long)){
  
  ##Random matrix to construct projection y 
  phi_norm_long <- matrix(data = rnorm(N * M_long[i], sd=1/sqrt(M_long[i])), nrow = M_long[i], ncol = N)
  phi_bin_long  <- matrix((2 * rbinom(M_long[i] * N, 1, 0.5) - 1) / (sqrt(M_long[i])), nrow= M_long[i], ncol = N)
  ##Defined Uniform matrix different than I showd in presentation:
  phi_uni_long  <- matrix(data = runif(N * M_long[i], min = -  sqrt(12) / (2*sqrt(M_long[i])), max = sqrt(12) / (2*sqrt(M_long[i]))), nrow = M_long[i], ncol = N)
  
  ##In the presentation I used an other matrix and said that we somehow get worse results in the uniform case
  ##Now with the new matrix it works as expected and we get a good reconstruction with the uniform matrix :)
  ##This was how I defined the uniform matrix before:
  ##phi_uni_long  <- matrix(data = runif(N * M_long[i], min = 0, max = 1), nrow = M_long[i], ncol = N)
  ##for(k in 1:M_long[i]){
    ##phi_uni_long[,k] <- phi_uni_long[, k] / norm(phi_uni_long[ , k], type = "2")
  ##}
  
  
  ##Random Projection vector y
  y_bin_long  <- phi_bin_long%*%vector_image_long
  y_norm_long <- phi_norm_long%*%vector_image_long
  y_uni_long  <- phi_uni_long%*%vector_image_long
  
  ##Solving linear program with lp
  program_dct_norm_long  <- lp(direction = "min", objective.in=rep(1, 2*N), 
                               const.mat = cbind(phi_norm_long%*%dct_basis_long, -phi_norm_long%*%dct_basis_long), 
                               const.dir = rep("==", M_long[i]), const.rhs = y_norm_long)
  program_haar_norm_long <- lp(direction = "min", objective.in = rep(1, 2*N), 
                               const.mat = cbind(phi_norm_long%*%haar_basis_long, -phi_norm_long%*%haar_basis_long),
                               const.dir = rep("==", M_long[i]), const.rhs = y_norm_long)
  program_dct_bin_long   <- lp(direction = "min", objective.in = rep(1, 2*N), 
                               const.mat = cbind(phi_bin_long%*%dct_basis_long, -phi_bin_long%*%dct_basis_long), 
                               const.dir = rep("==", M_long[i]), const.rhs=y_bin_long)
  program_haar_bin_long  <- lp(direction="min", objective.in=rep(1,2*N), 
                               const.mat = cbind(phi_bin_long%*%haar_basis_long, -phi_bin_long%*%haar_basis_long),
                               const.dir = rep("==", M_long[i]), const.rhs = y_bin_long)
  program_dct_uni_long   <- lp(direction = "min", objective.in = rep(1, 2*N), 
                               const.mat = cbind(phi_uni_long%*%dct_basis_long, -phi_uni_long%*%dct_basis_long), 
                               const.dir = rep("==", M_long[i]), const.rhs = y_uni_long)
  program_haar_uni_long  <- lp(direction = "min", objective.in = rep(1, 2*N), 
                               const.mat = cbind(phi_uni_long%*%haar_basis_long, -phi_uni_long%*%haar_basis_long),
                               const.dir = rep("==", M_long[i]), const.rhs = y_uni_long)
  
  
  ##We construct coefficients vector
  coef_dct_norm_long[ , i]  <- program_dct_norm_long$solution[1:N]  - program_dct_norm_long$solution[(N+1):(2*N)]
  coef_haar_norm_long[ , i] <- program_haar_norm_long$solution[1:N] - program_haar_norm_long$solution[(N+1):(2*N)]
  coef_dct_bin_long[ , i]   <- program_dct_bin_long$solution[1:N]   - program_dct_bin_long$solution[(N+1):(2*N)]
  coef_haar_bin_long[ , i]  <- program_haar_bin_long$solution[1:N]  - program_haar_bin_long$solution[(N+1):(2*N)]
  coef_dct_uni_long[ , i]   <- program_dct_uni_long$solution[1:N]   - program_dct_uni_long$solution[(N+1):(2*N)]
  coef_haar_uni_long[ , i]  <- program_haar_uni_long$solution[1:N]  - program_haar_uni_long$solution[(N+1):(2*N)]
  
  ##We reconstruct the signal
  newsignal_haar_norm_long[ , i] <- haar_basis_long%*%coef_haar_norm_long[ , i]
  newsignal_dct_norm_long[ , i]  <- dct_basis_long%*%coef_dct_norm_long[ , i]
  newsignal_haar_bin_long[ , i]  <- haar_basis_long%*%coef_haar_bin_long[ , i]
  newsignal_dct_bin_long[ , i]   <- dct_basis_long%*%coef_dct_bin_long[ , i]
  newsignal_haar_uni_long[ , i]  <- haar_basis_long%*%coef_haar_uni_long[ , i]
  newsignal_dct_uni_long[ , i]   <- dct_basis_long%*%coef_dct_uni_long[ , i]
}


##We look at the mean error of reconstruction vs original image, pixel wise
mean_error_dct_norm_long  <- rep(NA, length(M_long))
mean_error_haar_norm_long <- rep(NA, length(M_long))
mean_error_dct_bin_long   <- rep(NA, length(M_long))
mean_error_haar_bin_long  <- rep(NA, length(M_long))
mean_error_dct_uni_long   <- rep(NA, length(M_long))
mean_error_haar_uni_long  <- rep(NA, length(M_long))

for(i in 1:length(M_long)){
  mean_error_haar_norm_long[i]   <- mean(abs(newsignal_haar_norm_long[ , i] - vector_image_long))
  mean_error_dct_norm_long[i]    <- mean(abs(newsignal_dct_norm_long[ , i]  - vector_image_long))
  mean_error_haar_bin_long[i]    <- mean(abs(newsignal_haar_bin_long[ , i]  - vector_image_long))
  mean_error_dct_bin_long[i]     <- mean(abs(newsignal_dct_bin_long[ , i]   - vector_image_long))
  mean_error_haar_uni_long[i]    <- mean(abs(newsignal_haar_uni_long[ , i]  - vector_image_long))
  mean_error_dct_uni_long[i]     <- mean(abs(newsignal_dct_uni_long[ , i]   - vector_image_long))
}

##Now we look at a plot of the mean errors
##First we do it for only the different M values, below we also do it for to times the same M,
##Where we get a lot of lines, to see the influence of randomness
par(mfrow = c(1, 1))
plot(M_long[seq(1, 38, by = 2)], mean_error_dct_norm_long[seq(1 , 38, by = 2)], main = "Plot of mean Errors", 
     ylab = "Mean Error", xlab = "M", type = "l", lwd = "2")
legend("topright", col = c("black", "red", "green"), lty = c(1, 2), 
       legend = c( "Normal & Cosine", "Binomial & Haar", "Uniform & Cosine", "Normal & Haar", "Binomial & Cosine", "Uniform & Haar"))

lines(M_long[seq(2, 38, by = 2)], mean_error_haar_norm_long[seq(2, 38, by = 2)], lty = 2, lwd = "2")
lines(M_long[seq(2, 38, by = 2)], mean_error_dct_bin_long[seq(2, 38, by = 2)], col = "red", lwd = "2")
lines(M_long[seq(2, 38, by = 2)], mean_error_haar_bin_long[seq(2, 38, by = 2)], lty = 2, col = "red", lwd = "2")
lines(M_long[seq(2, 38, by = 2)], mean_error_dct_uni_long[seq(2, 38, by = 2)], col = "green", lwd = "2")
lines(M_long[seq(2, 38, by = 2)], mean_error_haar_uni_long[seq(2, 38, by = 2)], lty = 2, col = "green", lwd = "2")

##Plot with all lines, we take each M twice
plot(M_long[seq(1, 38, by = 2)], mean_error_dct_norm_long[seq(1 , 38, by = 2)], main = "Plot of mean Errors", 
     ylab = "Mean Error", xlab = "M", type = "l", lwd = "2")
legend("topright", col = c("black", "red", "green"), lty = c(1, 2), 
       legend = c( "Normal & Cosine", "Binomial & Haar", "Uniform & Cosine", "Normal & Haar", "Binomial & Cosine", "Uniform & Haar"))
lines(M_long[seq(2, 38, by = 2)], mean_error_dct_norm_long[seq(2, 38, by = 2)], lwd = "2")

##Normal & Haar
lines(M_long[seq(2, 38, by = 2)], mean_error_haar_norm_long[seq(2, 38, by = 2)], lty = 2, lwd = "2")
lines(M_long[seq(1, 38, by = 2)], mean_error_haar_norm_long[seq(1, 38, by = 2)], lty = 2, lwd = "2")

##Binomial & Cosine
lines(M_long[seq(1, 38, by = 2)], mean_error_dct_bin_long[seq(1, 38, by = 2)], col = "red", lwd = "2")
lines(M_long[seq(2, 38, by = 2)], mean_error_dct_bin_long[seq(2, 38, by = 2)], col = "red", lwd = "2")

##Binomial & Haar
lines(M_long[seq(1, 38, by = 2)], mean_error_haar_bin_long[seq(1, 38, by = 2)], lty = 2, col = "red", lwd = "2")
lines(M_long[seq(2, 38, by = 2)], mean_error_haar_bin_long[seq(2, 38, by = 2)], lty = 2, col = "red", lwd = "2")

##Uniform & Cosine
lines(M_long[seq(1, 38, by = 2)], mean_error_dct_uni_long[seq(1, 38, by = 2)], col = "green", lwd = "2")
lines(M_long[seq(2, 38, by = 2)], mean_error_dct_uni_long[seq(2, 38, by = 2)], col = "green", lwd = "2")

##Uniform & Haar
lines(M_long[seq(1, 38, by = 2)], mean_error_haar_uni_long[seq(1, 38, by = 2)], lty = 2, col = "green", lwd = "2")
lines(M_long[seq(2, 38, by = 2)], mean_error_haar_uni_long[seq(2, 38, by = 2)], lty = 2, col = "green", lwd = "2")

