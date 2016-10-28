create_meth_data <- function(N = 300, pi.c = c(0.45, 0.35, 0.2), max_L = 25,
                            xmin = -100, xmax=100, fmin = -1, fmax = 1){
  set.seed(1)
  # Create a list to store data for each methylation region
  X       <- list()
  # A vector for storing corresponding gene expression data
  Y       <- vector(mode = "numeric", length = N)

  # For each of the N objects
  for (i in 1:N){
    # L is the number of CpGs found in the ith region
    L <- rbinom(n = 1, size = max_L, prob = .8)
    X[[i]] <- matrix(0, nrow = L, ncol = 3)
    # Randomly sample locations for the CpGs
    obs <- sort(sample(xmin:xmax, L))
    # Scale them, so the data lie in the (fmin, fmax) range
    X[[i]][ ,1] <- minmax_scaling(data = obs,
                                  xmin = xmin,
                                  xmax = xmax,
                                  fmin = fmin,
                                  fmax = fmax)

    if (i < N * pi.c[1]){   # First methylation profile
      lb <- round(L / 4)

      X[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X[[i]][1:lb,3] <- rbinom(lb, 14, .9)
        if(all(X[[i]][1:lb,2] > X[[i]][1:lb,3]))
          break
      }

      X[[i]][(lb + 1):L,2] <- rbinom(L - lb, 20, .9)
      repeat{
        X[[i]][(lb + 1):L,3] <- rbinom(L - lb, 2, .9)
        if (all(X[[i]][(lb + 1):L,2] > X[[i]][(lb + 1):L,3]))
          break
      }
      Y[i] <- rpois(1, lambda=200)
    }else if (i < (N * pi.c[2] + N * pi.c[1])){ # Second methylation profile
      lb <- round(L / 1.5)

      X[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X[[i]][1:lb,3] <- rbinom(lb, 2, .8)
        if(all(X[[i]][1:lb,2] > X[[i]][1:lb,3]))
          break
      }

      X[[i]][(lb + 1):L,2] <- rbinom(L - lb, 20, .9)
      repeat{
        X[[i]][(lb + 1):L,3] <- rbinom(L-lb, 14, .9)
        if (all(X[[i]][(lb + 1):L,2] > X[[i]][(lb + 1):L,3]))
          break
      }
      Y[i] <- rpois(1, lambda=100)
    }else{                  # Third methylation profile
      lb <- round(L / 2.5)
      mb <- round(L / 3.5)

      X[[i]][1:lb,2] <- rbinom(lb, 20, .9)
      repeat{
        X[[i]][1:lb,3] <- rbinom(lb, 2, .9)
        if(all(X[[i]][1:lb,2] > X[[i]][1:lb,3]))
          break
      }

      X[[i]][(lb + 1):(lb + mb),2] <- rbinom(mb, 20, .9)
      repeat{
        X[[i]][(lb + 1):(lb + mb),3] <- rbinom(mb, 14, .9)
        if (all(X[[i]][(lb + 1):(lb + mb),2] > X[[i]][(lb + 1):(lb + mb),3]))
          break
      }

      X[[i]][(lb + 1 + mb):L,2] <- rbinom(L - mb - lb, 20, .9)
      repeat{
        X[[i]][(lb + 1 + mb):L,3] <- rbinom(L - mb - lb, 2, .9)
        if (all(X[[i]][(lb + 1 + mb):L,2] > X[[i]][(lb + 1 + mb):L],3))
          break
      }
      Y[i] <- rpois(1, lambda=50)
    }
  }
  return(list(X = X, Y = Y))
}

set.seed(1)
bpr <- create_meth_data(N=600)
meth_data <- bpr$X
gex_data <- bpr$Y
devtools::use_data(meth_data, gex_data, overwrite = TRUE)
