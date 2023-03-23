## Set up the packages and input the data

library(deSolve)
library(ggpubr)
library(ggplot2)
library(readxl)
library(readr)
library(writexl)

setwd("/Users/u2091432/Documents/MIBTP/African horse sickness")

# dataset used for horse population numbers, vaccination % and vaccination protection factor
AHS_dataset <- read_excel("large_grid_data.xlsx")

new_AHS_dataset <- data.frame()
for (row in 1:nrow(AHS_dataset)){
  if (AHS_dataset[row,2] >1){
    new_AHS_dataset <- rbind(new_AHS_dataset, AHS_dataset[row,])
  }
}

AHS_dataset <- new_AHS_dataset

# dataset derived from IDW interpolation from each month of c. imicola counts (from K.Labuschagne), then a periodic gaussian function fitted to the
# data to enable an estimation of daily midge counts throughout the year. Fitting was done with least sum of squares, and the best
# fit was when the maximum midge count was fixed to be the same as the maximum midges of the data. This CSV file details the daily midge estimation.
new_midge_df <- read_csv("midge_df_gaus_large_grid.csv")
new_midge_df <- as.data.frame(new_midge_df)


# satellite temperatures (from WorldClim) were downloaded as a raster layer and mean temperature extracted for all cells in each month. Similar
# to the midge fitting, a sine wave was fitted using least sum of squares so that daily temperature could be established. This CSV file details
# the inputs for the sine wave equation for each cell.
temp_data <- read_csv("temp_sine_wave_df_large_grid.csv")

## next create the functions for temperatures and midges for each grid cell

# time for simulation is 360 days but the waves need to repeat for 720 days so that the simulation
# that starts in month 12 has the full run of the wave

t <- seq(1,720,1)

# create a list for the temperature functions for each cell
all_temp_waves <- function(times){return(temp_data[[cellrow,2]]*sin((times/temp_data[[cellrow,3]])+temp_data[[cellrow,4]])+temp_data[[cellrow,5]])}

# create a list for the midge functions for each cell
midge_df <- vector("list", 153)
for (j in 1:153){
  midge_fun <- splinefun(t, predict(loess(new_midge_df[1:720,(j+1)] ~ t, data = new_midge_df, span = 0.05)), method = "periodic")
  midge_df[[j]] <- midge_fun
}

## Determine which grid cells are neighbours of each other

# Define the maximum number of grid squares.
n_x <- 15
n_y <- 16
n_max <- n_x * n_y

# Define which grid squares k are neighbours of gridsquares l.
gs_weights <- matrix(0, nrow=n_max, ncol=n_max)
for (k in seq(1, n_max))
{
  # Given the grid square i, what are its co-ordinates
  # in the overall grid (using matrix notation).
  j <- (k-1) %% n_x + 1
  i <- (k-1) %/% n_x + 1
  
  # Create a vector to store the neighbours of (i, j)
  neighb <- c()
  
  # If the grid is not a row-vector...
  if (n_y > 1) {
    
    # Determine the grid squares that can exist above
    # and below grid square k.
    if (i == 1) {
      neighb <- c(neighb, i*n_x + j)
    } else if (i == n_y) {
      neighb <- c(neighb, (i-2)*n_x + j)
    } else {
      neighb <- c(neighb, (i-2)*n_x + j, i*n_x + j)
    }
  }
  
  # If the grid square is not a column-vector...
  if (n_x > 1)
  {
    # Determine the grid squares that are to the left
    # and to the right of grid square k.
    if (j == 1) {
      neighb <- c(neighb, (i-1)*n_x + j + 1)
    } else if (j == n_x) {
      neighb <- c(neighb, (i-1)*n_x + j - 1)
    } else {
      neighb <- c(neighb, (i-1)*n_x + j + 1, (i-1)*n_x + j - 1)
    }
  }
  
  # Set neighbour grid squares to 1 in the matrix.
  gs_weights[k, neighb] <- 1
}

# create list of cell ID numbers and assign IDs

ID <- AHS_dataset$id
column_df <- data.frame(rep(0,240))

for (i in 1:length(ID)){
  col <- ID[i]
  column_df <- cbind(column_df, gs_weights[,col])
}

column_df <- column_df[,-1]

neighbours_df <- data.frame(matrix(nrow = 1, rep(0, 153)))
names(neighbours_df) <- ID

for (i in 1:length(ID)){
  row <- ID[i]
  new_row <- column_df[row,]
  names(new_row) <- ID
  neighbours_df <- rbind(neighbours_df, new_row)
}

neighbours_df <- neighbours_df[-1,]

colnames(neighbours_df) <- ID
rownames(neighbours_df) <- ID

neighbours_list <- list()

for (x in 1:nrow(neighbours_df)){
  neighbours_vec <- c()
  for (y in 1:ncol(neighbours_df)){
    if (neighbours_df[x,y] == 1){
      neighbours_vec <- c(neighbours_vec, colnames(neighbours_df)[y])
      neighbours_vec <- as.numeric(neighbours_vec)
    }
    neighbours_list[[x]] <- neighbours_vec
  }
}

# Set up root and event functions so that when the latent and infected midges and
# latent and infected horses are cumulatively under 1, the root function triggers the event function.
# In this function, a root is found at =<0, whereas otherwise it returns 1, and the event function is
# not triggered.

rootfun <- function(times,x,pars){
  
  l <- pars$l #no. of stages for latent host class
  n1 <- pars$n1 #no. of stages for 1st infectious host class
  n2 <- pars$n2 #no. of stages for 2nd infectious host class
  k <- pars$k#no. of stages for incubating vector class
  starttime <- pars$starttime #day on which the outbreak starts - each month is 30 days
  n <- pars$n #number of cells/patches in the model
  
  #initial conditions
  Sh <- x[1:n]
  Eh <- x[(n+1) : (6*n)]
  I1h <- x[((6*n)+1) : (17*n)]
  I2h <- x[((17*n)+1) : (30*n)]
  Rh <- x[((30*n) +1) : (31*n)]
  Dh <- x[((31*n) +1) : (32*n)]
  Sv <- x[((32*n) +1) : (33*n)]
  Ev <- x[((33*n) +1): (43*n)]
  Iv <- x[((43*n) +1) : (44*n)]
  Vh <- x[((44*n) +1) : (45*n)]
  
  #cut off time set to average EIP + horse infectious period + horse incubation period 10.9+8.7+4.6 (Emma paper)
  if (times>(starttime+24)){
    return((ceiling(sum(Eh))-1) + (ceiling(sum(I1h))-1) + (ceiling(sum(I2h))-1) +
             (ceiling(sum(Ev))-1) + (ceiling(Iv)-1))
  }else{
    return(1)
  }
  
}


## Event function is triggered when the root function returns a root. The event is to assign 0 to the latent/infected
## horses/midges compartments, which in turn stops the outbreak from continuing on and prevents a delayed (unrealistic) wave.

eventfun <- function(times,x,pars){
  
  l <- pars$l #no. of stages for latent host class
  n1 <- pars$n1 #no. of stages for 1st infectious host class
  n2 <- pars$n2 #no. of stages for 2nd infectious host class
  k <- pars$k #no. of stages for incubating vector class
  n <- pars$n #no. of cells/patches
  
  Sh <- x[1:n]
  Eh <- x[(n+1) : (6*n)]
  I1h <- x[((6*n)+1) : (17*n)]
  I2h <- x[((17*n)+1) : (30*n)]
  Rh <- x[((30*n) +1) : (31*n)]
  Dh <- x[((31*n) +1) : (32*n)]
  Sv <- x[((32*n) +1) : (33*n)]
  Ev <- x[((33*n) +1): (43*n)]
  Iv <- x[((43*n) +1) : (44*n)]
  Vh <- x[((44*n) +1) : (45*n)]
  
  # Indices which require being set to 0
  indices <- c((n+1):(30*n), ((33*n)+1):(44*n))
  
  newx <- x
  newx[indices]<-rep(0, 40*n)
  
  
  return(x <- newx)
}

#########################################################################

## Define required empty data frames and lists

results_df <- data.frame()
summary_graph_list <- list()
outbreak_graph_list <- list()
midges_graph_list <- list()
summary_combined <- data.frame()

## Define how many grid cells the model must iterate over

max_cells <- 153

## Open a for-loop for the grid cell iteration

for (cell in 44){
  
  # Define the cell ID and row numbers of the neighbouring cells
  
  neighboursID <- neighbours_list[[cell]]
  neighbours <- c(cell)
  for (rep in 1:length(neighboursID)){
    row <- which(AHS_dataset$id == neighboursID[rep], arr.ind=TRUE)
    neighbours <- c(neighbours,row)} #vector containing the data frame row numbers of all the cells, starting with the central cell
  
  
  #Define the model's function
  
  AHS.coupling.model <- function(times, x, pars){ 
    
    # find the number of patches
    n <- pars$n
    
    # find the parameters
    epsilon <- pars$epsilon #latent period for hosts (days) = 1/epsilon
    l <- pars$l #no. of stages for latent host class
    gamma1 <- pars$gamma1 #overall recovery rate from 1st infectious host class
    n1 <- pars$n1 #no. of stages for 1st infectious host class
    gamma2 <- pars$gamma2 #overall recovery rate from 2nd infectious host class
    n2 <- pars$n2 #no. of stages for 2nd infectious host class
    mh <- pars$mh #mortality hosts
    k <- pars$k #no. of stages for incubating vector class
    ph <- pars$ph #transmission probability host to vector
    pv <- pars$pv #transmission probability vector to host
    pf <- pars$pf #vaccine protective factor
    r <- pars$r  # force of infection between neighbouring cells
    
    
    # find the current states for each compartment in each patch
    Sh <- x[1:n]
    Eh <- x[(n+1) : (6*n)]
    I1h <- x[((6*n)+1) : (17*n)]
    I2h <- x[((17*n)+1) : (30*n)]
    Rh <- x[((30*n) +1) : (31*n)] 
    Dh <- x[((31*n) +1) : (32*n)] 
    Sv <- x[((32*n) +1) : (33*n)] 
    Ev <- x[((33*n) +1): (43*n)] 
    Iv <- x[((43*n) +1) : (44*n)] 
    Vh <- x[((44*n) +1) : (45*n)] 
    
    # initialise vectors to store the derivatives
    dSh <- c()
    dEh <- c()
    dI1h <- c()
    dI2h <- c()
    dRh <- c()
    dDh <- c()
    dSv <- c()
    dEv <- c()
    dIv <- c()
    dVh <- c()
    
    #iterate over all the cells/patches in the model, starting with the central cell, and working on to the neighbours
    
    for (patch in 1:n){
      
      # determine which row in the list of functions the patch correlates to
      cellrow <- neighbours[patch]
      
      # determine the temperature wave using every time step
      
      temp <-  all_temp_waves(times)
      
      # calculate the temperature dependent variable vectors for every time step
      
      alpha <- (0.015*temp)-0.125        
      upsilon <- (0.0085*temp) - 0.0821 
      muv <- 0.015*exp(0.063*temp)
      
      # define the derivative of the midge wave so that the change in midges can be added at each time step
      
      omega <- midge_df[[cellrow]](times, deriv = 1)
      
      # calculate the derivatives, where there is a force of infection between neighbours. 
      # This assumes an infected midge could cross the border and bite a susceptible/vaccinated horse.
      # Equally a susceptible midge could cross the border and bite an infected horse
      # No movement of horses as assumed immediate travel lockdown when a case occurs.
      
      dSh[patch] <- -sum(((alpha*pv) * ((Iv * Sh[patch]) / (Sh[patch] + sum(Eh[(l*(patch-1)+1):(l*patch)]) +
                    sum(I1h[(n1*(patch-1)+1):(n1*patch)]) + sum(I2h[(n2*(patch-1)+1):(n2*patch)]) + Rh[patch] + Vh[patch])))*r[patch,])
      
      dEh[((patch-1)*l)+1] <- sum(((alpha*pv) * ((Iv * Sh[patch]) / (Sh[patch] + sum(Eh[(l*(patch-1)+1):(l*patch)]) +
                              sum(I1h[(n1*(patch-1)+1):(n1*patch)]) + sum(I2h[(n2*(patch-1)+1):(n2*patch)]) + Rh[patch] + Vh[patch])))*r[patch,]) +
                              sum(((alpha*pv*(1-pf[patch])) * ((Iv*Vh[patch]) / (Sh[patch] + sum(Eh[(l*(patch-1)+1):(l*patch)]) +
                              sum(I1h[(n1*(patch-1)+1):(n1*patch)]) + sum(I2h[(n2*(patch-1)+1):(n2*patch)]) + Rh[patch] + Vh[patch])))*r[patch,]) - 
                              (l*epsilon*Eh[(l*(patch-1)+1)]) 
      
      for (i in 2:l){
        dEh[((patch-1)*l)+i] <- (l*epsilon*Eh[i-1+(l*(patch-1))]) - (l*epsilon*Eh[i+(l*(patch-1))])
      }
      
      
      
      dI1h[((patch-1)*n1)+1] <- (l*epsilon*Eh[(l*patch)]) - (n1*gamma1*I1h[(n1*(patch-1)+1)]) 
      
      for (i in 2:n1){
        dI1h[((patch-1)*n1)+i] <- (n1*gamma1*I1h[i-1+(n1*(patch-1))]) - (n1*gamma1*I1h[i+(n1*(patch-1))]) 
      }
      
      dI2h[((patch-1)*n2)+1] <- ((1-mh)*n1*gamma1*I1h[(n1*patch)]) - (n2*gamma2*I2h[(n2*(patch-1)+1)]) 
      
      for (i in 2:n2){
        dI2h[((patch-1)*n2)+i] <- (n2*gamma2*I2h[i-1+(n2*(patch-1))]) - (n2*gamma2*I2h[i+(n2*(patch-1))]) 
      }
      
      dRh[patch] <- n2 * gamma2 * I2h[n2*patch]
      
      dDh[patch] <- mh * n1 * gamma1 * I1h[n1*patch]
      
      dSv[patch] <- (muv * (Sv[patch]+sum(Ev[(k*(patch-1)+1):(k*patch)])+Iv[patch])) - 
                    sum(((alpha*ph*Sv[patch]) * ((colSums(matrix(I1h,nrow=n1,ncol=n)) + colSums(matrix(I2h,nrow=n2,ncol=n))) / (Sh[patch] + sum(Eh[(l*(patch-1)+1):(l*patch)]) +
                    sum(I1h[(n1*(patch-1)+1):(n1*patch)]) + sum(I2h[(n2*(patch-1)+1):(n2*patch)]) + Rh[patch] + Vh[patch])))*r[patch,]) - 
                    (muv*Sv[patch]) + (omega*(Sv[patch]/(Sv[patch]+sum(Ev[(k*(patch-1)+1):(k*patch)])+Iv[patch])))
      
      
      dEv[((patch-1)*k)+1] <- sum(((alpha * ph * Sv[patch]) * ((colSums(matrix(I1h,nrow=n1,ncol=n)) + colSums(matrix(I2h,nrow=n2,ncol=n))) / (Sh[patch] + sum(Eh[(l*(patch-1)+1):(l*patch)]) +
                              sum(I1h[(n1*(patch-1)+1):(n1*patch)]) + sum(I2h[(n2*(patch-1)+1):(n2*patch)]) + Rh[patch] + Vh[patch])))*r[patch,]) - 
                              ((k*upsilon) + muv) * Ev[(k*(patch-1)+1)] + (omega * (Ev[(k*(patch-1)+1)] / (Sv[patch] + sum(Ev[(k*(patch-1)+1):(k*patch)]) + Iv[patch])))
      for (i in 2:k){
        
        
      dEv[((patch-1)*k)+i] <- (k * upsilon * Ev[i-1+(k*(patch-1))]) - (((k*upsilon)+muv) * Ev[i+(k*(patch-1))]) + 
          (omega * (Ev[i+(k*(patch-1))] / (Sv[patch] + sum(Ev[(k*(patch-1)+1):(k*patch)]) + Iv[patch])))
      }
      
      
      dIv[patch] <- (k * upsilon * Ev[(k*patch)]) - (muv * Iv[patch]) + (omega * (Iv[patch] / (Sv[patch] + sum(Ev[(k*(patch-1)+1):(k*patch)]) + Iv[patch])))
      
      
      dVh[patch] <- -sum(((alpha*pv*(1-pf[patch])) * ((Iv* Vh[patch]) / (Sh[patch] + sum(Eh[(l*(patch-1)+1):(l*patch)]) +
                      sum(I1h[(n1*(patch-1)+1):(n1*patch)]) + sum(I2h[(n2*(patch-1)+1):(n2*patch)]) + Rh[patch] + Vh[patch])))*r[patch,])
      
    }
    
    # end of patch iteration loop
    
    # return the final derivatives containing results for all patches inside a list
    derivatives <- list(c(dSh,
                          dEh,
                          dI1h,
                          dI2h,
                          dRh,
                          dDh,
                          dSv,
                          dEv,
                          dIv,
                          dVh))
    return(derivatives)
  }
  
  #set up a data frame to store the outputs of all 12 outputs (each output starts one month later, which has been simplified to 30 days)
  outcombinedtotal <- data.frame()
  
  #iterate through the 12 months
  start_month_max <- 12
  
  for(start_month in 1:start_month_max){
    
    outcombinedmonth <- data.frame()
    
    # define the parameters and inputs for the function
    
    # loop through the neighbours to set up the parameters and initial conditions
    
    l <- 5 # number of stages in Eh
    n1 <- 11 # number of stages in I1h
    n2 <- 13 # number of stages in I2h
    k <- 10 # number of stages in Ev
    n <- length(neighbours) #number of patches in total including central patch and neighbours
    epsilon <- 1/4.6 #latent period for horses (days) = 1/epsilon
    gamma1 <- 0.26 #overall recovery rate from 1st infectious horse class
    gamma2 <- 1.25 #overall recovery rate from 2nd infectious horse class
    mh <- 0.86 #mortality horses
    ph <- 0.52 #transmission probability host to vector
    pv <- 0.77 #transmission probability vector to host
    
    
    pf <- c()
    Sh <- c()
    Sv <- c()
    Ev <- c()
    Iv <- c()
    Vh <- c()
    
    #iterate through the cells again. Used the term "area" since used "patch" extensively in the model
    
    for (area in 1:n){
      
      cellrow <- neighbours[area]
      
      #pre-define some parameters outside of the list to make it all work
      
      pf[area] <- AHS_dataset[[cellrow,5]] #protection factor of the vaccine
      
      Sh[area] <- (AHS_dataset[[cellrow,2]] - (AHS_dataset[[cellrow,2]] * AHS_dataset[[cellrow,4]]))  
      
      Eh <- c(1,rep.int(0,(l*n)-1))
      
      I1h <- rep.int(0,n1*n)
      
      I2h <- rep.int(0,n2*n)
      
      Rh <- rep.int(0,n)
      
      Dh <- rep.int(0,n)
      
      Sv[area] <- (midge_df[[cellrow]](1 + ((start_month * 30) - 30)))  
      
      Ev <- rep.int(0,k*n)
      
      Iv <- rep.int(0,n)
      
      Vh[area] <- ((AHS_dataset[[cellrow,2]] * AHS_dataset[[cellrow,4]])) 
      
    }
    
    # Define the matrix for the force of infection rho (based on Gaussian kernel as per Cmille Szmaragd and Simon Gubbin's work - https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007741)
    # Neighbour to neighbour has a rate based on distance. In a landlocked cell, the average distance from neighbour centroid to neighbour centroid
    # is 33.3km (i.e. (30km + 30km + 40km)/3)
    # Any neighbour cell centroid to central cell centroid has a rate based on 20km distance.
    # Any cell to itself has a rate of 1
    

    
    r <- matrix(data= (0.034/sqrt(pi))*exp((-0.034^2)*(33.3^2)), nrow=n, ncol=n, dimnames = list(neighbours,neighbours))
    for (i in 1:n){
      for (j in 1:n){
        if (i == j){
          r[i,j] <- 1
        }
        else if (i == 1 || j == 1){
          r[i,j] <- (0.034/sqrt(pi))*exp((-0.034^2)*(20^2))
        }
      }
    }
    
    
    # Define the initial conditions, parameters and times to compute the solution  
    pars <- list(n=n,
                 epsilon = epsilon,  
                 l = l, 
                 gamma1 = gamma1,
                 n1 = n1,
                 gamma2 = gamma2,
                 n2 = n2,
                 mh = mh,
                 k = k,
                 ph = ph,
                 pv = pv,
                 pf = pf, 
                 r = r,
                 starttime = 1+((start_month-1)*30)) 
    
    init <- c(Sh,Eh,I1h,I2h,Rh,Dh,Sv,Ev,Iv,Vh)
    
    # sort out names...
    
    init_numeric <- as.numeric(unlist(init))
    names(init_numeric) <- c(rep("Sh", n),rep(c("Eh1","Eh2","Eh3","Eh4","Eh5"), n), 
                             rep(c("I1h1","I1h2","I1h3","I1h4","I1h5","I1h6","I1h7","I1h8","I1h9","I1h10","I1h11"),n),
                             rep(c("I2h1","I2h2","I2h3","I2h4","I2h5","I2h6","I2h7","I2h8","I2h9","I2h10","I2h11","I2h12","I2h13"), n),
                             rep("Rh", n), rep("Dh",n), rep("Sv", n), rep(c("Ev1","Ev2","Ev3","Ev4","Ev5","Ev6","Ev7","Ev8","Ev9","Ev10"), n),
                             rep("Iv", n), rep("Vh", n))
    
    
    times <- seq(from = 1+((start_month-1)*30), to = 360+((start_month-1)*30), by = 1) 
    
    
    
    # run model using lsoda. 
    

    out <- lsoda(func = AHS.coupling.model, y = init_numeric, parms = pars, times = times, maxsteps = 50000,
                 events = list(func = eventfun, root = TRUE),
                 rootfun = rootfun)
    
    
    out <- as.data.frame(out) 
    
    #create additional columns to summarise output dataframe neatly
    
    for (i in 1:nrow(out)){
      
      out$sumSh[i] <- sum(out[i,2:(n+1)])
      out$sumEh[i] <- sum(out[i,(n+2) : ((6*n)+1)])
      out$sumIh[i] <- sum(out[i,((6*n)+2) : ((30*n)+1)])
      out$sumRh[i] <- sum(out[i,((30*n) +2) : ((31*n)+1)])
      out$sumVh[i] <- sum(out[i,((44*n) +2) : ((45*n)+1)])
      out$sumDh[i] <- sum(out[i,((31*n) +2) : ((32*n)+1)])
      out$sumSv[i] <- sum(out[i,((32*n) +2) : ((33*n)+1)])
      out$sumEv[i] <- sum(out[i,((33*n) +2): ((43*n)+1)])
      out$sumIv[i] <- sum(out[i,((43*n) +2) : ((44*n)+1)])
      
      out$cumulative[i] <- out$sumSh[i]+out$sumEh[i]+out$sumIh[i]+out$sumRh[i]+out$sumVh[i]+out$sumDh[i]
      out$totalinfv[i] <- out$sumEv[i]+out$sumIv[i]
      out$cumulativevec[i] <- out$sumSv[i]+out$sumEv[i]+out$sumIv[i]
      
    }
    
    #create data frame for each month start time output to bind onto the combined dataframe
    outcombinedmonth <- data.frame(Time = out$time, Susceptible = out$sumSh,
                                   Vaccinated = out$sumVh, Latent = out$sumEh, Infectious = out$sumIh,
                                   Recovered = out$sumRh, Dead = out$sumDh,
                                   Sv = out$sumSv, sumEv = out$sumEv,
                                   Iv = out$sumIv, sumvec = out$cumulativevec, Start_month = start_month)
    
    outcombinedmonth$day <- seq(1,360,1)
    
    outcombinedtotal <- rbind(outcombinedtotal, outcombinedmonth)
    
  }
  
  # Go through the combined output dataframe to subset each month start time.
  
  rows <- nrow(outcombinedtotal)
  subset1 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==1){
      subset1 <- rbind(subset1,outcombinedtotal[i,1:13])
    }
  }
  
  
  rows <- nrow(outcombinedtotal)
  subset2 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==2){
      subset2 <- rbind(subset2,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset3 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==3){
      subset3 <- rbind(subset3,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset4 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==4){
      subset4 <- rbind(subset4,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset5 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==5){
      subset5 <- rbind(subset5,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset6 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==6){
      subset6 <- rbind(subset6,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset7 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==7){
      subset7 <- rbind(subset7,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset8 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==8){
      subset8 <- rbind(subset8,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset9 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==9){
      subset9 <- rbind(subset9,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset10 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==10){
      subset10 <- rbind(subset10,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset11 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==11){
      subset11 <- rbind(subset11,outcombinedtotal[i,1:13])
    }
  }
  
  rows <- nrow(outcombinedtotal)
  subset12 <- data.frame()
  for (i in 1:rows){
    if (outcombinedtotal[i,12]==12){
      subset12 <- rbind(subset12,outcombinedtotal[i,1:13])
    }
  }
  
  # Summarise each subset to create a summary graph of the epidemiological outputs of interest - how long is the outbreak?
  # how many horses have been affected? How many horses were infected at the peak of the outbreak?
 
  summary <- data.frame()
  
  subset1$daily_change <- c(0,diff(subset1$Infectious))
  outbreaklength <- subset1$day[which((subset1$Infectious < 1)&(subset1$Latent < 1)&(subset1$daily_change < 0)&(subset1$day > 24))[1]]
  peakinfections <- subset1$day[which.max(subset1$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset1$Dead[which(subset1$day==outbreaklength)]
  recoveredhorses <- subset1$Recovered[which(subset1$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset1$Infectious[peakinfections]
  overwinter <- sum(subset1$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 1, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset2$daily_change <- c(0,diff(subset2$Infectious))
  outbreaklength <- subset2$day[which((subset2$Infectious < 1)&(subset2$Latent < 1)&(subset2$daily_change < 0)&(subset2$day > 24))[1]]
  peakinfections <- subset2$day[which.max(subset2$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset2$Dead[which(subset2$day==outbreaklength)]
  recoveredhorses <- subset2$Recovered[which(subset2$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset2$Infectious[peakinfections]
  overwinter <- sum(subset2$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 2, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset3$daily_change <- c(0,diff(subset3$Infectious))
  outbreaklength <- subset3$day[which((subset3$Infectious < 1)&(subset3$Latent < 1)&(subset3$daily_change < 0)&(subset3$day > 24))[1]]
  peakinfections <- subset3$day[which.max(subset3$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset3$Dead[which(subset3$day==outbreaklength)]
  recoveredhorses <- subset3$Recovered[which(subset3$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset3$Infectious[peakinfections]
  overwinter <- sum(subset3$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 3, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset4$daily_change <- c(0,diff(subset4$Infectious))
  outbreaklength <- subset4$day[which((subset4$Infectious < 1)&(subset4$Latent < 1)&(subset4$daily_change < 0)&(subset4$day > 24))[1]]
  peakinfections <- subset4$day[which.max(subset4$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset4$Dead[which(subset4$day==outbreaklength)]
  recoveredhorses <- subset4$Recovered[which(subset4$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset4$Infectious[peakinfections]
  overwinter <- sum(subset4$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 4, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset5$daily_change <- c(0,diff(subset5$Infectious))
  outbreaklength <- subset5$day[which((subset5$Infectious < 1)&(subset5$Latent < 1)&(subset5$daily_change < 0)&(subset5$day > 24))[1]]
  peakinfections <- subset5$day[which.max(subset5$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset5$Dead[which(subset5$day==outbreaklength)]
  recoveredhorses <- subset5$Recovered[which(subset5$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset5$Infectious[peakinfections]
  overwinter <- sum(subset5$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 5, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset6$daily_change <- c(0,diff(subset6$Infectious))
  outbreaklength <- subset6$day[which((subset6$Infectious < 1)&(subset6$Latent < 1)&(subset6$daily_change < 0)&(subset6$day > 24))[1]]
  peakinfections <- subset6$day[which.max(subset6$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset6$Dead[which(subset6$day==outbreaklength)]
  recoveredhorses <- subset6$Recovered[which(subset6$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset6$Infectious[peakinfections]
  overwinter <- sum(subset6$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 6, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset7$daily_change <- c(0,diff(subset7$Infectious))
  outbreaklength <- subset7$day[which((subset7$Infectious < 1)&(subset7$Latent < 1)&(subset7$daily_change < 0)&(subset7$day > 24))[1]]
  peakinfections <- subset7$day[which.max(subset7$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset7$Dead[which(subset7$day==outbreaklength)]
  recoveredhorses <- subset7$Recovered[which(subset7$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset7$Infectious[peakinfections]
  overwinter <- sum(subset7$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 7, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset8$daily_change <- c(0,diff(subset8$Infectious))
  outbreaklength <- subset8$day[which((subset8$Infectious < 1)&(subset8$Latent < 1)&(subset8$daily_change < 0)&(subset8$day > 24))[1]]
  peakinfections <- subset8$day[which.max(subset8$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset8$Dead[which(subset8$day==outbreaklength)]
  recoveredhorses <- subset8$Recovered[which(subset8$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset8$Infectious[peakinfections]
  overwinter <- sum(subset8$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 8, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset9$daily_change <- c(0,diff(subset9$Infectious))
  outbreaklength <- subset9$day[which((subset9$Infectious < 1)&(subset9$Latent < 1)&(subset9$daily_change < 0)&(subset9$day > 24))[1]]
  peakinfections <- subset9$day[which.max(subset9$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset9$Dead[which(subset9$day==outbreaklength)]
  recoveredhorses <- subset9$Recovered[which(subset9$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset9$Infectious[peakinfections]
  overwinter <- sum(subset9$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 9, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset10$daily_change <- c(0,diff(subset10$Infectious))
  outbreaklength <- subset10$day[which((subset10$Infectious < 1)&(subset10$Latent < 1)&(subset10$daily_change < 0)&(subset10$day > 24))[1]]
  peakinfections <- subset10$day[which.max(subset10$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset10$Dead[which(subset10$day==outbreaklength)]
  recoveredhorses <- subset10$Recovered[which(subset10$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset10$Infectious[peakinfections]
  overwinter <- sum(subset10$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 10, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset11$daily_change <- c(0,diff(subset11$Infectious))
  outbreaklength <- subset11$day[which((subset11$Infectious < 1)&(subset11$Latent < 1)&(subset11$daily_change < 0)&(subset11$day > 24))[1]]
  peakinfections <- subset11$day[which.max(subset11$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset11$Dead[which(subset11$day==outbreaklength)]
  recoveredhorses <- subset11$Recovered[which(subset11$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset11$Infectious[peakinfections]
  overwinter <- sum(subset11$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 11, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  subset12$daily_change <- c(0,diff(subset12$Infectious))
  outbreaklength <- subset12$day[which((subset12$Infectious < 1)&(subset12$Latent < 1)&(subset12$daily_change < 0)&(subset12$day > 24))[1]]
  peakinfections <- subset12$day[which.max(subset12$Infectious[1:outbreaklength])[1]]
  deadhorses <- subset12$Dead[which(subset12$day==outbreaklength)]
  recoveredhorses <- subset12$Recovered[which(subset12$day==outbreaklength)]
  totalaffected <- deadhorses + recoveredhorses
  numberpeak <- subset12$Infectious[peakinfections]
  overwinter <- sum(subset12$Infectious[outbreaklength:360] >1, na.rm = TRUE) > 0
  
  summary <- rbind(summary, data.frame(Month = 12, peakinfections, outbreaklength, deadhorses, recoveredhorses, totalaffected, numberpeak, overwinter))
  
  
  summary$cellrow <- rep(cell,12)
  
  summary_combined <- rbind(summary_combined, summary)
  
}

############## end of loop #############

#write_xlsx(summary_combined, "summarytotal_240522_121to153.xlsx")




# Explore all the data by developing graphs but must first adjust loop code so it is only run for one cell
# p graphs are horse outbreak graphs
# y graphs are infected/latent midge graphs
# z graphs look at the sum of all the vectors

# Explore outbreak starting in January

p1 <- ggplot() +
  geom_line(data = subset1, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset1, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset1, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset1, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset1, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset1, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y1 <- ggplot() +
  geom_line(data = subset1, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset1, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors")+
  ylim(-2,3100)


z1 <- ggplot() +
  geom_line(data = subset1, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 1))


# Explore outbreak starting in February

p2 <- ggplot() +
  geom_line(data = subset2, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset2, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset2, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset2, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset2, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset2, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")

y2 <- ggplot() +
  geom_line(data = subset2, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset2, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ylim(-2,3100)


z2 <- ggplot() +
  geom_line(data = subset2, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 2))


# Explore outbreak starting in March

p3 <- ggplot() +
  geom_line(data = subset3, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset3, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset3, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset3, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset3, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset3, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y3 <- ggplot() +
  geom_line(data = subset3, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset3, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors")+
  ylim(-2,3100)


z3 <- ggplot() +
  geom_line(data = subset3, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 3))


# Explore outbreak starting in April

p4 <- ggplot() +
  geom_line(data = subset4, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset4, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset4, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset4, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset4, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset4, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y4 <- ggplot() +
  geom_line(data = subset4, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset4, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ylim(-2,3100)


z4 <- ggplot() +
  geom_line(data = subset4, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 4))



# Explore outbreak starting in May

p5 <- ggplot() +
  geom_line(data = subset5, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset5, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset5, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset5, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset5, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset5, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y5 <- ggplot() +
  geom_line(data = subset5, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset5, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors")+
  ylim(-2,3100)


z5 <- ggplot() +
  geom_line(data = subset5, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 5))



# Explore outbreak starting in June

p6 <- ggplot() +
  geom_line(data = subset6, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset6, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset6, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset6, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset6, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset6, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y6 <- ggplot() +
  geom_line(data = subset6, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset6, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors")+
  ylim(-2,3100)


z6 <- ggplot() +
  geom_line(data = subset6, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 6))



# Explore outbreak starting in July

p7 <- ggplot() +
  geom_line(data = subset7, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset7, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset7, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset7, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset7, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset7, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y7 <- ggplot() +
  geom_line(data = subset7, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset7, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors")+
  ylim(-2,3100)


z7<- ggplot() +
  geom_line(data = subset7, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 7))


# Explore outbreak starting in August

p8 <- ggplot() +
  geom_line(data = subset8, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset8, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset8, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset8, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset8, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset8, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y8 <- ggplot() +
  geom_line(data = subset8, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset8, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ylim(-2,3100)


z8 <- ggplot() +
  geom_line(data = subset8, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 8))



# Explore outbreak starting in September

p9 <- ggplot() +
  geom_line(data = subset9, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset9, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset9, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset9, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset9, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset9, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y9 <- ggplot() +
  geom_line(data = subset9, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset9, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors")+
  ylim(-2,3100)


z9 <- ggplot() +
  geom_line(data = subset9, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 9))



# Explore outbreak starting in October

p10 <- ggplot() +
  geom_line(data = subset10, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset10, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset10, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset10, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset10, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset10, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y10 <- ggplot() +
  geom_line(data = subset10, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset10, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ylim(-2,3100)


z10  <- ggplot() +
  geom_line(data = subset10, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 10))




# Explore outbreak starting in November

p11 <- ggplot() +
  geom_line(data = subset11, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset11, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset11, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset11, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset11, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset11, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y11 <- ggplot() +
  geom_line(data = subset11, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset11, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ylim(-2,3100)


z11 <- ggplot() +
  geom_line(data = subset11, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 11))



# Explore outbreak starting in December

p12 <- ggplot() +
  geom_line(data = subset12, mapping = aes(x = day, y = Susceptible, color = "Susceptible"))+
  geom_line(data = subset12, mapping = aes(x = day, y = Latent, color = "Latent")) +
  geom_line(data = subset12, mapping = aes(x = day, y = Infectious, color = "Infectious")) +
  geom_line(data = subset12, mapping = aes(x = day, y = Recovered, color = "Recovered"))+
  geom_line(data = subset12, mapping = aes(x = day, y = Dead, color = "Dead")) +
  geom_line(data = subset12, mapping = aes(x = day, y = Vaccinated, color = "Vaccinated")) +
  scale_color_manual(values = c("Susceptible" = "green", "Latent" = "pink", "Infectious" = "red", "Recovered" = "purple", "Dead" = "black", "Vaccinated" = "orange"), limits=c("Susceptible", "Vaccinated", "Latent", "Infectious", "Recovered", "Dead")) +
  labs(color = "", x = "Time (days)", y = "Number of Horses")+
  ylim(-5,120)+
  ggtitle(" ")


y12 <- ggplot() +
  geom_line(data = subset12, mapping = aes(x = day, y = sumEv, color = "Latent")) +
  geom_line(data = subset12, mapping = aes(x = day, y = Iv, color = "Infectious"))+
  scale_color_manual(values = c("Latent" = "orange", "Infectious" = "red"), limits = c("Latent", "Infectious")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ylim(-2,3100)


z12 <- ggplot() +
  geom_line(data = subset12, mapping = aes(x = day, y = sumvec, color = "Sum Vectors")) +
  scale_color_manual(values = c("Sum Vectors" = "black")) +
  labs(color = "", x = "Time (days)", y = "Number of Vectors") +
  ggtitle(paste0("Outbreak starting in month ", 12))


# summary_combined <- read_excel(results_grid_gaus_fixed.xlsx)


s1 <- ggplot() +
  geom_line(data=summary_combined, mapping = aes(x = Month, y=totalaffected, color="Total Affected"))+
  scale_color_manual(values = c("Total Affected" = "red"))+
  labs(color = "", x = "Month of Start", y = "Number of Horses") +
  ggtitle("Number of horses affected")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))

s2 <- ggplot() +
  geom_line(data=summary_combined, mapping = aes(x = Month, y=peakinfections, color="Peak Infection Day"))+
  geom_line(data=summary_combined, mapping = aes(x = Month, y=outbreaklength, color="Outbreak Length"))+
  scale_color_manual(values = c("Peak Infection Day" = "red", "Outbreak Length" = "blue"))+
  labs(color = "", x = "Month of Start", y = "Days") +
  ggtitle("Outbreak length and time to peak infections")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))


summary_graph <- ggarrange(s1, s2,
                           labels = c(1,2),
                           ncol = 1, nrow = 2)



outbreak_graph <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
                            labels = c("A","B","C","D","E","F","G","H","I","J","K","L"),
                            ncol = 4, nrow = 3, common.legend = TRUE,
                            legend = "bottom")



midges_graph <- ggarrange(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12,
                          labels = c("A","B","C","D","E","F","G","H","I","J","K","L"),
                          ncol = 4, nrow = 3, common.legend = TRUE,
                          legend = "bottom")


