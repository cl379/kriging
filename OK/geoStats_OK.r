#################################################################################
# Title: Model-based ordinary kriging functions
# Authors: Jonas Alcaina-Mateos and Carla Lancelotti 
#           CaSEs Research Group - Department of HUmanities
#           Universitat Pompeu Fabra (UPF), Barcelona
#           j.alcaina.m@gmail.com / carla.lancelotti@upf.edu
# Comments: The present script include functions purposedly created by the author
#           for OK interpoletion methods, and functions called from the package
#           gstat (Pebesma 2004). In the script these last function are clearly
#           so that the user is aware of when they are called.
# Date: 25/03/2015
# url: https://github.com/cl379/kriging
#################################################################################
# ===============================================================================
cat('\n *** loading dependences and native functions ***\n\n')
require(sp)
require(rgdal)
require(raster)
require(gstat)
require(ggplot2)
require(gridExtra)
# ===============================================================================

#### Function that returns semivariogram models ####
.vgm <- function(model, nugg, sill, rang) {
    ## Linear model
    if (model == 'Lin') {
        return(function(x) {
            ifelse(x == 0, 
                0, 
                ifelse(x < rang, 
                    nugg + sill*(x/rang), 
                    nugg + sill))
        })
    } else
    ## Circular model
    if (model == 'Cir') {
        return(function(x) {
            ifelse(x == 0, 
                0, 
                ifelse(x < rang, 
                    nugg + sill*((2*x)/(pi*rang)*sqrt(1-(x/rang)^2)+(2/pi)*asin(x/rang)), 
                    nugg + sill))
        })
    } else
    # Spherical model
    if (model == 'Sph') {
        return(function(x) {
            ifelse(x == 0, 
                0, 
                ifelse(x < rang, 
                    nugg + sill*((3/2)*(x/rang) - (1/2)*(x/rang)^3), 
                    nugg + sill))
        })
    } else
    # Exponential model
    if (model == 'Exp') {
        return(function(x) {
            ifelse(x == 0, 
                0, 
                nugg + sill*(1 - exp(-x/rang)))
        })
    } else
    # Gaussian model
    if (model == 'Gau') {
        return(function(x) {
            ifelse(x == 0, 
                0, 
                nugg + sill*(1 - exp(-(x/rang)^2)))
        })
    } else
    # Potencial model
    if (model == 'Pow') {
        return(function(x) {
            ifelse(x == 0, 
                0, 
                nugg + sill*(x^rang))
        })
    } else 
    # Hole-effect model
    if (model == 'Hol') {
        return(function(x) {
            ifelse(x == 0, 
                0, 
                nugg + sill*(1 - rang * sin(x/rang) / x))
        })
    }
}

#### Ordiary kriging functions #### 

# Points estimation - returns a dataframe with a sinle row containing estimations and variance
.estimOk <- function (spdf, self.model, coords) {
    # A matrix
    A <- unname(as.matrix(dist(spdf@coords)))
    A <- self.model(A)
    A <- cbind(A, rep(1, nrow(A))) # Lagrange
    A <- rbind(A, c(rep(1, nrow(A)), 0)) # Lagrange
    # B vector
    B <- unname(as.matrix(dist(rbind(coords, spdf@coords))))[1,][-1]
    B <- self.model(B)
    B <- append(B, 1) # Lagrange
    # W vector
    W <- solve(A, B)
    # estimation
    return(data.frame(var1.pred = sum(spdf@data[,1] * W[1:nrow(spdf@data)]), 
                var1.var = sum(W[1:nrow(spdf@data)]*B[1:nrow(spdf@data)]) + W[nrow(spdf@data)+1]))
}

# Total estmation - returns a spatial pixel dataframe with the estimations and variance for each point in the grid
.Ok <- function (spdf, grd, model, nugg=0, sill=0, rang=0) {
    self.model <- .vgm(model, nugg, sill, rang)
    Ok <- Reduce(function(...) rbind(...), apply(grd@coords, 1, function(x) .estimOk(spdf, self.model, x)))
    return(SpatialPixelsDataFrame(grd, Ok))
}

# Leave-one-out validation - returns a data frame with the estimation and variance for each point in the spdf
.loo <- function(spdf, model, nugg=0, sill=0, rang=0) {
    self.model <- .vgm(model, nugg, sill, rang)
    estim <- data.frame()
    for (i in 1:nrow(spdf@data)) {
        estim <- rbind(estim, .estimOk(spdf[-i,], self.model, spdf[i,]@coords))
    }
    return(estim)
}

#### Fit and statistics ####

# function to calculate the fit statistics 
.fit.stats <- function(semvar, spdf, model, nugg=0, sill=0, rang=0) {
    # build objects
    self.model <- .vgm(model, nugg, sill, rang)
    loo <- .loo(spdf, model, nugg, sill, rang)
    # stats - returns a dtaframe with the value of different statististics of the model fitting
    value.1 <- sum((semvar$gamma - self.model(semvar$dist))^2) # least squares
    value.2 <- sum((semvar$gamma - self.model(semvar$dist))^2 / semvar$np) # least squares weighted
    value.3 <- sum((semvar$gamma - self.model(semvar$dist))^2 / (self.model(semvar$dist)^2/semvar$np)) # Cressie (1985)
    value.4 <- nrow(semvar) * log(sum((semvar$gamma - self.model(semvar$dist))^2)) + 6 # Akaike
    value.5 <- sum((spdf@data[,1] - loo[,1])^2) / nrow(spdf@data) # leave-one-out | (v-p)^2 / n
    value.6 <- sum((spdf@data[,1] - loo[,1])^2 / loo[,2]) / nrow(spdf@data) # leave-one-out variance 
    return(data.frame(ss=value.1, ss.w=value.2, Cressie=value.3, aic=value.4, loo=value.5, loo.var=value.6))
}

# function to minimize (least squares). To pass in the .fit function
.op <- function(par, data, model) {
    self.model <- .vgm(model, par[1], par[2], par[3])
    return(with(data, sum((gamma - self.model(dist))^2)))
}

# function to calculate the model parameters - return a dataframe with the parameters for the model specified
.fit <- function(semvar, model, spdf) {
    op <- optim(par=c(0, var(spdf@data[,1]), max(semvar$dist)/3), 
            .op, data = semvar, model = model, method = 'L-BFGS-B', lower=0, upper = max(semvar$dist, semvar$gamma))
    params <- data.frame(model = model, nugg = op$par[1], sill = op$par[2], rang = op$par[3])
    return(params)
}

# Model selection function - returns a datarame with the Cir, Sph, Exp and Gau models (parameters and fit stats)
.selection.model <- function(semvar, spdf) {
    stats <- data.frame()
    for (i in c('Cir', 'Sph', 'Exp', 'Gau')) {
        fit <- .fit(semvar, i, spdf)
        fit.stats <- .fit.stats(semvar, spdf, i, fit$nugg, fit$sill, fit$rang)
        stats <- rbind(stats, cbind(fit, fit.stats))
    }
    return(stats)
}

#### Model plots ####

# Plot a model. The function runs with a semivariogram dataframe, or ggplot 
.plot.model <- function(Object, model, nugg, sill, rang, Color='blue') {
    self.model <- .vgm(model, nugg, sill, rang)
    # Object: semivariogram (dataframe) / ggplot
    if (class(Object)[1] == 'gg') {
        self.plot <- Object
        self.plot <- self.plot + stat_function(fun=self.model, colour = Color)
    } else {
        semvar <- Object
        self.plot <- ggplot(data=semvar) + geom_point(aes(x=dist, y=gamma, alpha=np))
        self.plot <- self.plot + stat_function(fun=self.model, colour = Color)
    }
    return(self.plot)
}

# Plot multiple models - (only Cir, Sph, Exp, Gau) - requieres a .selection.model result  
.plot.sel <- function(stats, semvar) {
    Colors <- c('blue', 'red', 'orange3', 'green4')
    self.plot <- .plot.model(semvar, as.character(stats$model[1]), stats$nugg[1], stats$sill[1], stats$rang[1], Colors[1])
    self.plot <- .plot.model(self.plot, as.character(stats$model[2]), stats$nugg[2], stats$sill[2], stats$rang[2], Colors[2])
    self.plot <- .plot.model(self.plot, as.character(stats$model[3]), stats$nugg[3], stats$sill[3], stats$rang[3], Colors[3])
    self.plot <- .plot.model(self.plot, as.character(stats$model[4]), stats$nugg[4], stats$sill[4], stats$rang[4], Colors[4])
    return(self.plot)
}

#### Semivariogram function ####

# semivariogram function. Needs an object of SpatialPointsDataFrame
.semvar <- function(spdf, dimension=1, treshold=3, Method='All', nh=20) {
    ### pairing
    # create the distance vector
    distMat <- unname(as.matrix(dist(spdf@coords))) # dist function
    distMat[upper.tri(distMat, diag=T)] <- NA
    distVec <- na.omit(as.numeric(distMat))
    # create de square diferences
    vals <- as.matrix(spdf@data[,dimension]) # extract values
    difMat <- outer(1:nrow(vals),1:nrow(vals), FUN = Vectorize(function(i, j) (vals[i,] - vals[j,])^2))
    difMat[upper.tri(difMat, diag=T)] <- NA
    difVec <- na.omit(as.numeric(difMat))
    # Build data.frame and order by distance
    pairs <- data.frame(pitag = distVec, sqdif = difVec)
    pairs <- pairs[order(pairs$pitag),]
    # Calculate the relevant distance and subsample the pairs
    relevDist <- sqrt((max(spdf@coords[,1]) - min(spdf@coords[,1]))^2 + (max(spdf@coords[,2]) - min(spdf@coords[,2]))^2) / treshold
    pairs <- pairs[which(pairs$pitag < relevDist),]
    # Split classes by distances
    if (Method == 'All') {distGroup <- split(pairs, pairs$pitag)} else
    if (Method == 'Smooth') {
        distGroup <- list()
        interval <- (max(pairs$pitag) - min(pairs$pitag)) / nh
        for (i in 1:nh) {
            distGroup <- append(distGroup, list(pairs[pairs$pitag < min(pairs$pitag) + interval*i & pairs$pitag >= min(pairs$pitag) + interval*(i-1),]))
        }
    }
    # Cumulatives loops through breakes
    np <- c()
    dists <- c()
    gamma <- c()
    # Loops through breaks
    for (i in distGroup) {
        # Calculate the semivariance and establish number of pairs and distances tags 
	    np <- append(np, nrow(i))
	    dists <- append(dists, mean(i$pitag))
	    gamma <- append(gamma, (1 / (2 * nrow(i))) * sum(i$sqdif))
    }
    # Return a list, which the first object is the semivariogram and the second the variogram cloud
    return(list(data.frame(np=np, dist=dists, gamma=gamma), data.frame(dist=pairs$pitag, gamma=pairs$sqdif)))
}

# ggplot functions for exploring the variogram (needs a semivariogram and a variogram cloud)
.semvar.expl <- function(semvar, cloud) {
    p1 <- ggplot(data=cloud) + geom_boxplot(aes(x=dist, y=gamma, group=(dist)), fill = 'grey80', outlier.shape=NA, width=0.1)
    p2 <- ggplot(data=cloud) + geom_point(aes(x=dist, y=gamma), alpha = 1/10) 
    p3 <- ggplot(data=semvar) + geom_point(aes(x=dist, y=gamma))
    p4 <- ggplot(data=semvar) + geom_point(aes(x=dist, y=gamma, alpha=np)) + geom_smooth(aes(x=dist, y=gamma))
    grid.arrange(p1, p2, p3, p4, ncol = 2)
}

# ===============================================================================
cat('\n *** Geo Stats loaded ***\n\n')
