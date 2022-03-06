
library(shiny)
library(mvtnorm)
library(tibble)
library(DT)
library(MASS)

# Start
EM.Mixture.NH = function(data, k, p=NULL, mu=NULL, sigma=NULL, epsilon = 1e-8, maxit = 1000) {
  
  data <- as.matrix(data)
  nc = ncol(data) ; nr = nrow(data)
  
  f = matrix(0, nrow=nr, ncol=k)
  
  obs.loglik = function(f) {
    sum(log(rowSums(f)))
  }  
  
  zhat <- matrix(0, nrow=nr, ncol=k)
  sum.z <- numeric(k)
  
  
  # If initial values are unknown : Estimate initial values using Kmeans
  
  if (is.null(p) == TRUE | is.null(mu) == TRUE | is.null(sigma) == TRUE) {
    nc <- ncol(data)
    kmean <- kmeans(data, centers = k, nstart=3)
    
    if (is.null(p) == TRUE) { p <- kmean$size/nr }
    if (is.null(mu) == TRUE) { mu <- kmean$centers }
    if (is.null(sigma) == TRUE) {
      sigma <- array(0,dim=c(nc,nc,k))
      pred.dat <- cbind(data,kmean$cluster)
      init.sigma <- function(data, k) {
        nc = ncol(data)
        
        if (nc==2){
          for ( i in 1:k ) {
            sigma[,,i] <- sqrt(var(data[data[,nc]==i,-nc]))
          }          
        }
        else if (nc!=2)
          for ( i in 1:k ) {
            sigma[,,i] <- cov(data[data[,nc]==i,-nc])
          }
        sigma
      }
      sigma <- init.sigma(pred.dat,k)
    }
    
  }
  
  # 1st loop in order to record error & observed log-likelihood
  
  if (nc != 1) {
    
    for (j in 1:k) {
      
      f[,j] <- matrix(p[j]*dmvnorm(data, mu[j,], sigma[,,j]))
      
    }  
    
    denominator = rowSums(f)
    
    zhat = f/denominator
    
    sum.z <- colSums(zhat)
    
    zsum <- sum(sum.z)
    
    obs_ll <- obs.loglik(f)
    
    # Update Mixing Probability
    p <- sum.z/zsum
    
    # Update Mean Matrix
    
    for (j in 1:k) {
      
      for (l in 1:nc) {
        
        mu[j,l] <- sum(zhat[,j]*data[,l])/sum.z[j] 
        
      }
      
    }
    
    # Update Variance-Covariance Matrix
    
    for (j in 1:k) {
      
      for (l in 1:nc) {
        
        for (v in 1:l) {
          
          sigma[v,l,j] <- sum(zhat[,j] * (data[,l] - mu[j,l]) * (data[,v] - mu[j,v])) / sum.z[j]
          
        }
        
      }
      
      sigma[,,j][lower.tri(sigma[,,j])] <- t(sigma[,,j])[lower.tri(sigma[,,j])]
      
    }
    
    # Assign Error & Iteration
    
    error = err.set <- 100 
    iter <-1 # 1st loop already done
    
    # After calculating 1st loop
    
    ### E-step
    while (error > epsilon & iter < maxit) {
      
      iter <- iter + 1
      
      for (j in 1:k) {
        
        f[,j] <- matrix(p[j]*dmvnorm(data, mu[j,], sigma[,,j]))
        
      }  
      
      denominator = rowSums(f)
      
      zhat = f/denominator
      
      sum.z <- colSums(zhat)
      
      zsum <- sum(sum.z)
      
      ### M-step
      
      # Update Mixing Probability
      p <- sum.z/zsum
      
      # Update Mean Matrix
      
      for (j in 1:k) {
        
        for (l in 1:nc) {
          
          mu[j,l] <- sum(zhat[,j]*data[,l])/sum.z[j] 
          
        }
        
      }
      
      # Update Variance-Covariance Matrix
      
      for (j in 1:k) {
        
        for (l in 1:nc) {
          
          for (v in 1:l) {
            
            sigma[v,l,j] <- sum(zhat[,j] * (data[,l] - mu[j,l]) * (data[,v] - mu[j,v])) / sum.z[j]
            
          }
          
        }
        
        sigma[,,j][lower.tri(sigma[,,j])] <- t(sigma[,,j])[lower.tri(sigma[,,j])]
        
      }
      
      # Stack Observed Log-Likelihood
      obs_ll <- c(obs_ll, obs.loglik(f))
      
      # Record Error
      error <- obs_ll[length(obs_ll)] - obs_ll[length(obs_ll)-1]
      
      # Stack Error
      err.set <- c(err.set, error)
    }
  }
  
  ### Univariate Case ###
  
  else if (nc == 1) {
    
    for (j in 1:k) {
      
      f[,j] <- matrix(p[j]*dnorm(data, mu[j,], sigma[,,j]))
      
    }  
    
    denominator = rowSums(f)
    
    zhat = f/denominator
    
    sum.z <- colSums(zhat)
    
    zsum <- sum(sum.z)
    
    obs_ll <- obs.loglik(f)
    
    # Update Mixing Probability
    p <- sum.z/zsum
    
    # Update Variance-Covariance Matrix
    
    sigma2 <- array(0,dim=c(nc,nc,k))
    
    for (j in 1:k) {
      
      sigma2[,,j] <- sum(zhat[,j] * (data - mu[j,])^2) / sum.z[j]
      
    }
    
    sigma <- sqrt(sigma2)
    
    # Update Mean Matrix
    
    for (j in 1:k) {
      
      mu[j,] <- sum(zhat[,j]*data)/sum.z[j] 
      
    }
    
    # Assign Error & Iteration
    
    error = err.set <- 100 
    iter <-1 # 1st loop already done
    
    # After calculating 1st loop
    
    ### E-step
    while (error > epsilon & iter < maxit) {
      
      iter <- iter + 1
      
      for (j in 1:k) {
        
        f[,j] <- matrix(p[j]*dnorm(data, mu[j,], sigma[,,j]))
        
      }  
      
      denominator = rowSums(f)
      
      zhat = f/denominator
      
      sum.z <- colSums(zhat)
      
      zsum <- sum(sum.z)
      
      ### M-step
      
      # Update Mixing Probability
      p <- sum.z/zsum
      
      # Update Variance-Covariance Matrix
      
      for (j in 1:k) {
        
        sigma2[,,j] <- sum(zhat[,j] * (data - mu[j,])^2) / sum.z[j]
        
      }
      
      sigma <- sqrt(sigma2)
      
      # Update Mean Matrix
      
      for (j in 1:k) {
        
        mu[j,] <- sum(zhat[,j]*data)/sum.z[j] 
        
      }
      
      # Stack Observed Log-Likelihood
      obs_ll <- c(obs_ll, obs.loglik(f))
      
      # Record Error
      error <- obs_ll[length(obs_ll)] - obs_ll[length(obs_ll)-1]
      
      # Stack Error
      err.set <- c(err.set, error)
    }
  }
  
  if (error > epsilon) { word <- "Warning : Not converge!" }
  else if (error < epsilon) { word <- "EM-algorithm converges!" }
  
  return(list(Iteration = iter, EM.p = p, EM.mu = mu, EM.sigma = sigma,
              Observed_Log.Lik = obs_ll, Increase = err.set[-1], converge = word))
  
}

ui <- fluidPage(
  tags$h1("Shiny App: EM-algorithm for Mixture Model"),
  tags$h2("Namhwa Lee (nlee098@ucr.edu)"),
  tabsetPanel(
    tabPanel("Univariate", 
             fluidRow(
               column(3, wellPanel(
                 tags$h4("Data Generation & Estimation"),
                 numericInput(inputId = "size",
                              label = "Sample Size",
                              value = 1000),
                 textInput('p.vec',
                           label = 'Enter a true probability vector (sum to be 1, comma delimited)',
                           value="0.5,0.5"),
                 textInput('mu.vec',
                           label = 'Enter a true mean vector (comma delimited)',
                           value="0,2"),
                 textInput('sigma.vec',
                           label = 'Enter a true sigma vector (comma delimited, should be positive)',
                           value="0.5,1"),
               ),),
               column(6, plotOutput("hist")),
               column(3,
                      tags$h5("Estimation Result"),
                      tableOutput("Estimate")),
               column(3, 
                      tags$h5("Does EM converge?"),
                      textOutput("converge")),
               column(3, 
                      tags$h5("Iteration until EM stops"),
                      textOutput("iter"))
             ),
             fluidRow(
               column(3, wellPanel(
                 selectInput(inputId = "num",
                             label = "Number of Component (for Estimation)",
                             choices = list(2,3,4,5),
                             selected = 2),
                 textInput('init.p.vec',
                           label = 'Enter an initial probability vector (sum to be 1, comma delimited)'),
                 textInput('init.mu.vec',
                           label = 'Enter an initial mean vector (comma delimited)'),
                 textInput('init.sigma.vec',
                           label = 'Enter an initial sigma vector (comma delimited, should be positive)'),
               ),),
               column(6, plotOutput("obslik")),
               column(3,
                      tags$h5("Does EM work properly?"),
                      textOutput("increase")),
               column(3, wellPanel(
                 tags$h5("Stopping Rule"),
                 numericInput(inputId = "eps",
                              label = "Tolerance",
                              value = 1e-08),
                 numericInput(inputId = "maxit",
                              label = "Max Iteration",
                              value = 1000),
                 actionButton(inputId="go", label="Update"),
               ),)
             )
    ),
    tabPanel("Bivariate",
             fluidRow(
               column(3, wellPanel(
                 tags$h4("Data Generation & Estimation"),
                 numericInput(inputId = "size1",
                              label = "Sample Size",
                              value = 1000),
                 textInput('p.vec1',
                           label = 'Enter a probability vector (sum to be 1, comma delimited)',
                           value="0.5,0.5"),
                 textInput('mu.vec1',
                           label = 'Enter a mean vector (comma delimited)',
                           value="0,0,2,2"),
                 textInput('sigma.vec1',
                           label = 'Enter a covariance matrix (comma delimited, should be positive definite)',
                           value="1,0,0,1,1,0.5,0.5,1"),
                 checkboxGroupInput(inputId = "choice1", label = "Plot", choices = c("Sample", "True", "Estimated")),
               ),),
               column(6, plotOutput("contour")),
               column(1,
                      tags$h5("Mixing Probability"),
                      tableOutput("Estimate.p")),
               column(2,
                      tags$h5("Component Mean"),
                      tableOutput("Estimate.mu")),
               column(3,
                      tags$h5("Component Cov.mat"),
                      tableOutput("Estimate.sigma")),
               column(2,
                      tags$h5("Does EM converge?"),
                      textOutput("converge1")),
               column(1,
                      tags$h5("Iteration"),
                      textOutput("iter1"))
             ),
             fluidRow(
               column(3, wellPanel(
                 selectInput(inputId = "num1",
                             label = "Number of Component", choices = list(2,3,4,5),
                             selected = 2),
                 textInput('init.p.vec1',
                           label = 'Enter an initial probability vector (sum to be 1, comma delimited)'),
                 textInput('init.mu.vec1',
                           label = 'Enter an initial mean vector (comma delimited)'),
                 textInput('init.sigma.vec1',
                           label = 'Enter an initial sigma vector (comma delimited, should be positive)'),
               ),),
               column(6, plotOutput("obslik1")),
               column(3,
                      tags$h5("Does EM work properly?"),
                      textOutput("increase1")),
               column(3, wellPanel(
                 tags$h5("Stopping Rule"),
                 numericInput(inputId = "eps1",
                              label = "Tolerance",
                              value = 1e-08),
                 numericInput(inputId = "maxit1",
                              label = "Max Iteration",
                              value = 1000),
                 actionButton(inputId="go1", label="Update"),
               ),)
             )
    )
  )
)

gen.data <- function(n, m, p.vec, mu.vec, sigma.vec){
  # n <- input$size ; m <- input$num
  Z <- sample(1:m, size=n, prob=p.vec, replace=TRUE)
  rv <- function(samp, mu, sigma){
    for (j in 1:m) {
      if (samp == j) {
        return(rnorm(1, mu[j], sigma[j]))
      }
    }
  }
  return(as.matrix(sapply(Z, FUN=rv, mu=mu.vec, sigma=sigma.vec)))
}

gen.mvtdata <- function(n, m, p.vec, mu.mat, sigma.mat){
  # n <- input$size ; m <- input$num
  Z <- sample(1:m, size=n, prob=p.vec, replace=TRUE)
  rv <- function(samp, mu, sigma){
    for (j in 1:m) {
      if (samp == j) {
        return(rmvnorm(1, mu[,j], sigma[,,j]))
      }
    }
  }
  return(t(sapply(Z, FUN=rv, mu=mu.mat, sigma=sigma.mat)))
}

server <- function(input, output) {
  # Univariate Data
  data <- eventReactive(input$go, { 
    gen.data(n=input$size, m=length(as.numeric(unlist(strsplit(input$p.vec, ",")))),
             p.vec = as.numeric(unlist(strsplit(input$p.vec, ","))),
             mu.vec = as.numeric(unlist(strsplit(input$mu.vec, ","))),
             sigma.vec = as.numeric(unlist(strsplit(input$sigma.vec, ","))))
  })
  # Bivariate Data
  data2 <- eventReactive(input$go1, { 
    k <- length(as.numeric(unlist(strsplit(input$p.vec1, ","))))
    gen.mvtdata(n=input$size1, m=k,
                p.vec = as.numeric(unlist(strsplit(input$p.vec1, ","))),
                mu.mat = matrix(as.numeric(unlist(strsplit(input$mu.vec1, ","))),
                                ncol=k,
                                byrow = FALSE),
                sigma.mat = array(as.numeric(unlist(strsplit(input$sigma.vec1, ","))),
                                  dim=c(2,2,k)))
  })
  # Estimation: Univariate
  res <- eventReactive(input$go, {
    inputdata <- data()
    k.var <- ncol(inputdata)
    # All assigned X
    if (input$init.p.vec == "" & input$init.mu.vec == "" & input$init.sigma.vec == "") {
      EM.Mixture.NH(data=as.matrix(inputdata), k=as.numeric(input$num),
                    epsilon = input$eps, maxit = input$maxit)
    }
    # only init.p assigned
    else if (input$init.p.vec != "" & input$init.mu.vec == "" & input$init.sigma.vec == "") {
      EM.Mixture.NH(data=as.matrix(inputdata), k=as.numeric(input$num),
                    p = as.numeric(unlist(strsplit(input$init.p.vec, ","))),
                    epsilon = input$eps, maxit = input$maxit)
    }
    # only init.mu assigned
    else if (input$init.p.vec == "" & input$init.mu.vec != "" & input$init.sigma.vec == "") {
      EM.Mixture.NH(data=as.matrix(inputdata), k=as.numeric(input$num),
                    mu = matrix(as.numeric(unlist(strsplit(input$init.mu.vec, ","))),
                                nrow=length(unlist(strsplit(input$init.mu.vec, ",")))),
                    epsilon = input$eps, maxit = input$maxit)
    }
    # only init.sigma assigned : not work
    else if (input$init.p.vec == "" & input$init.mu.vec == "" & input$init.sigma.vec != "") {
      EM.Mixture.NH(data=as.matrix(inputdata), k=as.numeric(input$num),
                    sigma = array(as.numeric(unlist(strsplit(input$init.sigma.vec, ","))),
                                  dim=c(k.var,k.var,length(unlist(strsplit(input$init.sigma.vec, ",")))/(k.var^2))),
                    epsilon = input$eps, maxit = input$maxit)
    }
    # init.p & init.mu assigned
    else if (input$init.p.vec != "" & input$init.mu.vec != "" & input$init.sigma.vec == "") {
      EM.Mixture.NH(data=as.matrix(inputdata), k=as.numeric(input$num),
                    p = as.numeric(unlist(strsplit(input$init.p.vec, ","))),
                    mu = matrix(as.numeric(unlist(strsplit(input$init.mu.vec, ","))),
                                nrow=length(unlist(strsplit(input$init.mu.vec, ",")))),
                    epsilon = input$eps, maxit = input$maxit)
    }
    # init.p & init.sigma assigned: not work
    else if (input$init.p.vec != "" & input$init.mu.vec == "" & input$init.sigma.vec != "") {
      EM.Mixture.NH(data=as.matrix(inputdata), k=as.numeric(input$num),
                    p = as.numeric(unlist(strsplit(input$init.p.vec, ","))),
                    sigma = array(as.numeric(unlist(strsplit(input$init.sigma.vec, ","))),
                                  dim=c(k.var,k.var,length(unlist(strsplit(input$init.sigma.vec, ",")))/(k.var^2))),
                    epsilon = input$eps, maxit = input$maxit)
    }
    # init.mu & init.sigma assigned
    else if (input$init.p.vec == "" & input$init.mu.vec != "" & input$init.sigma.vec != "") {
      EM.Mixture.NH(data=as.matrix(inputdata), k=as.numeric(input$num),
                    mu = matrix(as.numeric(unlist(strsplit(input$init.mu.vec, ","))),
                                nrow=length(unlist(strsplit(input$init.mu.vec, ",")))),
                    sigma = array(as.numeric(unlist(strsplit(input$init.sigma.vec, ","))),
                                  dim=c(k.var,k.var,length(unlist(strsplit(input$init.sigma.vec, ",")))/(k.var^2))),
                    epsilon = input$eps, maxit = input$maxit)
    }
    # All assigned
    else if (input$init.p.vec != "" & input$init.mu.vec != "" & input$init.sigma.vec != "") {
      EM.Mixture.NH(data=as.matrix(inputdata), k=as.numeric(input$num),
                    p = as.numeric(unlist(strsplit(input$init.p.vec, ","))),
                    mu = matrix(as.numeric(unlist(strsplit(input$init.mu.vec, ","))),
                                nrow=length(unlist(strsplit(input$init.mu.vec, ",")))),
                    sigma = array(as.numeric(unlist(strsplit(input$init.sigma.vec, ","))),
                                  dim=c(k.var,k.var,length(unlist(strsplit(input$init.sigma.vec, ",")))/(k.var^2))),
                    epsilon = input$eps, maxit = input$maxit)
    }
  })
  # Estimation : Bivariate
  res2 <- eventReactive(input$go1, {
    inputdata1 <- data2()
    k.var <- ncol(inputdata1)
    # All assigned X
    if (input$init.p.vec1 == "" & input$init.mu.vec1 == "" & input$init.sigma.vec1 == "") {
      EM.Mixture.NH(data=inputdata1, k=as.numeric(input$num1),
                    epsilon = input$eps1, maxit = input$maxit1)
    }
    # only init.p assigned
    else if (input$init.p.vec1 != "" & input$init.mu.vec1 == "" & input$init.sigma.vec1 == "") {
      EM.Mixture.NH(data=as.matrix(inputdata1), k=as.numeric(input$num1),
                    p = as.numeric(unlist(strsplit(input$init.p.vec1, ","))),
                    epsilon = input$eps1, maxit = input$maxit1)
    }
    # only init.mu assigned
    else if (input$init.p.vec1 == "" & input$init.mu.vec1 != "" & input$init.sigma.vec1 == "") {
      EM.Mixture.NH(data=as.matrix(inputdata1), k=as.numeric(input$num1),
                    mu = matrix(as.numeric(unlist(strsplit(input$init.mu.vec1, ","))),
                                ncol=length(unlist(strsplit(input$init.mu.vec1, ",")))/as.numeric(input$num1),
                                byrow=FALSE),
                    epsilon = input$eps1, maxit = input$maxit1)
    }
    # only init.sigma assigned : not work
    else if (input$init.p.vec1 == "" & input$init.mu.vec1 == "" & input$init.sigma.vec1 != "") {
      EM.Mixture.NH(data=as.matrix(inputdata1), k=as.numeric(input$num1),
                    sigma = array(as.numeric(unlist(strsplit(input$init.sigma.vec1, ","))),
                                  dim=c(k.var,k.var,length(unlist(strsplit(input$init.sigma.vec1, ",")))/4)),
                    epsilon = input$eps1, maxit = input$maxit1)
    }
    # init.p & init.mu assigned
    else if (input$init.p.vec1 != "" & input$init.mu.vec1 != "" & input$init.sigma.vec1 == "") {
      EM.Mixture.NH(data=as.matrix(inputdata1), k=as.numeric(input$num1),
                    p = as.numeric(unlist(strsplit(input$init.p.vec1, ","))),
                    mu = matrix(as.numeric(unlist(strsplit(input$init.mu.vec1, ","))),
                                ncol=length(unlist(strsplit(input$init.mu.vec1, ",")))/as.numeric(input$num1),
                                byrow=FALSE),
                    epsilon = input$eps1, maxit = input$maxit1)
    }
    # init.p & init.sigma assigned:
    else if (input$init.p.vec1 != "" & input$init.mu.vec1 == "" & input$init.sigma.vec1 != "") {
      EM.Mixture.NH(data=as.matrix(inputdata1), k=as.numeric(input$num1),
                    p = as.numeric(unlist(strsplit(input$init.p.vec1, ","))),
                    sigma = array(as.numeric(unlist(strsplit(input$init.sigma.vec1, ","))),
                                  dim=c(k.var,k.var,length(unlist(strsplit(input$init.sigma.vec1, ",")))/4)),
                    epsilon = input$eps1, maxit = input$maxit1)
    }
    # init.mu & init.sigma assigned
    else if (input$init.p.vec1 == "" & input$init.mu.vec1 != "" & input$init.sigma.vec1 != "") {
      EM.Mixture.NH(data=as.matrix(inputdata1), k=as.numeric(input$num1),
                    mu = matrix(as.numeric(unlist(strsplit(input$init.mu.vec1, ","))),
                                ncol=length(unlist(strsplit(input$init.mu.vec1, ",")))/as.numeric(input$num1),
                                byrow=FALSE),
                    sigma = array(as.numeric(unlist(strsplit(input$init.sigma.vec1, ","))),
                                  dim=c(k.var,k.var,length(unlist(strsplit(input$init.sigma.vec1, ",")))/4)),
                    epsilon = input$eps1, maxit = input$maxit1)
    }
    # All assigned
    else if (input$init.p.vec1 != "" & input$init.mu.vec1 != "" & input$init.sigma.vec1 != "") {
      EM.Mixture.NH(data=as.matrix(inputdata1), k=as.numeric(input$num1),
                    p = as.numeric(unlist(strsplit(input$init.p.vec1, ","))),
                    mu = matrix(as.numeric(unlist(strsplit(input$init.mu.vec1, ","))),
                                ncol=length(unlist(strsplit(input$init.mu.vec1, ",")))/as.numeric(input$num1),
                                byrow=FALSE),
                    sigma = array(as.numeric(unlist(strsplit(input$init.sigma.vec1, ","))),
                                  dim=c(k.var,k.var,length(unlist(strsplit(input$init.sigma.vec1, ",")))/4)),
                    epsilon = input$eps1, maxit = input$maxit1)
    }
  })
  # Bivariate contour plot #
  output$contour <- renderPlot({
    bv.data <- data2()
    x <- seq(min(bv.data[,1])-1, max(bv.data[,1])+1, length.out=100)
    y <- seq(min(bv.data[,2])-1, max(bv.data[,2])+1, length.out=100)
    plot(NULL,
         xlim=c(min(bv.data[,1])-1, max(bv.data[,1])+1),
         ylim=c(min(bv.data[,2])-1, max(bv.data[,2])+1), ylab="y", xlab="x",
         main="Contour Plot")
    k.var <- ncol(bv.data)
    bivn.kde <- kde2d(bv.data[,1], bv.data[,2], n = 50)
    if ("Sample" %in% input$choice1 == TRUE) { 
      contour(bivn.kde, col="black", main="Contour Plot of Samples", lty=3, lwd=2, add=TRUE) 
    } 
    #contour(bivn.kde, col="black", main="Contour Plot of Samples", lty=3, lwd=2) 
    true.p <- isolate({as.numeric(unlist(strsplit(input$p.vec1, ",")))})
    true.mu <- isolate({matrix(as.numeric(unlist(strsplit(input$mu.vec1, ","))),
                               ncol=length(true.p),
                               byrow = FALSE)})
    true.sigma <- isolate({array(as.numeric(unlist(strsplit(input$sigma.vec1, ","))),
                                 dim=c(k.var,k.var,length(true.p)))})
    if (length(true.p) == 2) { # True: 2 component, 
      f <- function(X,Y, p=true.p, mu=true.mu, sigma=true.sigma) {
        return(p[1]*dmvnorm(cbind(X,Y), mean=mu[,1], sigma=sigma[,,1])
               +p[2]*dmvnorm(cbind(X,Y), mean=mu[,2], sigma=sigma[,,2]))
      }
      z <- outer(X=x, Y=y, f) 
      if ("True" %in% input$choice1 == TRUE) {
        contour(x=x, y=y, z=z, col="violetred", add=TRUE, lty=2, lwd=2)
      }
      if (length(res2()$EM.p) == 2) { # Est: 2 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])+p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 3) { # Est: 3 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 4) { # Est: 4 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3])
                 +p[4]*dmvnorm(cbind(X,Y), mean=mu[4,], sigma=sigma[,,4]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 5) { # Est: 5 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3])
                 +p[4]*dmvnorm(cbind(X,Y), mean=mu[4,], sigma=sigma[,,4])
                 +p[5]*dmvnorm(cbind(X,Y), mean=mu[5,], sigma=sigma[,,5]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
    }
    if (length(true.p) == 3) { # True: 3 component, 
      f <- function(X,Y, p=true.p, mu=true.mu, sigma=true.sigma) {
        return(p[1]*dmvnorm(cbind(X,Y), mean=mu[,1], sigma=sigma[,,1])
               +p[2]*dmvnorm(cbind(X,Y), mean=mu[,2], sigma=sigma[,,2])
               +p[3]*dmvnorm(cbind(X,Y), mean=mu[,3], sigma=sigma[,,3]))
      }
      z <- outer(X=x, Y=y, f) # when I want to draw true contour
      if ("True" %in% input$choice1 == TRUE) {
        contour(x=x, y=y, z=z, col="violetred", add=TRUE, lty=2, lwd=2)
      }
      if (length(res2()$EM.p) == 2) { # Est: 2 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])+p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 3) { # Est: 3 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 4) { # Est: 4 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3])
                 +p[4]*dmvnorm(cbind(X,Y), mean=mu[4,], sigma=sigma[,,4]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 5) { # Est: 5 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3])
                 +p[4]*dmvnorm(cbind(X,Y), mean=mu[4,], sigma=sigma[,,4])
                 +p[5]*dmvnorm(cbind(X,Y), mean=mu[5,], sigma=sigma[,,5]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
    }
    if (length(true.p) == 4) { # True: 4 component, 
      f <- function(X,Y, p=true.p, mu=true.mu, sigma=true.sigma) {
        return(p[1]*dmvnorm(cbind(X,Y), mean=mu[,1], sigma=sigma[,,1])
               +p[2]*dmvnorm(cbind(X,Y), mean=mu[,2], sigma=sigma[,,2])
               +p[3]*dmvnorm(cbind(X,Y), mean=mu[,3], sigma=sigma[,,3])
               +p[4]*dmvnorm(cbind(X,Y), mean=mu[,4], sigma=sigma[,,4]))
      }
      z <- outer(X=x, Y=y, f) # when I want to draw true contour
      if ("True" %in% input$choice1 == TRUE) {
        contour(x=x, y=y, z=z, col="violetred", add=TRUE, lty=2, lwd=2)
      }
      if (length(res2()$EM.p) == 2) { # Est: 2 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])+p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 3) { # Est: 3 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 4) { # Est: 4 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3])
                 +p[4]*dmvnorm(cbind(X,Y), mean=mu[4,], sigma=sigma[,,4]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 5) { # Est: 5 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3])
                 +p[4]*dmvnorm(cbind(X,Y), mean=mu[4,], sigma=sigma[,,4])
                 +p[5]*dmvnorm(cbind(X,Y), mean=mu[5,], sigma=sigma[,,5]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
    }
    if (length(true.p) == 5) { # True: 5 component, 
      f <- function(X,Y, p=true.p, mu=true.mu, sigma=true.sigma) {
        return(p[1]*dmvnorm(cbind(X,Y), mean=mu[,1], sigma=sigma[,,1])
               +p[2]*dmvnorm(cbind(X,Y), mean=mu[,2], sigma=sigma[,,2])
               +p[3]*dmvnorm(cbind(X,Y), mean=mu[,3], sigma=sigma[,,3])
               +p[4]*dmvnorm(cbind(X,Y), mean=mu[,4], sigma=sigma[,,4])
               +p[5]*dmvnorm(cbind(X,Y), mean=mu[,5], sigma=sigma[,,5]))
      }
      z <- outer(X=x, Y=y, f) # when I want to draw true contour
      if ("True" %in% input$choice1 == TRUE) {
        contour(x=x, y=y, z=z, col="violetred", add=TRUE, lty=2, lwd=2)
      }
      if (length(res2()$EM.p) == 2) { # Est: 2 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])+p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 3) { # Est: 3 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 4) { # Est: 4 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3])
                 +p[4]*dmvnorm(cbind(X,Y), mean=mu[4,], sigma=sigma[,,4]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
      else if (length(res2()$EM.p) == 5) { # Est: 5 Component
        est.f <- function(X,Y, p=res2()$EM.p, mu=res2()$EM.mu, sigma=res2()$EM.sigma) {
          return(p[1]*dmvnorm(cbind(X,Y), mean=mu[1,], sigma=sigma[,,1])
                 +p[2]*dmvnorm(cbind(X,Y), mean=mu[2,], sigma=sigma[,,2])
                 +p[3]*dmvnorm(cbind(X,Y), mean=mu[3,], sigma=sigma[,,3])
                 +p[4]*dmvnorm(cbind(X,Y), mean=mu[4,], sigma=sigma[,,4])
                 +p[5]*dmvnorm(cbind(X,Y), mean=mu[5,], sigma=sigma[,,5]))
        }
        z.est <- outer(X=x, Y=y, est.f)
        if ("Estimated" %in% input$choice1 == TRUE) {
          contour(x=x, y=y, z=z.est, col="steelblue", add=TRUE, lwd=2)
        }
      }
    }
  })
  output$hist <- renderPlot({
    dens <- density(data())$y
    hist(data(), freq=FALSE,  isolate({nclass=input$size/50}),
         main="Data and Estimation", xlab="data", ylim=c(0,max(dens)+0.1))
    true.p <- isolate({as.numeric(unlist(strsplit(input$p.vec, ",")))})
    true.mu <- isolate({as.numeric(unlist(strsplit(input$mu.vec, ",")))})
    true.sigma <- isolate({as.numeric(unlist(strsplit(input$sigma.vec, ",")))})
    if (length(true.p) == 2) { # True: 2 component, 
      curve(true.p[1]*dnorm(x, true.mu[1], true.sigma[1])
            +true.p[2]*dnorm(x, true.mu[2], true.sigma[2]),
            from=min(data()) - 1, to=max(data()) + 1,
            col="violetred", lty=2, n=1000, add=TRUE)
      if (length(res()$EM.p) == 2) { # Est: 2 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 3) { # Est: 3 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 4) { # Est: 4 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3])
              +res()$EM.p[4]*dnorm(x, res()$EM.mu[4,], res()$EM.sigma[,,4]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 5) { # Est: 5 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3])
              +res()$EM.p[4]*dnorm(x, res()$EM.mu[4,], res()$EM.sigma[,,4])
              +res()$EM.p[5]*dnorm(x, res()$EM.mu[5,], res()$EM.sigma[,,5]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
    }
    if (length(true.p) == 3) { # True: 3 component, 
      curve(true.p[1]*dnorm(x, true.mu[1], true.sigma[1])
            +true.p[2]*dnorm(x, true.mu[2], true.sigma[2])
            +true.p[3]*dnorm(x, true.mu[3], true.sigma[3]),
            from=min(data()) - 1, to=max(data()) + 1,
            col="violetred", lty=2, n=1000, add=TRUE)
      if (length(res()$EM.p) == 2) { # Est: 2 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 3) { # Est: 3 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 4) { # Est: 4 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3])
              +res()$EM.p[4]*dnorm(x, res()$EM.mu[4,], res()$EM.sigma[,,4]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 5) { # Est: 5 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3])
              +res()$EM.p[4]*dnorm(x, res()$EM.mu[4,], res()$EM.sigma[,,4])
              +res()$EM.p[5]*dnorm(x, res()$EM.mu[5,], res()$EM.sigma[,,5]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
    }
    if (length(true.p) == 4) { # True: 4 component, 
      curve(true.p[1]*dnorm(x, true.mu[1], true.sigma[1])
            +true.p[2]*dnorm(x, true.mu[2], true.sigma[2])
            +true.p[3]*dnorm(x, true.mu[3], true.sigma[3])
            +true.p[4]*dnorm(x, true.mu[4], true.sigma[4]),
            from=min(data()) - 1, to=max(data()) + 1,
            col="violetred", lty=2, n=1000, add=TRUE)
      if (length(res()$EM.p) == 2) { # Est: 2 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 3) { # Est: 3 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 4) { # Est: 4 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3])
              +res()$EM.p[4]*dnorm(x, res()$EM.mu[4,], res()$EM.sigma[,,4]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 5) { # Est: 5 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3])
              +res()$EM.p[4]*dnorm(x, res()$EM.mu[4,], res()$EM.sigma[,,4])
              +res()$EM.p[5]*dnorm(x, res()$EM.mu[5,], res()$EM.sigma[,,5]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
    }
    if (length(true.p) == 5) { # True: 5 component, 
      curve(true.p[1]*dnorm(x, true.mu[1], true.sigma[1])
            +true.p[2]*dnorm(x, true.mu[2], true.sigma[2])
            +true.p[3]*dnorm(x, true.mu[3], true.sigma[3])
            +true.p[4]*dnorm(x, true.mu[4], true.sigma[4])
            +true.p[5]*dnorm(x, true.mu[5], true.sigma[5]),
            from=min(data()) - 1, to=max(data()) + 1,
            col="violetred", lty=2, n=1000, add=TRUE)
      if (length(res()$EM.p) == 2) { # Est: 2 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 3) { # Est: 3 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 4) { # Est: 4 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3])
              +res()$EM.p[4]*dnorm(x, res()$EM.mu[4,], res()$EM.sigma[,,4]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
      else if (length(res()$EM.p) == 5) { # Est: 5 Component
        curve(res()$EM.p[1]*dnorm(x, res()$EM.mu[1,], res()$EM.sigma[,,1])
              +res()$EM.p[2]*dnorm(x, res()$EM.mu[2,], res()$EM.sigma[,,2])
              +res()$EM.p[3]*dnorm(x, res()$EM.mu[3,], res()$EM.sigma[,,3])
              +res()$EM.p[4]*dnorm(x, res()$EM.mu[4,], res()$EM.sigma[,,4])
              +res()$EM.p[5]*dnorm(x, res()$EM.mu[5,], res()$EM.sigma[,,5]),
              from=min(data()) - 1, to=max(data()) + 1,
              col="steelblue", n=1000, add=TRUE)
      }
    }
  })
  # Univariate
  output$Estimate <- renderTable({
    est <- tibble(
      p = round(c(res()$EM.p),4),
      mu = round(c(res()$EM.mu),4),
      sigma = round(c(res()$EM.sigma),4)
    )
    est$p <- formatC(est$p, digits = 4)
    est$mu <- formatC(est$mu, digits = 4)
    est$sigma <- formatC(est$sigma, digits =4)
    est
  })
  # Bivariate : How print different dimension output?
  output$Estimate.p <- renderTable({
    formatC(res2()$EM.p, digits=4)
  })
  output$Estimate.mu <- renderTable({
    formatC(res2()$EM.mu, digits=4)
  })
  output$Estimate.sigma <- renderTable({
    formatC(res2()$EM.sigma, digits=4)
  })
  # Tab1
  output$converge <- renderText({
    res()$converge
  })
  output$obslik <- renderPlot( {
    plot(1:res()$Iteration, res()$Observed_Log.Lik, type="l",
         main="Observed Log Likelihood", ylab="log-likelihood", xlab="Iteration")
  })
  output$iter <- renderText( {
    res()$Iteration
  })
  output$increase <- renderText({
    if (all(res()$Error > 0 ) == TRUE) {
      "Observed log-likelihood is non-decreasing."
    }
    else {
      "Warning: Observed lig-likelihood decreases. Your algorithm is wrong!"
    }
  })
  # Tab2
  output$converge1 <- renderText({
    res2()$converge
  })
  output$obslik1 <- renderPlot( {
    plot(1:res2()$Iteration, res2()$Observed_Log.Lik, type="l",
         main="Observed Log Likelihood", ylab="log-likelihood", xlab="Iteration")
  })
  output$iter1 <- renderText( {
    res2()$Iteration
  })
  output$increase1 <- renderText({
    if (all(res2()$Error > 0 ) == TRUE) {
      "Observed log-likelihood is non-decreasing."
    }
    else {
      "Warning: Observed lig-likelihood decreases. Your algorithm is wrong!"
    }
  })
}


shinyApp(ui = ui, server = server)