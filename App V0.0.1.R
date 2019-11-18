rm(list=ls())

#install.packages('DT', repos='http://cran.rstudio.com/')
#install.packages('MCMCpack', repos='http://cran.rstudio.com/')
#install.packages('rhandsontable', repos='http://cran.rstudio.com/')
#install.packages('shiny', repos='http://cran.rstudio.com/')
#install.packages('shinydashboard', repos='http://cran.rstudio.com/')
library(shiny) 
library(DT)
library(MCMCpack)
library(rhandsontable)


library(shinydashboard)


#### Input tables ####
#### 1. Prior ####

PriorE <- PriorC <- data.frame(success=c(0.5,0.5,1), failure=c(0.5,0.5,1), total=c(1,1,2))
colnames(PriorE) <- colnames(PriorC) <- c("success 2", "failure 2", "total 1")
rownames(PriorE) <- rownames(PriorE) <- c("success 1", "failure 1", "total 2")
#### 2. Data ####
DataE <- DataC <- data.frame(success=c(1,1,2), failure=c(1,1,2), total=c(2,2,4))
colnames(DataE) <- colnames(DataC) <- c("success 2", "failure 2", "total 1")
rownames(DataE) <- rownames(DataC) <- c("success 1", "failure 1", "total 2")

#### Global functions ####
#### 1. Sample Dirichlet #### 

SampleDir <- function (n, alpha){
  k = length(alpha)
  z = array(0, dim = c(n, k))
  for (i in 1:k) {z[, i] = rgamma(n,alpha[i])}#Callrgamma(n, as.double(alpha[i]),1/1.0)}
  s <- rowSums(z)
  z <- apply(z,2,"/",s)
  return(z)}

#### 2. Compute Posterior probability for 2 outcomes ####
PoP <- function(x11A, x10A, x01A, x00A, 
                x11B, x10B, x01B, x00B, 
                p11A=0.5, p10A=0.5, p01A=0.5, p00A=0.5,
                p11B=0.5, p10B=0.5, p01B=0.5, p00B=0.5){
  
  # Compute parameters Posterior distribution
  PostAlphaA <- c(x11A,x10A,x01A,x00A)+c(p11A,p10A,p01A,p00A)
  PostAlphaB <- c(x11B,x10B,x01B,x00B)+c(p11B,p10B,p01B,p00B)
  
  # Draw sample from Posterior distribution
  PostA <- SampleDir(5e4,PostAlphaA)
  PostB <- SampleDir(5e4,PostAlphaB)
  
  # Compute Posterior probability per decision rule
  theta1 <- mean(rowSums(PostA[,c(1,2)])-rowSums(PostB[,c(1,2)])>0)
  theta2 <- mean(rowSums(PostA[,c(1,3)])-rowSums(PostB[,c(1,3)])>0)
  all <- mean(rowSums(PostA[,c(1,2)])-rowSums(PostB[,c(1,2)])>0 &
                rowSums(PostA[,c(1,3)])-rowSums(PostB[,c(1,3)])>0)
  any <- mean(rowSums(PostA[,c(1,2)])-rowSums(PostB[,c(1,2)])>0 |
                rowSums(PostA[,c(1,3)])-rowSums(PostB[,c(1,3)])>0)
  comp <- mean((rowSums(PostA[,c(1,2)])+rowSums(PostA[,c(1,3)]))-
                 (rowSums(PostB[,c(1,2)])+rowSums(PostB[,c(1,3)]))>0)
  
  return(c(theta1,theta2,all,any,comp))
}

#### 3. Plot treatment distributions ####
PlotPost <- function(x11A, x10A, x01A, x00A, 
                     x11B, x10B, x01B, x00B, 
                     p11A=0.5, p10A=0.5, p01A=0.5, p00A=0.5,
                     p11B=0.5, p10B=0.5, p01B=0.5, p00B=0.5,
                     diff=FALSE){
  
  # Compute parameters Posterior distribution
  PostAlphaA <- c(x11A,x10A,x01A,x00A)+c(p11A,p10A,p01A,p00A)
  PostAlphaB <- c(x11B,x10B,x01B,x00B)+c(p11B,p10B,p01B,p00B)
  
  # Draw sample from Posterior distribution
  PostA <- SampleDir(5e4,PostAlphaA)
  PostB <- SampleDir(5e4,PostAlphaB)
  
  # Plot treatment distributions
  if(diff==FALSE){
    layout(matrix(c(1,1,2,3,4,4), nrow=3,ncol=2, byrow=TRUE), heights=c(32,16,4))
    
    # Two distributions in one plot 
    par(mar=c(5,5,2,0), mgp=c(3.7,1.5,0))
      contour(table(cut(rowSums(PostA[,c(1,2)]),breaks=seq(0,1,0.02)),
                  cut(rowSums(PostA[,c(1,3)]),breaks=seq(0,1,0.02))), 
            col="blue", frame.plot=FALSE,
            main="Experimental and control treatments",
            xlab="Outcome 1", 
            ylab="Outcome 2",
            xlim=c(0,1), ylim=c(0,1), las=1, lwd=2.50, nlevels=6,
            cex.main=1.5, cex.lab=1.5, cex.axis=1.5, drawlabels=FALSE)
    contour(table(cut(rowSums(PostB[,c(1,2)]),breaks=seq(0,1,0.02)),
                  cut(rowSums(PostB[,c(1,3)]),breaks=seq(0,1,0.02))),
            col="red", lwd=2.50, nlevels=6, drawlabels=FALSE, add=TRUE)
    
    # Experimental distribution in separate plot
    contour(table(cut(rowSums(PostA[,c(1,2)]),breaks=seq(0,1,0.02)),
                  cut(rowSums(PostA[,c(1,3)]),breaks=seq(0,1,0.02))), 
            col="blue", frame.plot=FALSE,
            main="Experimental treatment",
            xlab="Outcome 1", 
            ylab="Outcome 2",
            xlim=c(0,1), ylim=c(0,1), las=1, lwd=1.75, nlevels=6,
            cex.main=1.5, cex.lab=1.5, cex.axis=1.5, drawlabels=FALSE)
    
    # Control distribution in separate plot
    contour(table(cut(rowSums(PostB[,c(1,2)]),breaks=seq(0,1,0.02)),
                  cut(rowSums(PostB[,c(1,3)]),breaks=seq(0,1,0.02))), 
            col="red", frame.plot=FALSE,
            main="Control treatment",
            xlab="Outcome 1", 
            ylab="Outcome 2",
            xlim=c(0,1), ylim=c(0,1), las=1, lwd=1.75, nlevels=6,
            cex.main=1.5, cex.lab=1.5, cex.axis=1.5, drawlabels=FALSE)
    par(mar=c(0.1,0.1,0.1,0.1))
    plot.new()
    legend("center", legend=c("Experimental","Control"),
           col=c("blue", "red"), lty=1, ncol=2, cex=1.5)}
  
  # Distribution of the treatment difference
  if(diff==TRUE){
    DiffAB <- PostA-PostB
    contour(x=seq(-0.495,0.5,0.05), 
            y=seq(-0.495,0.5,0.05),
            z=table(cut(rowSums(DiffAB[,c(1,2)]),breaks=seq(-0.5,0.5,0.05)),
                    cut(rowSums(DiffAB[,c(1,3)]),breaks=seq(-0.5,0.5,0.05))), 
            col="#911eb4", frame.plot=FALSE, lwd=1.75, nlevels=6, drawlabels=FALSE,
            xlab="Difference outcome 1 (experimental-control)", 
            ylab="Difference outcome 2 (experimental-control)",
            xlim=c(-1,1), ylim=c(-1,1), las=1)
    abline(h=0, lty=2)
    abline(v=0, lty=2)
  }
}


#### 4. Compute correlations ####
#### 4.1 Correlation data ####
CorrData <- function(x11, x10, x01, x00){
  X <- c(x11, x10, x01, x00)
  n <- sum(X)
  Cov <- (X[1]*X[4]-X[2]*X[3])/
    sum(X)
  Var1 <- (X[1]+X[2])*(X[3]+X[4])/sum(X)
  Var2 <- (X[1]+X[3])*(X[2]+X[4])/sum(X)
  return(Cov/(sqrt(Var1*Var2)))
}
#### 4.2 Correlation prior ####
CorrPrior <- function(p11, p10, p01, p00){
  Alpha <- c(p11, p10, p01, p00)
  Cov <- (Alpha[1]*Alpha[4]-Alpha[2]*Alpha[3])
  Var1 <- ((Alpha[1]+Alpha[2])*(Alpha[3]+Alpha[4]))
  Var2 <- ((Alpha[1]+Alpha[3])*(Alpha[2]+Alpha[4]))
  return(Cov/(sqrt(Var1*Var2)))  
}


#### 4.3 Correlation posterior ####
CorrPost <- function(x11, x10, x01, x00,
                     p11=0.5, p10=0.5, p01=0.5, p00=0.5){
  PostAlpha <- c(x11, x10, x01, x00)+c(p11, p10, p01, p00)
  Cov <- (PostAlpha[1]*PostAlpha[4]-PostAlpha[2]*PostAlpha[3])
  Var1 <- ((PostAlpha[1]+PostAlpha[2])*(PostAlpha[3]+PostAlpha[4]))
  Var2 <- ((PostAlpha[1]+PostAlpha[3])*(PostAlpha[2]+PostAlpha[4]))
  return(Cov/(sqrt(Var1*Var2)))  
}



#### Server ####
server <- shinyServer(function(input, output, session) {

#### 1. Tab data ####
  # Table data experimental treatment
  previousXE <- reactive({DataE})
  
  changeXE <- reactive({
    if(is.null(input$hot_xe)){return(previousXE())}
    else{#if(!identical(previousXE(),input$hot_xe)){
       dataE <- as.data.frame(hot_to_r(input$hot_xe))
      # Perform marginal computations
      dataE[c(1,2),3] <- rowSums(dataE[c(1,2),c(1,2)])
      dataE[3,c(1,2)] <- colSums(dataE[c(1,2),c(1,2)])
      dataE[3,3] <- sum(dataE[c(1,2),3])
      dataE
    }
  })
  
  output$hot_xe <- renderRHandsontable({rhandsontable(changeXE(), width=375, 
                                                      rowHeaders=c("Success outcome 2", "Failure outcome 2", "Total outcome 1"),
                                                      colHeaders=c("Success outcome 1", "Failure outcome 1", "Total outcome 2")) %>%
      hot_table(rowHeaderWidth=150) %>%
      hot_cols(colWidths = c(75,75,75), format="0", allowInvalid=FALSE)%>%
      hot_col(3, readOnly=TRUE)%>%
      hot_row(3, readOnly=TRUE)})
 
    # Table control treatment
  previousXC <- reactive({DataC})
  
  changeXC <- reactive({
    if(is.null(input$hot_xc)){return(previousXC())}
    else{#if(!identical(previousXC(),input$hot_xc)){
      dataC <- as.data.frame(hot_to_r(input$hot_xc))
      # Perform marginal computations
      dataC[c(1,2),3] <- rowSums(dataC[c(1,2),c(1,2)])
      dataC[3,c(1,2)] <- colSums(dataC[c(1,2),c(1,2)])
      dataC[3,3] <- sum(dataC[c(1,2),3])
      dataC
    }
  })
  
  output$hot_xc <- renderRHandsontable({rhandsontable(changeXC(), width=375, 
                                                      rowHeaders=c("Success outcome 2", "Failure outcome 2", "Total outcome 1"),
                                                      colHeaders=c("Success outcome 1", "Failure outcome 1", "Total outcome 2")) %>%
      hot_table(rowHeaderWidth=150) %>%
      hot_cols(colWidths = c(75,75,75), format="0", allowInvalid=FALSE)%>%
      hot_col(3, readOnly=TRUE)%>%
      hot_row(3, readOnly=TRUE)})
  
  # Observed correlation
  output$Corr_data <- renderTable({
    
    if(is.null(input$hot_xe)){dataE.df <- DataE}
    else{dataE.df <- as.data.frame(hot_to_r(input$hot_xe))}
    if(is.null(input$hot_xc)){dataC.df <- DataC}
    else{dataC.df <- as.data.frame(hot_to_r(input$hot_xc))}
    
    data_corr_e <- CorrData(dataE.df[1,1], dataE.df[2,1], 
                              dataE.df[1,2], dataE.df[2,2])
    data_corr_c <- CorrData(dataC.df[1,1], dataC.df[2,1], 
                              dataC.df[1,2], dataC.df[2,2])
    data.frame("Treatment"=c("Experimental", "Control"), 
               "Correlation"=c(data_corr_e, data_corr_c))
  })
  
#### 2. Tab prior ####
  # Table prior experimental treatment
  previousPE <- reactive({PriorE})
  
  changePE <- reactive({
    if(is.null(input$hot_pe)){return(previousPE())}
    else #if(!identical(previousPE(),input$hot_pe))
      {
      priorE <- as.data.frame(hot_to_r(input$hot_pe))
      # Perform marginal computations
      priorE[c(1,2),3] <- rowSums(priorE[c(1,2),c(1,2)])
      priorE[3,c(1,2)] <- colSums(priorE[c(1,2),c(1,2)])
      priorE[3,3] <- sum(priorE[c(1,2),3])
      priorE
    }
  })
  
  output$hot_pe <- renderRHandsontable({rhandsontable(changePE(), width=375, 
                                                      rowHeaders=c("Success outcome 2", "Failure outcome 2", "Total outcome 1"),
                                                      colHeaders=c("Success outcome 1", "Failure outcome 1", "Total outcome 2")) %>%
      hot_table(rowHeaderWidth=150) %>%
      hot_cols(colWidths = c(75,75,75), format="0", allowInvalid=FALSE)%>%
      hot_col(3, readOnly=TRUE)%>%
      hot_row(3, readOnly=TRUE)})
  
  # Table control treatment
  previousPC <- reactive({PriorC})
  
  changePC <- reactive({
    if(is.null(input$hot_pc)){return(previousPC())}
    else{ #if(!identical(previousPC(),input$hot_pc)){
      priorC <- as.data.frame(hot_to_r(input$hot_pc))
      # Perform marginal computations
      priorC[c(1,2),3] <- rowSums(priorC[c(1,2),c(1,2)])
      priorC[3,c(1,2)] <- colSums(priorC[c(1,2),c(1,2)])
      priorC[3,3] <- sum(priorC[c(1,2),3])
      priorC
    }
  })
  
  output$hot_pc <- renderRHandsontable({rhandsontable(changePC(), width=375, 
                                                      rowHeaders=c("Success outcome 2", "Failure outcome 2", "Total outcome 1"),
                                                      colHeaders=c("Success outcome 1", "Failure outcome 1", "Total outcome 2")) %>%
      hot_table(rowHeaderWidth=150) %>%
      hot_cols(colWidths = c(75,75,75), format="0", allowInvalid=FALSE)%>%
      hot_col(3, readOnly=TRUE)%>%
      hot_row(3, readOnly=TRUE)})
  
  # Prior correlation
  output$Corr_prior <- renderTable({
    if(is.null(input$hot_pe)){priorE.df <- PriorE}
    else{priorE.df <- as.data.frame(hot_to_r(input$hot_pe))}
    if(is.null(input$hot_pc)){priorC.df <- PriorC}
    else{priorC.df <- as.data.frame(hot_to_r(input$hot_pc))}
    
    prior_corr_e <- CorrPrior(priorE.df[1,1], priorE.df[2,1], 
                              priorE.df[1,2], priorE.df[2,2])
    prior_corr_c <- CorrPrior(priorC.df[1,1], priorC.df[2,1], 
                              priorC.df[1,2], priorC.df[2,2])
    data.frame("Treatment"=c("Experimental", "Control"), 
               "Correlation"=c(prior_corr_e, prior_corr_c))
  })
  
  
  
#### 3. Tab posterior ####
  # Plot posterior treatment distributions 
  output$PlotPost = renderPlot({
    if(is.null(input$hot_pe)){priorE.df <- PriorE}
    else{priorE.df <- as.data.frame(hot_to_r(input$hot_pe))}
    if(is.null(input$hot_pc)){priorC.df <- PriorC}
    else{priorC.df <- as.data.frame(hot_to_r(input$hot_pc))}
    if(is.null(input$hot_xe)){dataE.df <- DataE}
    else{dataE.df <- as.data.frame(hot_to_r(input$hot_xe))}
    if(is.null(input$hot_xc)){dataC.df <- DataC}
    else{dataC.df <- as.data.frame(hot_to_r(input$hot_xc))}
    
    PlotPost(dataE.df[1,1],dataE.df[2,1],dataE.df[1,2],dataE.df[2,2],
             dataC.df[1,1],dataC.df[2,1],dataC.df[1,2],dataC.df[2,2],
             priorE.df[1,1],priorE.df[2,1],priorE.df[1,2],priorE.df[2,2],
             priorC.df[1,1],priorC.df[2,1],priorC.df[1,2],priorC.df[2,2],
             diff=FALSE)
  },
  width=400, height=650)  
  
 # Posterior correlation
  output$Corr_post <- renderTable({
    if(is.null(input$hot_pe)){priorE.df <- PriorE}
    else{priorE.df <- as.data.frame(hot_to_r(input$hot_pe))}
    if(is.null(input$hot_pc)){priorC.df <- PriorC}
    else{priorC.df <- as.data.frame(hot_to_r(input$hot_pc))}
    if(is.null(input$hot_xe)){dataE.df <- DataE}
    else{dataE.df <- as.data.frame(hot_to_r(input$hot_xe))}
    if(is.null(input$hot_xc)){dataC.df <- DataC}
    else{dataC.df <- as.data.frame(hot_to_r(input$hot_xc))}
    
    post_corr_e <- CorrPost(dataE.df[1,1],dataE.df[2,1],dataE.df[1,2],dataE.df[2,2],
                            priorE.df[1,1],priorE.df[2,1],priorE.df[1,2],priorE.df[2,2])
    post_corr_c <- CorrPost(dataC.df[1,1],dataC.df[2,1],dataC.df[1,2],dataC.df[2,2],
                            priorC.df[1,1],priorC.df[2,1],priorC.df[1,2],priorC.df[2,2])
    data.frame("Treatment"=c("Experimental", "Control"), 
               "Correlation"=c(post_corr_e, post_corr_c))
  })
  
#### 4. Tab posterior difference ####
   # Plot distribution treatment difference
  output$PlotDiff = renderPlot({
    if(is.null(input$hot_pe)){priorE.df <- PriorE}
    else{priorE.df <- as.data.frame(hot_to_r(input$hot_pe))}
    if(is.null(input$hot_pc)){priorC.df <- PriorC}
    else{priorC.df <- as.data.frame(hot_to_r(input$hot_pc))}
    if(is.null(input$hot_xe)){dataE.df <- DataE}
    else{dataE.df <- as.data.frame(hot_to_r(input$hot_xe))}
    if(is.null(input$hot_xc)){dataC.df <- DataC}
    else{dataC.df <- as.data.frame(hot_to_r(input$hot_xc))}
    
    PlotPost(dataE.df[1,1],dataE.df[2,1],dataE.df[1,2],dataE.df[2,2],
             dataC.df[1,1],dataC.df[2,1],dataC.df[1,2],dataC.df[2,2],
             priorE.df[1,1],priorE.df[2,1],priorE.df[1,2],priorE.df[2,2],
             priorC.df[1,1],priorC.df[2,1],priorC.df[1,2],priorC.df[2,2],
             diff=TRUE)
  }, width=400, height=400)
  
  # Posterior probabilities of treatment difference
  output$pop <- renderTable({
    if(is.null(input$hot_pe)){priorE.df <- PriorE}
    else{priorE.df <- as.data.frame(hot_to_r(input$hot_pe))}
    if(is.null(input$hot_pc)){priorC.df <- PriorC}
    else{priorC.df <- as.data.frame(hot_to_r(input$hot_pc))}
    if(is.null(input$hot_xe)){dataE.df <- DataE}
    else{dataE.df <- as.data.frame(hot_to_r(input$hot_xe))}
    if(is.null(input$hot_xc)){dataC.df <- DataC}
    else{dataC.df <- as.data.frame(hot_to_r(input$hot_xc))}
   
    POP <- PoP(dataE.df[1,1],dataE.df[2,1],dataE.df[1,2],dataE.df[2,2],
               dataC.df[1,1],dataC.df[2,1],dataC.df[1,2],dataC.df[2,2],
               priorE.df[1,1],priorE.df[2,1],priorE.df[1,2],priorE.df[2,2],
               priorC.df[1,1],priorC.df[2,1],priorC.df[1,2],priorC.df[2,2])
    
    data.frame("Decision rule"=c("Outcome 1: ", "Outcome 2: ",
                                 "All: ", "Any: ", "Compensatory: "), 
               "Posterior probability"=POP)
  })
  

  
  #### 5. Tab user guide ####
output$freq_table <- renderPlot({
  par(mar=c(0.1,0.1,0.1,0.1), fin=c(4.35, 2), pin=c(4.31, 1.95))
  plot(NULL, xlim=c(0,1), ylim=c(0.5,1), frame.plot=FALSE, xaxt="n", yaxt="n")
#  clip(x1=0, x2=1, y1= 0.5, y2=1)
  segments(x0=0, x1=1, y0=0.50, y1=0.50)
  segments(x0=0, x1=1, y0=0.80, y1=0.80)
  segments(x0=0, x1=1, y0=0.90, y1=0.90)
  segments(x0=0.25, x1=0.25, y0=0, y1=0)
  text(x=0.375, y=0.75, cex=1.25, labels=expression(x["11"]))
  text(x=0.375, y=0.65, cex=1.25, labels=expression(x["10"]))
  text(x=0.625, y=0.75, cex=1.25, labels=expression(x["01"]))
  text(x=0.625, y=0.65, cex=1.25, labels=expression(x["00"]))
  text(x=0.375, y=0.55, cex=1.25, labels=expression(x["1"]))
  text(x=0.625, y=0.55, cex=1.25, labels=expression({n-x["1"]}))
  text(x=0.875, y=0.75, cex=1.25, labels=expression(x["2"]))
  text(x=0.875, y=0.65, cex=1.25, labels=expression({n-x["2"]}))
  text(x=0.875, y=0.55, cex=1.25, labels=expression(n))
  text(x=0.375, y=0.85, labels="Successes")
  text(x=0.625, y=0.85, labels="Failures")
  text(x=0.02, y=0.75, pos=4, labels="Successes")
  text(x=0.02, y=0.65, pos=4, labels="Failures")
  text(x=0.02, y=0.55, pos=4,labels="Total")
  text(x=0.875, y=0.85, labels="Total")
  text(x=0.25, y=0.95, pos=4, labels="Outcome 1")
  text(x=0, y=0.85, pos=4, labels="Outcome 2")
  },
  width=400, height=150)


output$text_data <- renderUI({
    withMathJax("$\\Large(x_{11}^{obs}, x_{10}^{obs}, x_{01}^{obs}, x_{00}^{obs})$")
  })

output$text_prior <- renderUI({
  withMathJax("$\\Large(x_{11}^{prior}, x_{10}^{prior}, x_{01}^{prior}, x_{00}^{prior})$")
})

})
  
#### User interface ####
 ui <-  shinyUI(fixedPage(
   withMathJax(),
   # section below allows in-line LaTeX via $ in mathjax. Replace less-than-sign with < 
   # and grater-than-sign with >
   tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                });
                </script>
                ")),
#### 1. Main panel for displaying outputs ####
      #mainPanel(
navbarPage(title=" ",
           #### Tab 1. Home ####
           tabPanel("Home",
                    fixedRow(h2("Bayesian Analysis of Multiple Binary Outcomes")),
                    fixedRow(h6("Version 0.0.1, created by Xynthia Kavelaars")),
                    fixedRow(h3(" ")),
                    fixedRow(h4("Disclaimer")),
                    fixedRow(helpText("")),
                    fixedRow(h5(" ")),
                    fixedRow(h4("Goal")),
                    fixedRow(helpText("This Shiny App is designed to help users analyze data from multiple binary outcomes in a Bayesian way.
                                      Users are asked to enter observed cell frequencies and specify their plausible prior observations.
rmally distributed priors for regression coefficients.
                                      The data is based on data described in <<INSERT REF HERE>>")),
                    fixedRow(column(width=12,  plotOutput("freq_table", height="150px"))),
                    fixedRow(h3("")),
                    fixedRow(h3("Output per tab")),
                    fixedRow(h4(" ")),
                    fixedRow(h4("Data")),
                    fixedRow(column(width=12, align="left", offset=0.25,
                                    "- Sums of total successes and failures for each outcome.")),
                    fixedRow(column(width=12, align="left", offset=0.25,
                                    "- Observed correlations between outcomes")),
                    fixedRow(h4(" ")),
                    fixedRow(h4("Prior")),
                    fixedRow(column(width=12, align="left", offset=0.25, 
                                    "- Sums of total successes and failures for each outcome.")),
                    fixedRow(column(width=12, align="left", offset=0.25,
                                    "- Prior correlations between outcomes")),
                    fixedRow(h4(" ")),
                    fixedRow(h4("Treatment distributions")),
                    fixedRow(column(width=12, align="left", offset=0.25, 
                                    "- Plot of posterior treatment distributions.")),
                    fixedRow(column(width=12, align="left", offset=0.25,
                                    "- Posterior correlations between outcomes.")),
                    fixedRow(h4(" ")),
                    fixedRow(h4("Treatment differences")),
                    fixedRow(column(width=12, align="left", 
                                    "- Plot of posterior treatment difference between experimental and control treatment.")),
                    fixedRow(column(width=12, align="left", 
                                    "- Posterior probabilities of the following decision rules:")),
                    fixedRow(column(width=1, " "),
                             column(width=2, "Single outcome 1:"),
                             column(width=9, "The experimental treatment performs better on outcome 1")),
                    fixedRow(column(width=1, " "),
                             column(width=2, offset=0.25, "Single outcome 2:"),
                             column(width=9, "The experimental treatment performs better on outcome 2")), 
                    fixedRow(column(width=1, " "),
                             column(width=2, offset=0.25, "All:"),
                             column(width=9, "The experimental treatment performs better on both outcomes")),
                    fixedRow(column(width=1, " "),
                             column(width=2, offset=0.25, "Any:"),
                             column(width=9, "The experimental treatment performs better on outcome 1, outcome 2, or both outcomes")), 
                    fixedRow(column(width=1, " "),
                             column(width=2, offset=0.25, "Compensatory:"), 
                             column(width=9, "The experimental treatment performs better on the sum of outcome 1 and 2")),
                    fixedRow(h3(" ")),
                    fixedRow(h3("Refere")),
                    fixedRow("<<INSERT REFERENCE TO BOOK CHAPTER HERE>>"),
                    fixedRow(h3(" "))
                    ),
           
           #### Tab 2. User guide ####
           tabPanel("User guide",
                    fixedRow(h2("Decision-making with multiple binary outcomes")),
                    fixedRow(h3(" ")),
                    fixedRow(h3("Input per tab")),
                    fixedRow(h4("Data")),
                    fixedRow(helpText("Enter the observed response frequencies $\\Large(x_{11}^{obs}, x_{10}^{obs}, x_{01}^{obs}, x_{00}^{obs})$
                                      for the experimental and control treatments in the Data tab as presented in the table below. $\\Large x_{1}^{obs}, x_{2}^{obs}, n-x_{1}^{obs}, n-x_{2}^{obs}$ and $\\Large n^{obs}$ are computed automatically and cannot be entered.")),
                    fixedRow(h5(" ")),
                    fixedRow(h4("Prior")),
                    fixedRow(helpText("Enter the prior response frequencies $\\Large(x_{11}^{prior}, x_{10}^{prior}, x_{01}^{prior}, x_{00}^{prior})$ for the experimental and control treatments in th Prior tab as presented in the table below. 
                                      $\\Large x_{1}^{prior}, x_{2}^{prior}, n-x_{1}^{prior}, n-x_{2}^{prior}$ and $\\Large n^{prior}$ are computed automatically and cannot be entered. Default is Jeffrey's prior with 0.5 observation per response category.")),
                    fixedRow(column(width=12,  plotOutput("freq_table", height="150px"))),
                    fixedRow(h3("")),
                    fixedRow(h3("Output per tab")),
                    fixedRow(h4(" ")),
                    fixedRow(h4("Data")),
                    fixedRow(column(width=12, align="left", offset=0.25,
                                    "- Sums of total successes and failures for each outcome.")),
                    fixedRow(column(width=12, align="left", offset=0.25,
                                    "- Observed correlations between outcomes")),
                    fixedRow(h4(" ")),
                    fixedRow(h4("Prior")),
                    fixedRow(column(width=12, align="left", offset=0.25, 
                                    "- Sums of total successes and failures for each outcome.")),
                    fixedRow(column(width=12, align="left", offset=0.25,
                                    "- Prior correlations between outcomes")),
                    fixedRow(h4(" ")),
                    fixedRow(h4("Treatment distributions")),
                    fixedRow(column(width=12, align="left", offset=0.25, 
                                    "- Plot of posterior treatment distributions.")),
                    fixedRow(column(width=12, align="left", offset=0.25,
                                    "- Posterior correlations between outcomes.")),
                    fixedRow(h4(" ")),
                    fixedRow(h4("Treatment differences")),
                    fixedRow(column(width=12, align="left", 
                                    "- Plot of posterior treatment difference between experimental and control treatment.")),
                    fixedRow(column(width=12, align="left", 
                                    "- Posterior probabilities of the following decision rules:")),
                    fixedRow(column(width=1, " "),
                             column(width=2, "Single outcome 1:"),
                             column(width=9, "The experimental treatment performs better on outcome 1")),
                    fixedRow(column(width=1, " "),
                             column(width=2, offset=0.25, "Single outcome 2:"),
                             column(width=9, "The experimental treatment performs better on outcome 2")), 
                    fixedRow(column(width=1, " "),
                             column(width=2, offset=0.25, "All:"),
                             column(width=9, "The experimental treatment performs better on both outcomes")),
                    fixedRow(column(width=1, " "),
                             column(width=2, offset=0.25, "Any:"),
                             column(width=9, "The experimental treatment performs better on outcome 1, outcome 2, or both outcomes")), 
                    fixedRow(column(width=1, " "),
                             column(width=2, offset=0.25, "Compensatory:"), 
                             column(width=9, "The experimental treatment performs better on the sum of outcome 1 and 2")),
                    fixedRow(h3(" ")),
                    fixedRow(h3("Disclaimer")),
                    fixedRow(""),
                    fixedRow(h3(" "))
                    ),
           #### Tab 2. Data ####
  tabPanel("Data",
           fixedRow(""),
           fixedRow(column(width=12, align="center", h3("Observed frequencies"))),
           fixedRow(h3("")),
           fixedRow(
             column(width = 6, align="center", h4("Experimental treatment")),
             column(width = 6, align="center", h4("Control treatment"))),
           fixedRow(
             column(width = 5, rHandsontableOutput("hot_xe", width="100%")),
             column(width = 2, " "),
             column(width = 5, rHandsontableOutput("hot_xc", width="100%"))),
           fixedRow(h1("")),
           fixedRow(column(width=12, align="center", h3("Observed correlations"))),
           fixedRow(""),
           fixedRow(column(width=12, align="center", tableOutput("Corr_data"))),
           fixedRow("")
           ),

#### Tab 3. Prior distributions ####       
 tabPanel("Prior",
          fixedRow(""),
          fixedRow(column(width=12, align="center", h3("Prior frequencies"))), 
          fixedRow(h3("")),
          fixedRow(column(width = 6, align="center", h4("Experimental treatment")),
                   column(width = 6, align="center", h4("Control treatment"))),
          fixedRow(column(width = 5, rHandsontableOutput("hot_pe")),
                   column(width = 2),
                   column(width = 5, rHandsontableOutput("hot_pc"))),
          fixedRow(h1("")), 
          fixedRow(column(width=12, align="center", h3("Prior correlations"))),
          fixedRow(""),
          fixedRow(column(width=12, align="center",tableOutput("Corr_prior"))),
          fixedRow("")
          ),
                             
                    
                             
#### Tab 4. Posterior Treatment distributions ####
                             tabPanel("Treatment distributions",
                                      fixedRow(column(width=12, align="center", h3("Posterior distributions"))),
                                      fixedRow(column(width=12, align="center", plotOutput("PlotPost", inline=TRUE))),
                                      fixedRow(h1("")),
                                      fixedRow(column(width=12, align="center", h3("Posterior correlations"))),
                                      fixedRow(""),
                                      fixedRow(""),
                                      fixedRow(""),
                                      fixedRow(column(width=12, align="center", tableOutput("Corr_post")))
                             ),
#### Tab 5. Treatment difference ####
 tabPanel("Treatment difference",
  fixedRow(column(width=12, align="center", h3("Treatment difference distribution"))),
  fixedRow(column(width=12, align="center", plotOutput("PlotDiff", inline=TRUE))),
  fixedRow(h1(" ")),
  fixedRow(column(width=12, align="center", h3("Posterior probabilities"))),
  fixedRow(column(width=12, align="center", tableOutput("pop")))
  )
                             



#### Close user interface ----         
                           
      ,
      fluid=FALSE)))
  
  
  
  
  
  
 
shinyApp(ui=ui, server=server)
