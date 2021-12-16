mirorenderer_resultOutput <- function(id, height = NULL, options = NULL, path = NULL){
    ns <- NS(id)
# set default height
    if(is.null(height)){
        height <- 400
    } 
    tagList( 
    #define rendererOutput function here
        fluidRow(
            #column(5,dataTableOutput(ns('dataT'))),
            column(12,plotOutput(ns('plot'),height=height))
        ),

        fluidRow(
            #column(5,dataTableOutput(ns('dataT'))),
            column(12,plotOutput(ns('plot2'),height=height))
        )

    ) 
}
 
 renderMirorenderer_result <- function(input, output, session, data, options = NULL, path = NULL, rendererEnv = NULL, views = NULL, outputScalarsFull = NULL, ...){
     # separate sir information from beta information
     dataSIR <- data[data$header=='sirDetail',]
     dataBeta <- data[data$header=='b',]

     S <-  dataSIR[dataSIR$state=='S',]$value
     I <-  dataSIR[dataSIR$state=='I',]$value
     R <-  dataSIR[dataSIR$state=='R',]$value
     
     beta <-  dataBeta[dataSIR$state=='R',]$value


     t <-  c(0)
     for (i in 1:4000){
         t <- c(t,i*0.01)
     }

     SIR <- data.frame(t=t, S=S, I=I, R=R)
     
     beta <- data.frame(t=t,beta=beta)
     #output$dataT<-DT::renderDataTable(datatable(SIR))

     output$plot <- renderPlot({
        plot(SIR$t,                              # Draw first time series
            SIR$S,
            type = "l",
            col = 2,
            main = "SIR population over time",
            ylim = c(0, 1),
            xlab = "time",
            ylab = "population ratio")
        lines(SIR$t,                             # Draw second time series
            SIR$I,
            type = "l",
            col = 3)
        lines(SIR$t,                             # Draw third time series
            SIR$R,
            type = "l",
            col = 4)
        legend("topright",                           # Add legend to plot
            c("S", "I", "R"),
            lty = 1,
            col = 2:4)

     })

     output$plot2 <- renderPlot({
        plot(beta$t,                              # Draw first time series
            beta$beta,
            type = "l",
            col = 2,
            main = "beta over time",
            ylim = c(0, 1),
            xlab = "time",
            ylab = "beta value")

     })
}
