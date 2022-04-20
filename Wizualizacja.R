library(plotly)


x <- c('Cvine<br>Optimal','Cvine<br>Report',
       'Dvine<br>Optimal', 'Dvine<br>Report',
       'Cb<br>Optimal', 'Cb<br>Report',
       'MixCb<br>Optimal','MixCb<br>Report',
       'Cort<br>Optimal','Cort<br>Report')
y <- c(-2, -1.8, -1.5,-1.3,-1,-0.7, -0.5, 0, 0.5, 1)
base <- c(-1.47,0.18,-1.55,0.17,-1.48,0.26,-1.48,0.29,-1.54,0.27)
revenue <- c(1.89,0.14,1.96,0.15,2.07,0.24,2.1,0.23,2.16,0.26)
costs <- c(0,0,0,0, 0, 0, -1, -1.2, -0.4, 0)
profit <- c(0,0,0,0, 0, 0, 0, 0, 0, 0.6)
data <- data.frame(x, base, revenue, costs, profit)

#The default order will be alphabetized unless specified as below:
data$x <- factor(data$x, levels = data[["x"]])
fig <- plot_ly(data, x = ~x, y = ~base,type = 'bar', marker = list(color = 'rgba(1,1,1, 0.0)'))
fig
fig <- fig %>% add_trace(y = ~revenue, marker = list(color = 'rgba(55, 128, 191, 0.7)',
                                                     line = list(color = 'rgba(55, 128, 191, 0.7)',
                                                                 width = 2)))

width = 2)))
fig
fig<-fig%>%add_text(x = ~x,
                    y= base, text =~as.character(base)
)
fig

fig<-fig%>%add_text(x = ~x,
                    y= base+revenue, text =~as.character(base+revenue)
)
fig


fig <- fig %>% layout(title = '',
                      xaxis = list(title = ""),
                      yaxis = list(title = "Diversification effect"),
                      autosize = TRUE,
                      barmode = 'stack',
                      paper_bgcolor = 'rgba(245, 246, 249, 1)',
                      plot_bgcolor = 'rgba(245, 246, 249, 1)',
                      showlegend = FALSE)
fig
fig <- fig %>% add_annotations(text = c(''),
                               x = x,
                               y = y,
                               xref = "x",
                               yref = "y",
                               font = list(family = 'Arial',
                                           size = 14,
                                           color = 'rgba(245, 246, 249, 1)'),
                               showarrow = FALSE)

fig

