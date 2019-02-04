reward <- function(g){
  r <- 1*(g >= 80 & g <= 140)
  return(r)
}

reward <- function(g){
  r <- exp(- (g - 100)^4 / (50)^4)
  return(r)
}

reward <- function(g){
  # taking ra close to 0, you get a quasi-uniform reward, equal to 1
  # taking ra very large, you get a quasi 0-1 reward, depending on whether you are within [a;b]
  a <- 80
  b <- 120
  ra <- 3
  r <- exp(- (g - a)^2 / (a/ra)^2) * (g < a) + 
    exp(- (g - b)^2 / (b/ra)^2) * (g > b) + 
    1 * (g >= a & g <= b)
  return(r)
}

# reward <- function(g){
#   r <- ((g-100)^2 / 100^2) * exp(- (abs(g - 100))^2 / (100)^2)
#   return(r)
# }

if(FALSE){ # to prevent to run when sourcing
  x <- seq(0, 300, 1)
  y <- reward(g = x)
  plot(x = x, y = y, type = 'o', pch = 20, cex = 0.5)
  for(u in c(60, 100, 140)){
    segments(x0 = u, y0 = 0, y1 = reward(g = u), lty = 2)
  }
  dfxy <- data.frame('x' = x, 'y' = y)
  m <- 20
  ggplot(data = dfxy, mapping = aes(x=x, y=y)) + 
    geom_line() + 
    geom_segment(aes(x=80, xend=80, y=0, yend=reward(80)), lty=2) +
    geom_segment(aes(x=120, xend=120, y=0, yend=reward(120)), lty=2) +
    xlab('Glucose Level') +
    ylab('Reward associated') +
    ggtitle('Reward level \n as a function of registered glucose.') +
    theme(
      title = element_text(size = 20, margin = margin(m,m,m,m)), 
      plot.margin = margin(m, m, m, m),
      axis.title.x = element_text(size = 20, margin = margin(m,m,m,m)),
      axis.title.y = element_text(size = 20, margin = margin(m,m,m,m), hjust = 1/2, angle = 90)
    )
  
} 

