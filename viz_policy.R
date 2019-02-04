# plot heatmap
plot_my_heatmap <- function(mat, r, c, te){
  require(plotly)
  m <- list(l = 300, r = 0, b = 100, t = 100)
  p <- plot_ly(
    y = r, x = c,
    text = te, 
    z = mat, type = "heatmap", hoverinfo = "text", width = 800, height = 800, colors=c("red","green")
  ) %>% layout(autosize = F, margin = m, title = 'Q(s,a) estimates')
  p
}