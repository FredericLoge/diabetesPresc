# select uniformly amongst maximizers
greedy <- function(x){
  v <- which(x == max(x))
  if(length(v)>1) v <- base::sample(x = v, size = 1)
  return(v)
}
greedy(x = c(0, 6, 2))
greedy(c(4,3,6))

# select uniformly amongst maximizers with proba 1-e, otherwise random
e_greedy <- function(x, e = 0.1){
  if(runif(n = 1) > e){
    greedy(x)
  }else{
    base::sample(x = 1:nb_actions_, size = 1)
  }
}

# select with UCB criterion
ucb <- function(q, n){
  y <- q + 2*log(n)/sqrt(0.01 + n)
  greedy(y)
}