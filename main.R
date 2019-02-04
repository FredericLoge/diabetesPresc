# R packages, crucial
library(arules) # discretize
library(nnet) # multinomial

# R viz packages 
library(ggplot2)
library(shiny)
library(shinydashboard)
# library(dplyr)

# source R files - the order is indicative of their use within this .R file
source('descriptive_app.R') # A shiny app and some ggplot, basic stat analysis of the data
source('read_and_parse_data.R') # Prepare the Diabetes data, collected within file 'datasets/'
source('construct_qlearning_ready_data.R') # Prepare the data for a RL format (S, A, R, S, A)
source('construct_dag_from_data.R') # Construct a Directed Acyclic Graph (dag) based on precedent data (which includes samplers)
source('simulators.R') # Data simulators built upon the dag
source('policies.R') # Implementing classic policies : Greedy, e-Greedy, UCB  
source('qlearning_and_variants.R') # QLearning algorithm and its variants (Batch) (Speedy) Q-Learning ; ideally would include Zap Q-Learning
source('reward.R') # Defining an appropriate reward function
source('evaluate_policy.R') # Policy evaluation approach based on simulators
source('viz_policy.R') # Visualize proposed policies
source('dynamic_prog.R')

# extract function names
get_function_names_from_R_file <- function(fn = 'descriptive_app.R'){
  # read lines from R file
  re <- readLines(con = fn)  
  # identify lines with text '<- function('
  re <- re[grep(pattern = '<- function(', x = re, fixed = TRUE)]
  # extract function name
  re <- sub(pattern = '<- function\\(.*', replacement = '', x = re)
  return(re)
}
get_function_names_from_R_file('construct_qlearning_ready_data.R')

# convention : global variables will have their names ending with '_'
# list of global variables throughout file:
#   tmp_; tm_;  dag_; levels_daytime_;  levels_action_; cuts_;  state_discretization_grid_ 

# load dataset
tmp_ <- read_and_parse_diabetes_data()

# GOTO descriptive_app.R for descriptive analytics
# GOTO "source('descriptive_analysis.R')" for nice plots 

# construct dataset appropriate for Q-Learning
tm_ <- construct_qlearning_ready_data()

# checking on action variable
table(tm_$action_0)

# removing rare action instances
tm_ <- tm_[tm_$action_0 %in% names(which(table(tm_$action_0) > 100)),]
tm_$action_0 <- factor(x = tm_$action_0, levels = sort(unique(tm_$action_0)))

# keep
levels_daytime_ <- levels(tm_$daytime_0)
levels_action_ <- levels(tm_$action_0)
nb_actions_ <- length(levels_action_)

# creating state grid for glucose in order to build state representation
## cuts_ <- discretize(x = tm_$glucose_0, method = 'cluster', categories = 4, onlycuts = TRUE)
cuts_ <- c(0, 50, 80, 120, 200, max(tm_$glucose_0)*1.01)
state_discretization_grid_ <- expand.grid(
  'daytime_0' = levels_daytime_, 
  'glucose_0_discretized' = levels(discretize(x = tm_$glucose_0, method = 'fixed', categories = cuts_)))
nb_states_ <- nrow(state_discretization_grid_)

# identify state ids
get_state <- function(daytime_0, glucose_0){
  glucose_0_discretized <- discretize(x = glucose_0, method = "fixed", categories = cuts_)
  to_match_1 <- paste0(daytime_0, glucose_0_discretized)
  to_match_2 <- apply(state_discretization_grid_, 1, paste0, collapse = '')
  to_match_1 <- as.numeric(factor(x = to_match_1, levels = to_match_2))
  return(to_match_1)
}
levels_glucose_ <- unique(state_discretization_grid_$glucose_0_discretized)
get_state_2 <- function(daytime_0, glucose_0){
  glucose_0_index <- which(glucose_0 < cuts_)[1] - 1
  index <- which(state_discretization_grid_$daytime_0 == daytime_0 &
              state_discretization_grid_$glucose_0_discretized == levels_glucose_[glucose_0_index])
  return(index)
}

#' @title Compute MDP elements from global dataset tm_, based on pre-specified Kernels
compute_kernelized_world_model <- function(tm_samp = tm_, index_reward_on_next_state = TRUE, use_kernel_for_rew = FALSE){
  
  # compute discretized versions
  state_tm <- get_state(daytime_0 = tm_samp$daytime_0, glucose_0 = tm_samp$glucose_0)
  state_tm_tp1 <- c(state_tm[-1], 1)
  
  # state x action x state transition kernel
  ktme <- expand.grid('s0' = 1:nb_states_, 'a0' = 1:nb_actions_)
  temppp <- array(0, dim = c(nrow(ktme), nb_states_))
  
  #
  for(i in 1:nrow(ktme)){
    s0_index <- ktme$s0[i]
    a0_index <- ktme$a0[i]
    for(j in 1:nb_states_){
      s1_index <- j
      Tt <- state_discretization_grid_$daytime_0[s0_index]
      Ttp1 <- state_discretization_grid_$daytime_0[s1_index]
      Ot <- state_discretization_grid_$glucose_0_discretized[s0_index]
      Otp1 <- state_discretization_grid_$glucose_0_discretized[s1_index]
      At <- levels_action_[a0_index]
      PROBA_1 <- sum(tm_samp$daytime_0 == Tt & tm_samp$daytime_1 == Ttp1) / sum(tm_samp$daytime_0 == Tt)
      PROBA_2 <- sum(state_discretization_grid_$glucose_0_discretized[state_tm] == Ot & 
                       state_discretization_grid_$glucose_0_discretized[state_tm_tp1] == Otp1 & 
                       tm_samp$action_0 == At) / sum(state_discretization_grid_$glucose_0_discretized[state_tm] == Ot & tm_samp$action_0 == At)
      temppp[i,j] <- PROBA_1 * PROBA_2
    }
  }

  ktme <- cbind(ktme, 's1_prob' = temppp)
  
  # adjacency matrix: putting aside never observed state x action pairs
  adj <- array(data = 1, dim = c(nb_states_, nb_actions_))
  index_not_ok <- rowSums(is.nan(as.matrix(ktme))) > 0
  index_not_ok <- ktme[index_not_ok, c(1,2)]
  for(i in 1:nrow(index_not_ok)){
    adj[index_not_ok[i,1], index_not_ok[i,2]] <- 0
  }
  
  # expected reward for given (state, action) pair -> could be kernelized as well, but a priori not necessary
  rew <- expand.grid('s0' = 1:nb_states_, 'a0' = 1:nb_actions_, 's1' = 1:nb_states_)
  rew$r0 <- NA
  for(i in 1:nb_states_){
    if(use_kernel_for_rew){
      j <- as.numeric(state_discretization_grid_$glucose_0_discretized[i])
      mu <- (cuts_[j] + cuts_[j+1])/2
      ratio <- (10 * (mu - tm_samp$glucose_1) / (mu))^2
      k4 <- exp( - ratio)
      tj <- sum(k4 * reward(tm_samp$glucose_1)) / sum(k4)
    }else{
      j <- as.numeric(state_discretization_grid_$glucose_0_discretized[i])
      condi <- (tm_samp$glucose_0 >= cuts_[j] & tm_samp$glucose_0 <= cuts_[j+1])
      tj <- mean(reward(tm_samp$glucose_0[condi]))
    }
    if(index_reward_on_next_state == TRUE){
      ind <- which(rew$s1 == i)
    }else{
      ind <- which(rew$s0 == i)
    }
    rew$r0[ind] <- tj
  }
  
  # return elements
  l <- list('tme' = ktme, 'adj' = adj, 'rew' = rew)
  return(l)
  
}

#' @title Compute MDP elements from global dataset tm_, crossing all modalities of state variables
compute_world_model <- function(tm_samp = tm_, index_reward_on_next_state = TRUE){

  # state x action x state transition matrix
  tme <- expand.grid('s0' = 1:nb_states_, 'a0' = 1:nb_actions_)
  temp <- array(0, dim = c(nrow(tme), nb_states_))
  a0 <- as.numeric(tm_samp$action_0)
  s0 <- apply(tm_samp[,c("daytime_0", "glucose_0")], 1, function(x) get_state(daytime_0 = x[1], glucose_0 = x[2]))
  s1 <- get_state(daytime_0 = tm_samp$daytime_1, glucose_0 = tm_samp$glucose_1)
  ni <- numeric(nb_states_)
  for(i in 1:nrow(tme)){
    cond <- ((s0 == tme$s0[i]) & (a0 == tme$a0[i]))
    temp[i,] <- as.numeric(table(x = factor(x = s1[cond], levels = 1:nb_states_)))
    ni[i] <- sum(temp[i,])
    temp[i,] <- temp[i,] / sum(temp[i,])
  }
  tme <- cbind(tme, 's1_prob' = temp)

  # adjacency matrix: putting aside never observed state x action pairs
  adj <- array(data = 1, dim = c(nb_states_, nb_actions_))
  index_not_ok <- which(ni < 1)
  index_not_ok <- tme[index_not_ok, c(1,2)]
  for(i in 1:nrow(index_not_ok)){
    adj[index_not_ok[i,1], index_not_ok[i,2]] <- 0
  }
  
  # expected reward for given (state, action) pair
  rew <- expand.grid('s0' = 1:nb_states_, 'a0' = 1:nb_actions_, 's1' = 1:nb_states_)
  rew$r0 <- NA
  for(i in 1:nb_states_){
    j <- as.numeric(state_discretization_grid_$glucose_0_discretized[i])
    condi <- (tm_samp$glucose_0 >= cuts_[j] & tm_samp$glucose_0 <= cuts_[j+1])
    tj <- mean(reward(tm_samp$glucose_0[condi]))
    ## tj <- mean(reward(runif(n = 1e4, min = cuts_[j], max = cuts_[j+1])))
    if(index_reward_on_next_state == TRUE){
      ind <- which(rew$s1 == i)
    }else{
      ind <- which(rew$s0 == i)
    }
    rew$r0[ind] <- tj
  }
  
  # return such values
  return(list('tme' = tme, 'adj' = adj, 'rew' = rew))
  
}

###


# estimate world model
wm0 <- compute_world_model(index_reward_on_next_state = TRUE)
wm1 <- compute_kernelized_world_model(use_kernel_for_rew = FALSE)

# nice plot, comparing transition matrices
new_tme <- rbind.data.frame(wm0$tme[order(wm0$tme$s0, wm0$tme$a0),], wm1$tme)
new_tme$model <- rep(c('classic', 'kernel'), each = nrow(wm0$tme))
new_tme$id <- paste0('(', new_tme$s0, ', ', new_tme$a0, ')')
new_tme$id <- factor(new_tme$id, levels = new_tme$id[1:nrow(wm0$tme)])
new_tme <- data.table::melt(data = new_tme[,-(1:2)], id.vars = c('id', 'model'), value.name = 'value')
new_tme$variable <- factor(sub('s1_prob.', '', new_tme$variable, fixed = TRUE), levels = 1:nb_states_)
m <- 20
ggplot(data = new_tme) +
  geom_tile(mapping = aes(x = variable, y = id, fill = value)) +
  facet_grid(~model) +
  xlab(expression(S[t+1])) +
  ylab(expression(paste('(', S[t], ', ', A[t], ')'))) +
  ggtitle('Transition matrix estimates based on classic estimates and their Kernel versions\n') +
  labs(fill = 'Transition probability   ') +
  theme(
    text = element_text(size = 10),
    legend.position = 'bottom',
    legend.title = element_text(size = 20),
    title = element_text(size = 25, margin = margin(m,m,m,m), vjust = 0.5), 
    plot.margin = margin(m, m, m, m),
    axis.title.x = element_text(size = 20, margin = margin(m,m,m,m)),
    axis.title.y = element_text(size = 20, margin = margin(m,m,m,m), hjust = 1/2, angle = 0)
  )

# reward | S_1
estim <- data.frame(
  state_discretization_grid_[,2],
  round(wm0$rew$r0[wm0$rew$s0 == 1 & wm0$rew$a0 == 1],3),
  round(wm1$rew$r0[wm1$rew$s0 == 1 & wm1$rew$a0 == 1],3))
unique(estim)
cat(paste0(apply(unique(estim), 1 , paste0, collapse = '  &  '), collapse = '\\\\ \n'))
              
# softmax, for representation
softmax <- function(x, beta){
  prop.table(exp(x*beta), 1)
}

# compute optimal policies for both world estimates
gamma_ <- 0.80
dp1_wm0 <- dynamic_programming(tme = wm0$tme, rew = wm0$rew, adj = wm0$adj)
dp1_wm1 <- dynamic_programming(tme = wm1$tme, rew = wm1$rew, adj = wm1$adj)
apply(dp1_wm0[[2]], 1, greedy)
apply(dp1_wm1[[2]], 1, greedy)
plot_my_heatmap(mat = dp1_wm0[[1]], r = NULL, c = NULL, te = dp1_wm0[[1]])
plot_my_heatmap(mat = dp1_wm1[[1]], r = NULL, c = NULL, te = dp1_wm1[[1]])
plot_my_heatmap(mat = softmax(x = dp1_wm0[[1]], beta = 10), r = NULL, c = NULL, te = dp1_wm0[[1]])
plot_my_heatmap(mat = softmax(x = dp1_wm1[[1]], beta = 10), r = NULL, c = NULL, te = dp1_wm1[[1]])

# comparative plot of recommendations
newdf <- rbind(dp1_wm0[[1]], dp1_wm1[[1]])
newdf <- softmax(newdf, beta = 5)
newdf <- as.data.frame(newdf)
colnames(newdf) <- paste0('a', 1:5)
newdf$model <- rep(c('classic', 'kernel'), each = nb_states_)
newdf$state <- rep(paste0('s', 1:nb_states_), times = 2)
newdf$state <- factor(newdf$state, levels = paste0('s', 1:nb_states_))
newdf <- data.table::melt(data = newdf, id.vars = c('state', 'model'))
newdf$value[newdf$value == 0] <- NA
newdf$rescaled_value <- newdf$value
for(l in unique(newdf$model)){
  cond <- (newdf$model == l)
  min_value <- min(newdf$value[cond], na.rm = TRUE)
  max_value <- max(newdf$value[cond], na.rm = TRUE)
  newdf$rescaled_value[cond] <- (newdf$rescaled_value[cond] - min_value) / (max_value - min_value)
}
m <- 20
ggplot(data = newdf) +
  geom_tile(mapping = aes(x = variable, y = state, fill = value)) +
  facet_grid(~model) +
  xlab(expression(A[t])) +
  ylab(expression(S[t])) +
  ggtitle('State-Action function estimates - normalized\n') +
  labs(fill = expression(paste(exp(beta*hat(Q)[sa]), " /", sum(exp(beta*hat(Q)[si]), i==1, 20) ))) +
  scale_fill_continuous(low = 'white', high = 'darkblue', na.value = 'black') +
  theme(
    text = element_text(size = 10),
    legend.position = 'bottom',
    legend.title = element_text(size = 20),
    title = element_text(size = 25, margin = margin(m,m,m,m), vjust = 0.5), 
    plot.margin = margin(m, m, m, m),
    axis.title.x = element_text(size = 20, margin = margin(m,m,m,m)),
    axis.title.y = element_text(size = 20, margin = margin(m,m,m,m), hjust = 1/2, angle = 0)
  )
              
###

# empirical policy
m0 <- get_empirical_policy()
plot_my_heatmap(mat = prop.table(m0, 1), r = NULL, c = NULL, te = m0)
newdf <- prop.table(m0, 1)
colnames(newdf) <- paste0('a', 1:5)
rownames(newdf) <- paste0('s', 1:20)
newdf <- as.data.frame(newdf)

# represent empirical policy
m <- 20
ggplot(data = newdf) +
  geom_tile(mapping = aes(x = Var2, y = Var1, fill = Freq)) +
  xlab(expression(A[t])) +
  ylab(expression(S[t])) +
  ggtitle('State-Action function counts from real data \n') +
  labs(fill = expression(paste(N[sa], " /", sum(N[si], i==1, 20), '  '))) +
  scale_fill_continuous(low = 'white', high = 'darkblue', na.value = 'black') +
  theme(
    text = element_text(size = 10),
    legend.position = 'bottom',
    legend.title = element_text(size = 20),
    title = element_text(size = 25, margin = margin(m,m,m,m), vjust = 0.5), 
    plot.margin = margin(m, m, m, m),
    axis.title.x = element_text(size = 20, margin = margin(m,m,m,m)),
    axis.title.y = element_text(size = 20, margin = margin(m,m,m,m), hjust = 1/2, angle = 0)
  )

###

# construct Directed Acyclic Graph (dag) from tm_ data 
dag_ <- construct_dag_from_data()
