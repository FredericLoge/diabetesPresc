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

  # state x action x state transition kernel
  ktme <- expand.grid('s0' = 1:nb_states_, 'a0' = 1:nb_actions_)
  temppp <- array(0, dim = c(nrow(ktme), nb_states_))
  for(i in 1:nrow(ktme)){
    
    k1 <- k2 <- k3 <- rep(NA, nrow(tm_samp))
    
    # evaluate first kernel
    state_ref <- state_discretization_grid_[ktme$s0[i],]
    ind <- as.numeric(state_ref$glucose_0_discretized)
    mu <- (cuts_[ind] + cuts_[ind+1])/2
    ratio <- (10 * (mu - tm_samp$glucose_0) / (mu))^2
    k1 <- exp( - ratio) * (tm_$daytime_0 == tm_)
    
    # evaluate second kernel
    temp <- levels_action_[ktme$a0[i]]
    if(TRUE){
      k2 <- 1*(tm_samp$action_0 == temp)
    }else{
      temp <- as.numeric(strsplit(x = temp, split = '')[[1]])
      wei <- rep(1, 3)
      wei <- wei / sum(wei)
      k2 <- wei[1] * (tm_samp$action_nph_0 == temp[1]) + wei[2] * (tm_samp$action_reg_0 == temp[2]) + wei[3] * (tm_samp$action_ult_0 == temp[3])
    }
    
    cum_sum_k3 <- 0
    for(s1_index in 1:nb_states_){
      
      # evaluate third kernel
      state1_ref <- state_discretization_grid_[s1_index,]
      ind <- as.numeric(state1_ref$glucose_0_discretized)
      mu <- (cuts_[ind] + cuts_[ind+1])/2
      ratio <- (10 * (mu - tm_samp$glucose_1) / (mu))^2
      k3 <- exp( - ratio)
      k3 <- k3 * ((as.numeric(tm_samp$daytime_0) - as.numeric(tm_samp$daytime_1)) %in% c(+1,-1))
      
      # store kernel value
      temppp[i,s1_index] <- sum(k1*k2*k3) / sum(k1*k2)
      cum_sum_k3 <- cum_sum_k3 + sum(k3)
    }
    
    temppp[i,] <- temppp[i,]/sum(temppp[i,])
    
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

# estimate world model
wm0 <- compute_world_model(index_reward_on_next_state = TRUE)
wm1 <- compute_kernelized_world_model(use_kernel_for_rew = FALSE)

#
image(x = t(as.matrix(wm0$tme[,-(1:2)]))>0)
sum(t(as.matrix(wm0$tme[,-(1:2)]))>0, na.rm=TRUE)
sum(is.na((as.matrix(wm0$tme[,-(1:2)]))))
dim(wm0$tme[,-(1:2)])
# look at differences in transition matrices
image(x = t(as.matrix(wm0$tme[,-(1:2)])))
image(x = (as.matrix(wm1$tme[,-(1:2)])))
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

estim <- data.frame(
  state_discretization_grid_[,2],
  round(wm0$rew$r0[wm0$rew$s0 == 1 & wm0$rew$a0 == 1],3),
  round(wm1$rew$r0[wm1$rew$s0 == 1 & wm1$rew$a0 == 1],3))
unique(estim)
cat(paste0(apply(unique(estim), 1 , paste0, collapse = '  &  '), collapse = '\\\\ \n'))
#
softmax <- function(x, beta){
  prop.table(exp(x*beta), 1)
}

# myopic optimization
gamma_ <- 0
dp0_wm0 <- dynamic_programming(tme = wm0$tme, rew = wm0$rew, adj = wm0$adj)
dp0_wm1 <- dynamic_programming(tme = wm1$tme, rew = wm1$rew, adj = wm1$adj)
apply(dp0_wm0[[2]], 1, greedy)
apply(dp0_wm1[[2]], 1, greedy)
plot_my_heatmap(mat = softmax(x = dp0_wm0[[1]], beta = 10), r = rn, c = NULL, te = dp0_wm0[[1]])
plot_my_heatmap(mat = softmax(x = dp0_wm1[[1]], beta = 10), r = rn, c = NULL, te = dp0_wm1[[1]])

# in this myopic criterion, it seems that action 3 is brought into the light thanks to
# kernelized estimates, which detect its good 1-step performance

# mid-term optimization
gamma_ <- 0.80
dp1_wm0 <- dynamic_programming(tme = wm0$tme, rew = wm0$rew, adj = wm0$adj)
dp1_wm1 <- dynamic_programming(tme = wm1$tme, rew = wm1$rew, adj = wm1$adj)
apply(dp1_wm0[[2]], 1, greedy)
apply(dp1_wm1[[2]], 1, greedy)
plot_my_heatmap(mat = dp1_wm0[[1]], r = NULL, c = NULL, te = dp1_wm0[[1]])
plot_my_heatmap(mat = dp1_wm1[[1]], r = NULL, c = NULL, te = dp1_wm1[[1]])
plot_my_heatmap(mat = softmax(x = dp1_wm0[[1]], beta = 10), r = NULL, c = NULL, te = dp1_wm0[[1]])
plot_my_heatmap(mat = softmax(x = dp1_wm1[[1]], beta = 10), r = NULL, c = NULL, te = dp1_wm1[[1]])

newdf <- rbind(dp1_wm0[[1]], dp1_wm1[[1]])
# normalize ?
newdf <- softmax(newdf, beta = 5)
#
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

# comparatives of recommendations given
summary(c(abs(softmax(x = dp1_wm0[[1]], beta = 10) - softmax(x = dp0_wm0[[1]], beta = 10))))
summary(c(abs(softmax(x = dp1_wm1[[1]], beta = 10) - softmax(x = dp0_wm1[[1]], beta = 10))))
summary(c(abs(softmax(x = dp0_wm0[[1]], beta = 10) - softmax(x = dp0_wm1[[1]], beta = 10))))
summary(c(abs(softmax(x = dp1_wm0[[1]], beta = 10) - softmax(x = dp1_wm1[[1]], beta = 10))))

#
m3 <- dynamic_programming_2()
cbind(m3[[2]], dp1_wm0[[1]])
cbind(m3[[1]], dp1_wm0[[2]])
table(apply(m3[[1]], 1, greedy), apply(dp1_wm0[[2]], 1, greedy))

#
library(dplyr)
wm0$tme %>% filter(a0 == 1 & s0 == 1)
wm0$rew %>% filter(a0 == 1 & s0 == 1)
mdp_learnt$transition %>% filter(a0 == 1 & s0 == 1)
mdp_learnt$reward %>% filter(a0 == 1 & s0 == 1)

# transition_matrix_estimate_ <- ktme
# state_action_adjacency_ <- array(data = 1, dim = c(nb_states_, nb_actions_))
# index_not_ok <- rowSums(is.nan(as.matrix(ktme))) > 0
# index_not_ok <- ktme[index_not_ok, c(1,2)]
# for(i in 1:nrow(index_not_ok)){
#   state_action_adjacency_[index_not_ok[i,1], index_not_ok[i,2]] <- 0
# }
dp1 <- dynamic_programming()
head(dp0[[1]])
head(dp1[[1]])
image(t(prop.table(dp0[[1]], 1)))
image(t(prop.table(dp1[[1]], 1)))
matt <- dp1[[1]]
for(i in 1:nrow(matt)){
  matt[index_not_ok[i,1], index_not_ok[i,2]] <- NA
}
plot_my_heatmap(mat = matt, r = rn, c = NULL, te = dp1[[1]])
plot_my_heatmap(mat = prop.table(exp(10*dp1[[1]]),1), r = rn, c = NULL, te = dp1[[1]])
plot_my_heatmap(mat = dp0[[1]], r = rn, c = NULL, te = dp0[[1]])
plot_my_heatmap(mat = prop.table(exp(10*dp0[[1]]),1), r = rn, c = NULL, te = dp0[[1]])
apply(dp0[[2]], 1, function(x) which(x==1))
apply(dp1[[2]], 1, function(x) which(x==1))

# approach 1. construct non-parametric model estimates and apply Dynammic Programming

# construct Directed Acyclic Graph (dag) from tm_ data 
dag_ <- construct_dag_from_data()

# TODO add dag representation


## TODO : add a simulator from category of GLUCOSE to continuous GLUCOSE.

## TODO : constrain the action space conditionally to the state information
## gros probleme : couples etats-action rares
## idee: contraindre l'espace d'etat pour apprendre nos controles

# for graphical representations
rn <- with(state_discretization_grid_, paste0('D: ', daytime_0, '; G: ', glucose_0_discretized))
cn <- NULL ## levels_action_

# model 0: empirical policy
m0 <- get_empirical_policy()
plot_my_heatmap(mat = prop.table(m0, 1), r = NULL, c = NULL, te = m0)

newdf <- prop.table(m0, 1)
colnames(newdf) <- paste0('a', 1:5)
rownames(newdf) <- paste0('s', 1:20)
newdf <- as.data.frame(newdf)

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

# model 1: speedy batch Q-Learning, sampling policy: uniform
m1 <- batch_qlearning(
  my_simulate_batch = function(N, reco, prior_init){
    simulate_batch(N = N, reco = reco, prior_init = prior_init)
  },
  prior_over_states = NULL,
  gamma = 0.8, batch_size = 10, 
  qsa_max_abs_diff_stopping = 1e-03, max_nb_runs = 20000,
  do_speedy_qlearning = TRUE, alpha_k_indexed_on_s_a = TRUE, omega = 1,
  sampling_policy = c('uniform', 'e_greedy_10pct', 'greedy', 'ucb'))

# average reward behavior throughout epochs
plot.ts(m1$avg_reward, main = 'Average return per iteration')
plot_my_heatmap(mat = m1$qsa, r = rn, c = NULL, te = m1$qsa)
softmax_qsa <- prop.table(exp(10*m1$qsa), 1)
plot_my_heatmap(mat = softmax_qsa, r = rn, c = NULL, te = m1$qsa)

# model 2: batch Q-Learning, sampling policy: uniform
m2 <- batch_qlearning(
  my_simulate_batch = function(N, reco, prior_init){
    simulate_batch(N = N, reco = reco, prior_init = prior_init)
  },
  prior_over_states = NULL,
  gamma = 0.9, batch_size = 10, 
  qsa_max_abs_diff_stopping = 1e-03, max_nb_runs = 2000,
  do_speedy_qlearning = FALSE, alpha_k_indexed_on_s_a = TRUE, omega = 0.6,
  sampling_policy = c('uniform', 'e_greedy_10pct', 'greedy', 'ucb'))

# average reward behavior throughout epochs
plot.ts(m2$avg_reward, main = 'Average return per iteration')
plot_my_heatmap(mat = m2$qsa, r = rn, c = NULL, te = m2$qsa)
softmax_qsa <- prop.table(exp(10*m2$qsa), 1)
plot_my_heatmap(mat = softmax_qsa, r = rn, c = NULL, te = m2$qsa)

### Comparing results from m1 and m2
### both results are coherent, so that's a first good news !
### it means Speedy just speeds it up, did not get us to 
### a completely different solution

### Let's look at other sampling policies

# model 3: batch Q-Learning, sampling policy: greedy
m3 <- batch_qlearning(
  my_simulate_batch = function(N, reco, prior_init){
    simulate_batch(N = N, reco = reco, prior_init = prior_init)
  }, initial_qsa = m1$qsa,
  prior_over_states = NULL,
  gamma = 0.8, batch_size = 10, 
  qsa_max_abs_diff_stopping = 1e-04, max_nb_runs = 20000,
  do_speedy_qlearning = TRUE, alpha_k_indexed_on_s_a = TRUE, omega = 0.6,
  sampling_policy = 'greedy')
m3uw <- batch_qlearning(
  my_simulate_batch = function(N, reco, prior_init){
    simulate_batch(N = N, reco = reco, prior_init = prior_init)
  }, initial_qsa = m1$qsa,
  prior_over_states = NULL,
  gamma = 0.9, batch_size = 10, 
  qsa_max_abs_diff_stopping = 1e-04, max_nb_runs = 20000,
  do_speedy_qlearning = TRUE, alpha_k_indexed_on_s_a = TRUE, omega = 0.6,
  sampling_policy = 'greedy')

# m3_mem <- list()
m3_mem[[length(m3_mem)+1]] <- m3$qsa
## m3_mem[[length(m3_mem)]] <- NULL

i <- i + 1
plot_my_heatmap(mat = m4_mem[[i]], r = rn, c = cn, te = '')


plot.ts(m1$avg_reward[-length(m1$avg_reward)], ylim = c(0,6), xlim = c(0, 5000))
lines(m2$avg_reward[-length(m2$avg_reward)], col = 'darkblue')
lines(m3$avg_reward[-length(m3$avg_reward)], col = 'lightblue')
lines(m3uw$avg_reward[-length(m3uw$avg_reward)], col = 'blue')

nb_experiments <- 100
length_traj <- 80
gamma <- 0.9
m3_reco <- apply(m3$qsa, 1, greedy)
m3_eval <- evaluate_reco(reco = m3_reco, nb_experiments = nb_experiments, length_traj = length_traj, gamma = gamma)
m3uw_reco <- apply(m3uw$qsa, 1, greedy)
m3uw_eval <- evaluate_reco(reco = m3uw_reco, nb_experiments = nb_experiments, length_traj = length_traj, gamma = gamma)
summary(m3_eval)
summary(m3uw_eval)

m4_reco <- apply(m4$qsa, 1, e_greedy, e = 0.1)
m4_eval <- evaluate_reco(reco = m4_reco, nb_experiments = nb_experiments, length_traj = length_traj, gamma = gamma)
m5_reco <- sapply(1:nb_states_, function(s) ucb(q = m5$qsa[s,], n = m5$nsa[s,]))
m5_eval <- evaluate_reco(reco = m5_reco, nb_experiments = nb_experiments, length_traj = length_traj, gamma = gamma)

# model 4: batch Q-Learning, sampling policy: greedy
m4 <- batch_qlearning(
  my_simulate_batch = function(N, reco, prior_init){
    simulate_batch(N = N, reco = reco, prior_init = prior_init)
  }, initial_qsa = m4$qsa,
  prior_over_states = NULL,
  gamma = 0.8, batch_size = 30, 
  qsa_max_abs_diff_stopping = 1e-04, max_nb_runs = 20000,
  do_speedy_qlearning = TRUE, alpha_k_indexed_on_s_a = TRUE, omega = 0.6,
  sampling_policy = 'e_greedy_10pct')

# m4_mem <- list()
m4_mem[[length(m4_mem)+1]] <- m4$qsa
## m4_mem[[length(m4_mem)]] <- NULL

i = 0
i <- i + 1
plot_my_heatmap(mat = softmax(m4_mem[[i]], beta=5), r = NULL, c = NULL, te = '')

# model 5: batch Q-Learning, sampling policy: greedy
m5 <- batch_qlearning(
  my_simulate_batch = function(N, reco, prior_init){
    simulate_batch(N = N, reco = reco, prior_init = prior_init)
  },
  prior_over_states = NULL,
  gamma = 0.9, batch_size = 10, 
  qsa_max_abs_diff_stopping = 1e-06, max_nb_runs = 20000,
  do_speedy_qlearning = TRUE, alpha_k_indexed_on_s_a = TRUE, omega = 0.7,
  sampling_policy = 'ucb')

# model 6: batch Q-Learning, e-greedy: 20\%
m6 <- batch_qlearning(
  my_simulate_batch = function(N, reco, prior_init){
    simulate_batch(N = N, reco = reco, prior_init = prior_init)
  },
  prior_over_states = NULL,
  gamma = 0.9, batch_size = 10, 
  qsa_max_abs_diff_stopping = 1e-7, max_nb_runs = 2000,
  do_speedy_qlearning = TRUE, alpha_k_indexed_on_s_a = TRUE, omega = 0.7,
  sampling_policy = 'e_greedy_20pct')
m7 <- batch_qlearning(
  my_simulate_batch = function(N, reco, prior_init){
    simulate_batch(N = N, reco = reco, prior_init = prior_init)
  }, 
  initial_qsa = m6$qsa, 
  prior_over_states = NULL,
  gamma = 0.9, batch_size = 100, 
  qsa_max_abs_diff_stopping = 1e-10, max_nb_runs = 2000,
  do_speedy_qlearning = TRUE, alpha_k_indexed_on_s_a = TRUE, omega = 0.7,
  sampling_policy = 'greedy')


# average reward behavior throughout epochs
plot.ts(m3$avg_reward, main = 'Average return per iteration')
plot.ts(m4$avg_reward, col = 'darkblue')
plot.ts(m5$avg_reward, col = 'blue')
plot.ts(m6$avg_reward)
plot.ts(m7$avg_reward, col = 'red')

# very different recommendations
softmat_coef <- 1
plot_my_heatmap(mat = prop.table(exp(softmat_coef*m1$qsa), 1), r = rn, c = NULL, te = m1$qsa)
plot_my_heatmap(mat = prop.table(exp(softmat_coef*m2$qsa), 1), r = rn, c = NULL, te = m2$qsa)
plot_my_heatmap(mat = prop.table(exp(softmat_coef*m3$qsa), 1), r = rn, c = NULL, te = m3$qsa)
plot_my_heatmap(mat = prop.table(exp(softmat_coef*m4$qsa), 1), r = rn, c = NULL, te = m4$qsa)
plot_my_heatmap(mat = prop.table(exp(softmat_coef*m6$qsa), 1), r = rn, c = NULL, te = m6$qsa)
plot_my_heatmap(mat = prop.table(exp(softmat_coef*m5$qsa), 1), r = rn, c = NULL, te = m5$qsa)

plot_my_heatmap(mat = m3$qsa, r = rn, c = NULL, te = m3$qsa)
plot_my_heatmap(mat = m4$qsa, r = rn, c = NULL, te = m4$qsa)

plot_my_heatmap(mat = m6$qsa, r = rn, c = NULL, te = m6$qsa)
plot_my_heatmap(mat = m7$qsa, r = rn, c = NULL, te = m7$qsa)


# # evaluate policy followed on real data
# get_empirical_policy <- function(){
#   # recover states throughout dataset
#   s0 <- get_state(daytime_0 = tm_$daytime_0, glucose_0 = tm_$glucose_0)
#   s0 <- factor(s0, levels = 1:nb_states_)
#   # cross states and actions from available dataset
#   a0 <- factor(tm_$action_0, levels = levels_action_)
#   # recover states throughout dataset
#   s1 <- c(s0[-1],NA)
#   # return Q matrix
#   prop.table(table(s1[a0 == '000'], s0[a0 == '000']))
#   return(nsa)
# }

#
nb_experiments <- 100
length_traj <- 80
gamma <- 0.9
m3_reco <- apply(m3$qsa, 1, greedy)
m3_eval <- evaluate_reco(reco = m3_reco, nb_experiments = nb_experiments, length_traj = length_traj, gamma = gamma)
m4_reco <- apply(m4$qsa, 1, e_greedy, e = 0.1)
m4_eval <- evaluate_reco(reco = m4_reco, nb_experiments = nb_experiments, length_traj = length_traj, gamma = gamma)
m5_reco <- sapply(1:nb_states_, function(s) ucb(q = m5$qsa[s,], n = m5$nsa[s,]))
m5_eval <- evaluate_reco(reco = m5_reco, nb_experiments = nb_experiments, length_traj = length_traj, gamma = gamma)

#
cbind(m3_reco, m4_reco, m5_reco)

#
model_eval <- data.frame('model' = -1, 'value' = c(m3_eval, m4_eval, m5_eval))
model_eval$model <- rep(x = 3:5, times = nb_experiments)
ggplot(data = model_eval) +
  geom_violin(mapping = aes(x = factor(model), y = value))
ggplot(data = model_eval) +
  geom_histogram(mapping = aes(col = factor(model), x = value))
by(model_eval$value, model_eval$model, summary)

#
sim <- simulate_batch(N = 1000, reco = base::sample(x = 1:nb_actions_, replace = TRUE, size = nb_states_))
sim$reward <- reward(sim$glucose_0)
plot.ts(cumsum(0.9^{1:1000 - 1} * sim$reward))

# SYNCHRONOUS Q-LEARNING
#
#

# continuous sampler of glucose value, conditionally to the fact it belongs in some range
clustered_ref <- arules::discretize(x = tm_$glucose_0, method = 'fixed', categories = cuts)
from_discrete_to_conti <- function(glucose_0_discretized){
  base::sample(x = tm_$glucose_0[clustered_ref == glucose_0_discretized], size = 1)
}


## TODO NEXT:
## > METTRE UNE ACTION REALISTE / DONE!
## > CHOISIR PARAMETRES (ALPHA, GAMMA) ASSURANT LA CONVERGENCE / DONE!
## > ADD FUNCTIONAL REPRESENTATIONS
## > SPEEDY Q-LEARNING / DOES NOT WORK BETTER THAN Q-LEARNING ... ISN'T THAT WEIRD?
## > CREATE SIMULATOR WITH DESIGNATED STARTING STATE TO PROPERLY IMPLEMENT Q-LEARNING
## > EVALUATION OFFLINE DE POLICY SAUVEGARDEES
## > DEEP Q-NETWORK
## > PROPOSER D'AUTRES DAG PLUS COMPLEXES ? PEUT-ETRE PAS BESOIN SI LE MODELE INITIAL EST VRAI.