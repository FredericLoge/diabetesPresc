# reco_mat <- m4$qsa
# Ti <- 200
# ss <- simulate_dataset(N = 20, Ti = 100, reco_mat = m4$qsa)

simulate_dataset <- function(N, Ti, reco_mat, e_greedy = TRUE){
    
  global_simu <- list()
  for(i in 1:N){
    
    # choose next actions via sub-optimal policy
    if(e_greedy){
      current_reco = apply(reco_mat, 1, e_greedy, e = 0.30)      
    }else{
      current_reco = apply(reco_mat, 1, greedy)
    }
    
    # simulate new batch of data
    simu <- my_simulate_batch(N = Ti, reco = current_reco, prior_init = NULL)
    
    # add lagged variables
    simu$daytime_1 <- simu$daytime_0
    simu$daytime_1[-Ti] <- simu$daytime_0[-1]
    simu$glucose_1 <- c(simu$glucose_0[-1], NA)
    global_simu[[i]] <- simu[1:(Ti-1),]
    
  }
  
  global_simu <- do.call(rbind.data.frame, global_simu)
  return(global_simu)
  
}

#' @title simulate data batch
#' @param N : nb simulation steps
#' @param dag : DAG from which to sample
simulate_batch <- function(N, reco = NULL, prior_init = NULL){
  
  # if NULL, uniform recommendations
  if(is.null(reco)){
    base::sample(x = 1:length(levels_action_), replace = TRUE, size = nb_states_)
  }
  
  # result data frame
  res <- data.frame('time' = 1:N, 
                    'daytime_0' = factor(x = NA, levels = levels_daytime_), 
                    'glucose_0' = NA,
                    'action_0' = NA)
  
  # init
  ## TODO: use prior_init in this section
  index <- base::sample(x = 1:nrow(tm_), size = 1)
  d0 <- tm_$daytime_0[index] # factor(x = base::sample(x = levels_daytime_, size = 1), levels = levels_daytime_)
  g0 <- tm_$glucose_0[index] # min(500, max(30, rnorm(n = 1, mean = 100, sd = 30)))
  s0 <- get_state_2(daytime_0 = d0, glucose_0 = g0)
  a0 <- reco[s0]

  # run
  for(i in 1:N){
    # save values
    res$daytime_0[i] <- d0 ; res$glucose_0[i] <- g0 ; res$action_0[i] <- a0
    # generate daytime_1 | daytime_0
    d1 <- dag_$sampler_1(daytime_0 = d0)
    # generate glucose_1 | glucose_0, action_0, daytime_0
    g1 <- dag_$sampler_3(glucose_0 = g0, daytime_1 = d1, action_0 = a0)
    # generate action_0 | glucose_0, daytime_0
    s1 <- get_state_2(daytime_0 = d1, glucose_0 = g1)
    a1 <- reco[s1]
    # change containers
    d0 <- d1 ; g0 <- g1 ; a0 <- a1
  }
  #
  res$action_0 <- factor(x = levels_action_[res$action_0], levels = levels_action_)
  return(res)
}

### TESTS
# res <- simulate_batch(N = 1000)
# plot.ts(res[1:100,])
# table(res$daytime_0)

#' @title simulate next state based on current state 
#' @param daytime_0
#' @param glucose_0
#' @param action_0
simulate_one_observation_from <- function(daytime_0, glucose_0, action_0){
  # generate daytime_1 | daytime_0
  d1 <- dag_$sampler_1(daytime_0 = daytime_0)
  # generate glucose_1 | glucose_0, action_0, daytime_0
  g1 <- dag_$sampler_3(glucose_0 = glucose_0, daytime_1 = d1, action_0 = action_0)
  return(list('daytime_0' = d1, 'glucose_0' = g1))
}

#' @title simulate batch with deciced policy, stored in reco vector
#' @param N : nb simulation steps
#' @param reco : same size as the number of states, guides the decision-making
simulate_batch_w_reco <- function(N, reco){
  
  # result data frame
  res <- data.frame('time' = 1:N, 
                    'daytime_0' = factor(x = NA, levels = levels_daytime_), 
                    'glucose_0' = NA,
                    'action_0' = factor(x = NA, levels = levels_action_))
  
  # init
  d0 <- factor(x = sample(x = levels_daytime_, size = 1), levels = levels_daytime_)
  g0 <- min(500, max(30, rnorm(n = 1, mean = 100, sd = 30)))
  s0 <- get_state(daytime_0 = d0, glucose_0 = g0)
  a0 <- levels_action_[reco[s0]]
  
  # run
  for(i in 1:N){
    # save values
    res$daytime_0[i] <- d0 ; res$glucose_0[i] <- g0 ; res$action_0[i] <- a0
    # generate daytime_1 | daytime_0
    d1 <- dag_$sampler_1(daytime_0 = d0)
    # generate glucose_1 | glucose_0, action_0, daytime_0
    g1 <- dag_$sampler_3(glucose_0 = g0, daytime_1 = d1, action_0 = a0)
    # generate state and take action based on recommendations
    s1 <- get_state(daytime_0 = d1, glucose_0 = g1)
    a1 <- levels_action_[reco[s1]]
    # change containers
    d0 <- d1 ; g0 <- g1 ; a0 <- a1
  }
  #
  return(res)
}
### TESTS
# res <- simulate_batch_w_reco(N = 1000, reco = sample(x = 1:length(levels_action_), replace = TRUE, size = nrow(state_discretization_grid)))
# plot.ts(res[1:100,])