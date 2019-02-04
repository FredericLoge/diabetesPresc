## TODO: 
## - ADD TIMER ESTIMATE
## - MASK FOR CONSTRAINED ACTION

#' @title Batch Q-Learning
#' @param my_simulate_batch simulator
#' @param nb_states
#' @param nb_actions
#' @param prior_over_states
#' @param gamma = 0.9, 
#' @param batch_size = 100, 
#' @param qsa_max_abs_diff_stopping = 1e-03, 
#' @param max_nb_runs = 2000,
#' @param  do_speedy_qlearning = TRUE, 
#' @param alpha_k_indexed_on_s_a = TRUE, 
#' @param omega = 1,
#' @param sampling_policy c('uniform', 'e_greedy_10pct', 'greedy', 'ucb')
batch_qlearning <- function(
  my_simulate_batch = function(N, reco, prior_init){ NULL },
  initial_qsa = array(data = 0, dim = c(nb_states_, nb_actions_)),
  prior_over_states,
  gamma = 0.9, batch_size = 100, 
  qsa_max_abs_diff_stopping = 1e-03, max_nb_runs = 2000,
  do_speedy_qlearning = TRUE, alpha_k_indexed_on_s_a = TRUE, omega = 1,
  sampling_policy = c('uniform', 'e_greedy_10pct', 'greedy', 'ucb')){
  
  # keep first element 
  sampling_policy <- sampling_policy[1]
  
  # initialize Q(s,a)
  qsa_old <- qsa <- nsa <- initial_qsa
  
  # store average rewards obtained
  avg_reward <- numeric(length = max_nb_runs)
  
  ## core algorithm
  k <- 1 # counter
  convergence_not_attained <- TRUE
  while(convergence_not_attained){
    
    # store old Q(s,a) matrices
    qsa_very_old <- qsa_old
    qsa_old <- qsa
    
    # choose next actions via sub-optimal policy
    ## current_reco <- get_reco_from_qtable(...)
    if(sampling_policy == 'uniform'){
      current_reco = base::sample(x = 1:nb_actions_, replace = TRUE, size = nb_states_)
    }else if(sampling_policy == 'greedy'){
      current_reco = apply(qsa, 1, greedy)
    }else if(sampling_policy == 'e_greedy_10pct'){
      current_reco = apply(qsa, 1, e_greedy, e = 0.1)
    }else if(sampling_policy == 'e_greedy_20pct'){
      current_reco = apply(qsa, 1, e_greedy, e = 0.2)
    }else if(sampling_policy == 'ucb'){
      current_reco = sapply(1:nb_states_, function(s) ucb(q = qsa[s,], n = nsa[s,]))
    }else{
      stop('Sorry, declared policy was not recognized.')
    }
    
    # store average reward under selected policy, weighted by prior over states
    if(!is.null(prior_over_states)){
      avg_reward[k] <- sum(prior_states * sapply(1:nb_states_, function(s) qsa[s,current_reco[s]])) / sum(prior_states)
    }else{
      avg_reward[k] <- mean(sapply(1:nb_states_, function(s) qsa[s,current_reco[s]]))
    }
    
    # simulate new batch of data
    simu <- my_simulate_batch(N = batch_size, reco = current_reco, prior_init = prior_over_states)
    
    # get state and rewards
    simu$state <- NA
    for(i in 1:nrow(simu)){
      simu$state[i] <- get_state(daytime_0 = simu$daytime_0[i], glucose_0 = simu$glucose_0[i])
    }
    simu$reward <- c(reward(simu$glucose_0[-1]), NA)
    
    # set learning rate, \alpha_k
    alpha_k <- 1/(k^omega)
    
    # for every transition:
    for(i in 1:(nrow(simu)-1)){
      
      # get (state, action, reward, state)
      s0 <- simu$state[i]
      a0 <- as.numeric(simu$action_0[i])
      r0 <- simu$reward[i]
      s1 <- simu$state[i+1]
      
      # add counts to nsa
      nsa[s0, a0] <- nsa[s0, a0] + 1 
      
      #
      if(alpha_k_indexed_on_s_a){
        alpha_k <- 1 / (nsa[s0,a0]+1)
      }
      
      # update rules
      if(do_speedy_qlearning){
        # update rule for Speedy Q-Learning
        update_0 <- r0 + gamma * max(qsa_very_old[s1,])
        update_1 <- r0 + gamma * max(qsa_old[s1,])
        qsa[s0,a0] <- qsa_old[s0,a0] + alpha_k * (update_0 - qsa_old[s0,a0]) + (1 - alpha_k) * (update_1 - update_0)
      }else{
        # update rule for regular Q-Learning
        qsa[s0,a0] <- qsa_old[s0,a0] + alpha_k * (r0 + gamma * max(qsa_old[s1,]) - qsa_old[s0,a0])
      }
      
    }
    
    # step counter
    k <- k + 1
    
    # stopping conditions check
    convergence_not_attained <- (max(abs(qsa - qsa_old)) > qsa_max_abs_diff_stopping)
    if(k > max_nb_runs) convergence_not_attained <- FALSE
    
  }
  
  # return convergence status, algorithm controls, qsa table
  l <- list('convergence' = !convergence_not_attained, 
            'algo_ctrl' = list(prior_over_states,
                               gamma, batch_size,
                               qsa_max_abs_diff_stopping, max_nb_runs,
                               do_speedy_qlearning, alpha_k_indexed_on_s_a, omega,
                               sampling_policy),
            'qsa' = qsa,
            'nsa' = nsa,
            'avg_reward' = avg_reward[1:k])
  return(l)
  
}

