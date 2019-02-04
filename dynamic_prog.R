# (implementation of algorithm 3) computing the state-action value functions for some policy pi in input
compute_qsa_based_on_pi <- function(pi_mat, tme, rew, adj){
  K <- 1e04
  # initialize qsa and vs for each moment order
  qsa <- array(data = 0, dim = c(nb_states_, nb_actions_))
  vs <- array(data = 0, dim = c(nb_states_))
  # boolean controlling the stopping of algo
  has_not_converged <- TRUE
  ##
  k <- 0
  while(has_not_converged){
    qsa_old <- qsa
    k = k + 1
    for(i in 1:nb_states_){
      vs[i] <- sum( pi_mat[i,] * qsa[i,] )
    }
    for(i in 1:nb_states_){
          for(a in 1:nb_actions_){
            if(adj[i,a] == 1){
              j <- which(tme$s0 == i & tme$a0 == a)
              ptr <- tme[j,-(1:2)]
              i_prim <- 1:nb_states_
              qsa[i,a] <- sum(ptr * (rew$r0[rew$s0==i & rew$a0==a] + gamma_ * vs[i_prim])) / sum(ptr)
            }
          }
      }
    has_not_converged <- ( max(qsa - qsa_old) > 1e-10 )
    if(k > K) has_not_converged <- FALSE
  }
  # return result
  return(qsa)
}

# (implementation of algorithm 4) computing the optimal action for each state based on evaluation of f
compute_optimalpi_based_on_qsa <- function(qsa_mat){
  #
  pi_mat <- array(data = 0, dim = c(nb_states_, nb_actions_))
  vsa <- array(data = 0, dim = c(nb_states_, nb_actions_))
  #
  for(i in 1:nb_states_){
      eval_i <- numeric(length = nb_actions_)
      for(a in 1:nb_actions_){
        if(adj[i,a]==1){
          eval_i[a] <- qsa_mat[i,a]
        }else{
          eval_i[a] <- (-1e10)
        }
      }
      vsa[i,] <- (eval_i)
      ai_star <- which.max(eval_i)
      pi_mat[i,ai_star] <- 1
  }
  #
  return( list('pi' = pi_mat, 'vsa' = vsa ) ) 
}

# initialize random policy
init_pi <- function(adj){
  #
  pi_mat <- array(data = 0, dim = c(nb_states_, nb_actions_))
  #
  for(i in 1:nb_states_){
    a_possible <- which(adj[i,] == 1)
    pi_mat[i, a_possible] <- 1/length(a_possible)
  }
  return(pi_mat)
}


dynamic_programming_2 <- function(){

  sta0 <- get_state(daytime_0 = tm_$daytime_0, glucose_0 = tm_$glucose_0)
  sta1 <- get_state(daytime_0 = tm_$daytime_1, glucose_0 = tm_$glucose_1)
  act0 <- tm_$action_0
  rew0 <- reward(tm_$glucose_0)
  mdp_learnt <- learn_mdp_model(reward_indep_next_state = TRUE, state = sta0, new_state = sta1, action = act0, reward = rew0, initial_model = NULL)
  ### m_test <- td_lookup(state = sta0, action = act0, new_state = sta1, reward = rew0)
  ## mdp_learnt$reward
  ns <-  max(mdp_learnt$reward$s0)
  na <- max(mdp_learnt$reward$a0)

  # policy
  softmax <- FALSE
  epsilon_par <- 0.0 # 0.05
  tau <- 1
  
  # vector of state value function v(s)
  vs <- rep(0, ns)

  # probability matrix p(a | s)
  psa <- matrix(data = rep(1/na, ns*na), nrow = ns, ncol = na)
  psa[1,] <- psa[ns,] <- 0
  
  # transition probability P(s' | s, a)
  tsa <- mdp_learnt$transition
  
  # reward function
  reward_vec <- mdp_learnt$reward
  
  # is state terminal ?
  isTerminal <- rep(FALSE, ns)
  
  # gamma parameter
  ## gamma <- 0.5

#
j <- 1
not_stopping <- TRUE
while(not_stopping){
  
  # save previous state value functions
  pvs <- vs
  
  # iterate over states
  for(i in 1:ns){
    if(!isTerminal[i]){
      # update state values for each probable output
      vs[i] <- 0
      for(k in 1:na){
        new_state <- tsa[tsa$s0 == i & tsa$a0 == k,]
        if(TRUE){
          vs[i] <- vs[i] + psa[i,k]*(sum(new_state$prob * (wm0$rew$r0[wm0$rew$s0==i][1] + gamma_  * pvs[new_state$s1])))
        }else{
          vs[i] <- vs[i] + psa[i,k]*(reward_vec$mean_value[reward_vec$s0 == i & reward_vec$a0 == k] + gamma_*(new_state$prob %*% pvs[new_state$s1]))
        }
      }
    }
  }
  
  # update policy
  for(i in 1:ns){
    for(k in 1:na){
      new_state <- tsa[tsa$s0 == i & tsa$a0 == k,]
      if(TRUE){
        psa[i,k] <- sum(new_state$prob * (wm0$rew$r0[wm0$rew$s0==i][1] + gamma_  * pvs[new_state$s1]))
      }else{
        psa[i,k] <- reward_vec$mean_value[reward_vec$s0 == i & reward_vec$a0 == k] + gamma_ * (new_state$prob %*% pvs[new_state$s1]) 
      }
    }
    if(softmax){
      # first possibility : softmax
      psa[i,] <- exp(psa[i,] / tau)
      psa[i,] <- psa[i,] / sum(psa[i,])
    }else{
      # second : \epsilon-greedy policy
      max_index <- which(psa[i,] == max(psa[i,]))
      psa[i,] <- 0
      if(runif(n = 1) < epsilon_par){
        max_index <- sample(x = 1:na, size = 1)
      }else{
        if(length(max_index)>1) max_index <- sample(x = max_index, size = 1)
      }
      psa[i, max_index] <- 1
    }
  }
  
  # print value functions
  ## cat('\n') cat(round(vs, 1))
  not_stopping <- (max(abs(pvs - vs)) > 1e-10)
  if(j > 20000) not_stopping <- FALSE
  
  # counter update
  j <- j + 1
  
}

qsa <- array(data = 0, dim = c(ns, na))
for(i in 1:ns){
  for(k in 1:na){
    new_state <- tsa[tsa$s0 == i & tsa$a0 == k,]
    if(TRUE){
      qsa[i,k] <- sum(new_state$prob * (wm0$rew$r0[wm0$rew$s0==i][1] + gamma_  * pvs[new_state$s1]))
    }else{
      qsa[i,k] <- reward_vec$mean_value[reward_vec$s0 == i & reward_vec$a0 == k] + gamma_ * (new_state$prob %*% pvs[new_state$s1]) 
    }
  }
}

return(list(psa, qsa, j))

}


# 
dynamic_programming <- function(tme, rew, adj){
  # solve optimization problem
  pi_mat <- init_pi(adj)
  qsa_mat <- compute_qsa_based_on_pi(pi_mat, tme , rew, adj)
  has_not_converged <- TRUE
  count <- 0
  while(has_not_converged){
    ## action_weights <- colSums(pi_mat) * as.numeric(rowSums(adj) %*% adj)
    old_qsa_mat <- qsa_mat
    ## temp <- compute_optimalpi_based_on_qsa(qsa_mat)
    temp <- apply(qsa_mat, 1, greedy)
    temp <- t(sapply(temp, function(x){ v <- rep(0,nb_actions_) ; v[x] <- 1 ; v}))
    pi_mat <- temp
    qsa_mat <- compute_qsa_based_on_pi(pi_mat, tme, rew, adj)
    has_not_converged <- ( max(qsa_mat - old_qsa_mat) > 1e-04 )
    count <- count + 1
  }
  return(list(qsa_mat, pi_mat))
}
## image(t(prop.table(old_qsa_mat, 1)))

