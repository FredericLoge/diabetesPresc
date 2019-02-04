# evaluate policy followed on real data
get_empirical_policy <- function(){
  # recover states throughout dataset
  sta <- get_state(daytime_0 = tm_$daytime_0, glucose_0 = tm_$glucose_0)
  # cross states and actions from available dataset
  nsa <- table(factor(sta, levels = 1:nb_states_), factor(tm_$action_0, levels = levels_action_))
   # return Q matrix
  return(nsa)
}

# evaluate a completely random policy
evaluate_random <- function(nb_experiments, length_traj, gamma = 1){
  v <- replicate(nb_experiments, {
    reco <- base::sample(x = 1:length(levels_action_), size = nb_states_, replace = TRUE)
    simu <- simulate_batch_w_reco(N = length_traj, reco = reco)
    simu$reward <- reward(simu$glucose_0)
    sum(gamma^{1:length_traj - 1} * simu$reward)
  })
}

# evaluate a policy encoded in a vector with state to single action assigment
evaluate_reco <- function(reco, nb_experiments, length_traj, gamma = 1){
  v <- replicate(nb_experiments, {
    simu <- simulate_batch_w_reco(N = length_traj, reco = reco)
    simu$reward <- reward(simu$glucose_0)
    sum(gamma^{1:length_traj - 1} * simu$reward)
  })
}
