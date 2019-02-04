set.seed(1234)
ss <- simulate_dataset(N = 2000, Ti = 200, reco_mat = m4$qsa)
na_indexes <- which(rowSums(is.na(ss))>0)
Np <- c(10, 25, 50)
gamma_ <- 0.80
result_wm0 <- list()
result_wm1 <- list()
nbe <- 10
k <- 0
for(i in 3){
  siz <- Np[i]
  for(j in 1:nbe){
    k <- k+1
    in_samp <- floor(seq(1, nrow(ss), length.out = 2000))
    in_samp <- sample(x = in_samp[-length(in_samp)], size = siz)
    in_samp <- unlist(sapply(in_samp, function(i) i+1:Ti, simplify = FALSE))
    ss_samp <- ss[in_samp,]
    wm0 <- compute_world_model(tm_samp = ss_samp, index_reward_on_next_state = TRUE)
    wm1 <- compute_kernelized_world_model(tm_samp = ss_samp, use_kernel_for_rew = FALSE)
    dp1_wm0 <- dynamic_programming(tme = wm0$tme, rew = wm0$rew, adj = wm0$adj)
    dp1_wm1 <- dynamic_programming(tme = wm1$tme, rew = wm1$rew, adj = wm1$adj)
    eval_wm0 <- evaluate_reco(nb_experiments = 100, length_traj = 20, reco_mat = dp1_wm0[[1]])
    eval_wm1 <- evaluate_reco(nb_experiments = 100, length_traj = 20, reco_mat = dp1_wm1[[1]])
    result_wm0[[k]] <- list('siz' = siz, 'qsa' = dp1_wm0[[1]], 'vec' = eval_wm0)
    result_wm1[[k]] <- list('siz' = siz, 'qsa' = dp1_wm1[[1]], 'vec' = eval_wm1)
  }
}

k <- length(result_wm0)+1
ss_samp <- ss
wm0 <- compute_world_model(tm_samp = ss_samp, index_reward_on_next_state = TRUE)
wm1 <- compute_kernelized_world_model(tm_samp = ss_samp, use_kernel_for_rew = FALSE)
dp1_wm0 <- dynamic_programming(tme = wm0$tme, rew = wm0$rew, adj = wm0$adj)
dp1_wm1 <- dynamic_programming(tme = wm1$tme, rew = wm1$rew, adj = wm1$adj)
eval_wm0 <- evaluate_reco(nb_experiments = 100, length_traj = 20, reco_mat = dp1_wm0[[1]])
eval_wm1 <- evaluate_reco(nb_experiments = 100, length_traj = 20, reco_mat = dp1_wm1[[1]])
result_wm0[[k]] <- list('siz' = nrow(ss), 'qsa' = dp1_wm0[[1]], 'vec' = eval_wm0)
result_wm1[[k]] <- list('siz' = nrow(ss), 'qsa' = dp1_wm1[[1]], 'vec' = eval_wm1)

res_wm0 <- t(sapply(result_wm0, function(x) c(x$siz, mean(x$vec), median(x$vec), var(x$vec))))
res_wm0 <- data.frame(res_wm0)
colnames(res_wm0) <- c('siz', 'mean', 'median', 'var')
res_wm0$model <- 'classic'

res_wm1 <- t(sapply(result_wm1, function(x) c(x$siz, mean(x$vec), median(x$vec), var(x$vec))))
res_wm1 <- data.frame(res_wm1)
colnames(res_wm1) <- c('siz', 'mean', 'median', 'var')
res_wm1$model <- 'kernel'

resres <- rbind.data.frame(res_wm0, res_wm1)
resres$siz[resres$siz==398000] <- 2000
ggplot(data = resres) +
  geom_boxplot(mapping = aes(x=factor(model), y=mean, fill=factor(siz))) +
  ggtitle(label = 'Cumulative discounted rewards on simulated data, \n following previously optimized strategies') +
  xlab('Estimation method') +
  ylab('Cumulative\n discounted\n reward') +
  labs('fill' = 'Sample sizes') +
  theme(
    text = element_text(size = 10),
    title = element_text(size = 15),
    axis.title.x = element_text(size = 15, margin = margin(m,m,m,m)),
    axis.title.y = element_text(size = 15, hjust = 1, angle = 0)
  )

# rem2 <- get_state(ss_samp$daytime_0, ss_samp$glucose_0)
# rem <- apply(ss_samp[,c("daytime_0", "glucose_0")], 1, function(x) get_state(daytime_0 = x[1], glucose_0 = x[2]))
# rem2[998]

# evaluate a policy encoded in a vector with state to single action assigment
evaluate_reco <- function(reco_mat, nb_experiments, length_traj){
  v <- replicate(nb_experiments, {
    simu <- simulate_dataset(N = 1, Ti = length_traj, reco_mat = reco_mat, e_greedy = FALSE)
    simu$reward <- reward(simu$glucose_1)
    sum((gamma_)^{1:(length_traj-1) - 1} * simu$reward)
  })
}

##
complete_wm0 <- compute_world_model(tm_samp = ss, index_reward_on_next_state = TRUE)
image(as.matrix(complete_wm0$tme[,-(1:2)]))
image(as.matrix(wm0$tme[,-(1:2)]))
image(as.matrix(wm1$tme[,-(1:2)]))


wm0 <- compute_world_model(tm_samp = ss, index_reward_on_next_state = TRUE)
wm1 <- compute_kernelized_world_model(tm_samp = ss, use_kernel_for_rew = FALSE)
image(as.matrix(wm0$tme[,-(1:2)]))
image(as.matrix(wm1$tme[,-(1:2)]))
dp1_wm0 <- dynamic_programming(tme = wm0$tme, rew = wm0$rew, adj = wm0$adj)
dp1_wm1 <- dynamic_programming(tme = wm1$tme, rew = wm1$rew, adj = wm1$adj)
eval_wm0 <- evaluate_reco(nb_experiments = 100, length_traj = 20, reco_mat = dp1_wm0[[1]])
eval_wm1 <- evaluate_reco(nb_experiments = 100, length_traj = 20, reco_mat = dp1_wm1[[1]])

