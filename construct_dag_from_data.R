# DAG is constructed according to the following relations :
# 
#  DAYTIME_{t+1} ~ DAYTIME_{t} (1)
#  ACTION_{t} ~ GLUCOSE_{t} + ACTION_{t} (2)
#  GLUCOSE_{t+1} | GLUCOSE_{t} + ACTION_{t} + DAYTIME_{t} (3)


#' @title building DAG (Directed Acyclic Graph) based on dataset
construct_dag_from_data <- function(){
  
  # construct dag
  dag <- list()
  
  # construct sample, arc 1 | multinomial distribution
  dag$fit_model_1 <- prop.table(table('daytime_0' = tm_$daytime_0, 'daytime_1' = tm_$daytime_1), 1)
  dag$sampler_1 <- function(daytime_0, prob_mat = dag$fit_model_1){
    # extract probability distribution
    index <- which(daytime_0 == levels_daytime_)
    prob <- as.numeric(prob_mat[index,])
    # sample over multinomial
    rand <- base::sample(x = levels_daytime_, size = 1, prob = prob)
    # return formatted version
    return(factor(x = rand, levels = levels_daytime_))
  }
  
  # construct sample, arc 2 | logistic model
  dag$fit_model_2 <- multinom(formula = action_0 ~ glucose_0 + daytime_0, data = tm_, family = 'multinomial', verbose = FALSE)
  dag$summary_model_2 <- summary(dag$fit_model_2)
  dag$sampler_2 <- function(glucose_0, daytime_0,
                            mc = dag$summary_model_2$coefficients,
                            si = NULL){
    # compute scalar product X\beta
    log_proba <- mc[,1] + mc[,2]*glucose_0
    index <- which(daytime_0 == levels_daytime_)
    if(index > 1) log_proba <- log_proba + mc[,2 - 1 + index]
    # add error and apply logit function with 0.5 threshold
    proba <- c(1, exp(as.numeric(log_proba)))
    rand <- base::sample(x = levels_action_, size = 1, prob = proba)
    # return formatted version
    return(rand)
  }
  
  # construct sample, arc 3 | 
  dag$fit_model_3 <- lm(formula = glucose_1 ~ glucose_0 + daytime_1 + action_0, data = tm_)
  dag$sampler_3 <- function(glucose_0, daytime_1, action_0, 
                            mc = as.numeric(dag$fit_model_3$coefficients), 
                            si = sqrt(mean((dag$fit_model_3$residuals)^2))){
    # compute scalar product X\beta
    rand <- mc[1] + mc[2]*glucose_0
    index_daytime <- which(daytime_1 == levels_daytime_)
    if(index_daytime > 1) rand <- rand + mc[2 + index_daytime - 1]
    index_action <- action_0
    if(length(action_0)==0) stop('length action 0 is 0')
    if(index_action > 1){
      rand <- rand + mc[2 + length(levels_daytime_)-1 + index_action - 1]
    }
    # add gaussian noise
    rand <- rand + rnorm(n = 1, mean = 0, sd = si)
    # return truncated observation
    return(min(500, max(30, rand)))
  }
  
  return(dag)
  
}