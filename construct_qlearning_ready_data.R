foo_action <- function(x){
  v <- x
  v[is.na(x)] <- 0
  v <- 1*(v > 0)
  return(v)
}

construct_qlearning_ready_data <- function(){
  
  # keep nb rows
  n <- nrow(tmp_)
  
  # learning dataset (tm was for transition matrix, but surely a better name would be suited)
  tm <- data.frame(
    'daytime_0' = tmp_$daytime_disc[-n],
    'daytime_1' = tmp_$daytime_disc[-1],
    'glucose_0' = tmp_$value_blood_glucose_measurement_[-n],
    'glucose_1' = tmp_$value_blood_glucose_measurement_[-1]
  )
  tm$action_nph_0 <- foo_action(tmp_$value_nph_insulin_dose_[-n])
  tm$action_reg_0 <- foo_action(tmp_$value_regular_insulin_dose_[-n])
  tm$action_ult_0 <- foo_action(tmp_$value_ultralente_insulin_dose_[-n])
  tm$action_0 <- apply(cbind(tm[,c('action_nph_0', 'action_reg_0', 'action_ult_0')]), 1, paste0, collapse = "")
  tm$action_0 <- factor(x = tm$action_0)
  ### tm$action_P <- tm$action_0[]
  is_ok <- (tmp_$id[-1] == tmp_$id[-n])
  tm <- tm[is_ok,]
  tm <- na.omit(tm)
  return(tm)
  
}
