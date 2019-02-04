# blood glucose histogram
m <- 20
ggplot(data = tmp) +
  geom_histogram(aes(x = value_blood_glucose_measurement_)) +
  xlim(0, 500) + 
  xlab('Blood Glucose Level') +
  ylab('Count') +
  ggtitle('Blood Glucose Level Histogram') +
  theme(
    title = element_text(size = 25, margin = margin(m,m,m,m)), 
    plot.margin = margin(m, m, m, m),
    axis.title.x = element_text(size = 25, margin = margin(m,m,m,m)),
    axis.title.y = element_text(size = 25, margin = margin(m,m,m,m), hjust = 1, angle = 0)
  )

# measurement time histogram
m <- 20
ggplot(data = tmp) +
  geom_histogram(aes(x = daytime)) +
  xlim(0, 24) + 
  xlab('Time of Day') +
  ylab('Count') +
  ggtitle('Measurement Time Histogram') +
  theme(
    title = element_text(size = 25, margin = margin(m,m,m,m)), 
    plot.margin = margin(m, m, m, m),
    axis.title.x = element_text(size = 25, margin = margin(m,m,m,m)),
    axis.title.y = element_text(size = 25, margin = margin(m,m,m,m), hjust = 1, angle = 0)
  )


# QUICK LOOK AT CURRENT POLICIES 
#
#

# analyze current policies, id per id
na_to_0 <- function(x){
  v <- x
  v[is.na(x)] <- 0
  return(v)
}
id = 3
par(mfrow = c(2,2))
boxplot(na_to_0(tmp$value_regular_insulin_dose_[tmp$id==id]) ~ tmp$daytime_disc[tmp$id==id], main = 'Regular insulin, daily distribution')
boxplot(na_to_0(tmp$value_nph_insulin_dose_[tmp$id==id]) ~ tmp$daytime_disc[tmp$id==id], ylim = c(0,20), main = 'NPH insulin, daily distribution')
boxplot(na_to_0(tmp$value_ultralente_insulin_dose_[tmp$id==id]) ~ tmp$daytime_disc[tmp$id==id], main = 'UltraLente insulin, daily distribution')
par(mfrow = c(1,1))
