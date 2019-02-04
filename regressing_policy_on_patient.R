ggplot(data = tmp_ %>% filter(value_regular_insulin_dose_ <90)) +
  geom_histogram(mapping = aes(x = value_regular_insulin_dose_))

boxplot(tmp_$value_regular_insulin_dose_, 
        tmp_$value_nph_insulin_dose_,
        tmp_$value_ultralente_insulin_dose_, ylim = c(0, 40))

m <- 20
ggplot(data = tmp_) +
  geom_histogram(aes(x = na_to_0(value_regular_insulin_dose_)), bins = 30) +
  xlim(0, 25) + 
  xlab('Regular Insulin Dose') +
  ylab('Count') +
  ggtitle('Regular Insulin Dose Histogram') +
  theme(
    title = element_text(size = 25, margin = margin(m,m,m,m)), 
    plot.margin = margin(m, m, m, m),
    axis.title.x = element_text(size = 25, margin = margin(m,m,m,m)),
    axis.title.y = element_text(size = 25, margin = margin(m,m,m,m), hjust = 1, angle = 0)
  )

x <- tmp_$id

y <- na_to_0(tmp_$value_regular_insulin_dose_)
# id_sel <- as.numeric(names(which(tapply(y, x, sum)>0)))
# y_sel <- y[x %in% id_sel]
# x_sel <- factor(x[x %in% id_sel])
lm0 <- lm(y~factor(x)-1)
summary(lm0)
hist(as.numeric(lm0$coefficients))

y <- na_to_0(tmp_$value_nph_insulin_dose_)
lm0 <- lm(y~factor(x)-1)
summary(lm0)
hist(as.numeric(lm0$coefficients), breaks = 20)

y <- na_to_0(tmp_$value_ultralente_insulin_dose_)
lm0 <- lm(y~factor(x)-1)
summary(lm0)
