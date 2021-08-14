library(tidyverse)
df <- read.csv("sim_results_all.csv")

df <- df %>%
  mutate(p = pmin(p_uc, p.value_l, p.value_u, na.rm=TRUE))

df <- df %>%
  filter(n %in% c(100, 500, 1000))

plot_pvalue_histogram <- function(df, test="new", distr="standard_normal"){
  df %>%
    filter(test == !!test & distr == !!distr) %>%
    ggplot(aes(x=p)) +
      facet_grid(rows=vars(quantile), cols=vars(n)) +
      geom_histogram(breaks = seq(0, 1, 0.1)) +
      labs(title=paste0("Test: ", test, ", distribution: ", distr))
}

plot_pvalue_histogram(df, test="old", distr="standard_normal")
plot_pvalue_histogram(df, test="old", distr="negbin_10")
plot_pvalue_histogram(df, test="old", distr="negbin_100")

plot_pvalue_histogram(df, test="new", distr="standard_normal")
plot_pvalue_histogram(df, test="new", distr="negbin_10")
plot_pvalue_histogram(df, test="new", distr="negbin_100")

plot_pvalue_histogram(df, test="new_bonf_ind", distr="standard_normal")
plot_pvalue_histogram(df, test="new_bonf_ind", distr="negbin_10")
plot_pvalue_histogram(df, test="new_bonf_ind", distr="negbin_50")
plot_pvalue_histogram(df, test="new_bonf_ind", distr="negbin_100")


false_positive_rates <- df %>%
  group_by(test, distr, n, quantile) %>%
  summarize(fp = mean(p < 0.05))

false_positive_rates %>%
  filter(test == "new_bonf_ind") %>%
  ggplot(aes(x = n, y = fp, color = test)) +
    geom_line() +
    facet_grid(rows=vars(distr), cols=vars(quantile)) +
    geom_hline(yintercept = 0.05, linetype="dashed")


### POWER FUNCTIONS
get_cr <- function(sample_size, p, alpha, alternative='less'){
  if(alternative == 'less'){
    C <- 0
    while(pbinom(C + 1, sample_size, p) < alpha) C <- C+1
    # C - 1
  } else if(alternative == 'greater'){
    C <- sample_size
    while (1 - pbinom(C-1, sample_size, p) < alpha) C <- C -1
    # C + 1
  }
  return(C)
}

df <- data.frame(p=seq(0, 0.25, 0.005))
df$y50 <- sapply(df$p, function(x){pbinom(get_cr(50, 0.1, 0.05), 50, x)})
df$y100 <- sapply(df$p, function(x){pbinom(get_cr(100, 0.1, 0.05), 100, x)})
df$y500 <- sapply(df$p, function(x){pbinom(get_cr(500, 0.1, 0.05), 500, x)})
df$y1000 <- sapply(df$p, function(x){pbinom(get_cr(1000, 0.1, 0.05), 1000, x)})

df <- pivot_longer(df, !p, names_to='n')
df$n <- factor(str_sub(df$n, start=2), levels=c("50", "100", "500", "1000"))

ggplot(df, aes(x=p, y=value)) + 
  geom_line(aes(color=n), size=1.05) +
  labs(y="g(p)") +
  geom_hline(yintercept=0.05, linetype='dashed', size=1.05)

df2 <- data.frame(q=seq(0, 0.4, 0.005))
df2$y50 <- sapply(df2$q, function(x){1 - pbinom(get_cr(50, 0.1, 0.025, "greater"), 50, x)})
df2$y100 <- sapply(df2$q, function(x){1 - pbinom(get_cr(100, 0.1, 0.025, "greater"), 100, x)})
df2$y500 <- sapply(df2$q, function(x){1 - pbinom(get_cr(500, 0.1, 0.025, "greater"), 500, x)})
df2$y1000 <- sapply(df2$q, function(x){1 - pbinom(get_cr(1000, 0.1, 0.025, "greater"), 1000, x)})

df2 <- pivot_longer(df2, !q, names_to='n')
df2$n <- factor(str_sub(df2$n, start=2), levels=c("50", "100", "500", "1000"))

ggplot(df2, aes(x=q, y=value)) + 
  geom_line(aes(color=n), size=1.05) +
  labs(y="g(q)") +
  geom_hline(yintercept=0.05, linetype='dashed', size=1.05)


#####
# combined

q = 0.1
alpha = 0.05
xmax = 0.4

df1 <- data.frame(p=seq(0, xmax, 0.005))
df1$y50 <- sapply(df1$p, function(x){pbinom(get_cr(50, q, alpha), 50, x)})
df1$y100 <- sapply(df1$p, function(x){pbinom(get_cr(100, q, alpha), 100, x)})
df1$y500 <- sapply(df1$p, function(x){pbinom(get_cr(500, q, alpha), 500, x)})
df1$y1000 <- sapply(df1$p, function(x){pbinom(get_cr(1000, q, alpha), 1000, x)})

df1 <- pivot_longer(df1, !p, names_to='n')
df1$n <- factor(str_sub(df1$n, start=2), levels=c("50", "100", "500", "1000"))
df1$test <- "left_tailed"

df2 <- data.frame(p=seq(0, xmax, 0.005))
df2$y50 <- sapply(df2$p, function(x){1 - pbinom(get_cr(50, q, alpha, "greater"), 50, x)})
df2$y100 <- sapply(df2$p, function(x){1 - pbinom(get_cr(100, q, alpha, "greater"), 100, x)})
df2$y500 <- sapply(df2$p, function(x){1 - pbinom(get_cr(500, q, alpha, "greater"), 500, x)})
df2$y1000 <- sapply(df2$p, function(x){1 - pbinom(get_cr(1000, q, alpha, "greater"), 1000, x)})

df2 <- pivot_longer(df2, !p, names_to='n')
df2$n <- factor(str_sub(df2$n, start=2), levels=c("50", "100", "500", "1000"))
df2$test <- "right_tailed"

df <- bind_rows(df1, df2)

ggplot(df, aes(x=p, y=value)) +
  facet_wrap("n") +
  geom_line(aes(color=test), size=1.05) +
  labs(y="g(p)") +
  geom_hline(yintercept=0.05, linetype='dashed', size=1.05)


#####
names(df2)[1] <- "q"

df1 <- df1 %>%
  filter(n == 1000) %>%
  select(-c(n, test))

df2 <- df2 %>%
  filter(n == 1000) %>%
  select(-c(n, test))

df <- merge(df1, df2, by = NULL)

df <-df %>%
  mutate(value = pmax(value.x, value.y))

df <-df %>%
  mutate(value = value.x + value.y - value.x*value.y)

ggplot(subset(df), aes(p, q, z = value)) + 
  geom_contour_filled() +
  coord_fixed()



ggplot(subset(df, q <= p), aes(p, q, z = value)) + 
  geom_contour_filled() +
  coord_fixed()

ggplot(subset(df, q <= p), aes(p, q, z = value)) + 
  geom_raster(aes(fill = value)) +
  geom_contour(colour = "white")

ggplot(subset(df, p <= 0.115 & p >= 0.05 & q >= 0.07 & q <= p), aes(p, q, z = value)) + 
  geom_contour_filled()


df <- df %>%
  arrange(p, q)

persp(df$p, df$q, df$value)

library(rgl)
persp3d(df$p, df$q, df$value)


library(plotly)

df2 <- subset(df, q <= p)
plot_ly() %>% add_trace(data = df2, x = df2$p, y = df2$q, z = df2$value, type="mesh3d")

library(akima)


s <- interp(df$p, df$q, df$value, duplicate = "mean")
plot_ly() %>% add_surface(x = s$x, y = s$y, z = t(s$z), type="mesh3d")


# a <- df %>%
#   filter(q <= p) %>%
#   mutate(v = value.x*value.y)

# all ns
names(df2)[1] <- "q"

df1 <- df1 %>%
  select(-test)

df2 <- df2 %>%
  select(-test)

df <- merge(df1, df2, by = NULL)

df <- df %>%
  filter(n.y == n.x) %>%
  mutate(n = n.x) %>%
  select(-c(n.x, n.y))

# df <-df %>%
#   mutate(value = pmax(value.x, value.y))

df <-df %>%
  mutate(value = value.x + value.y - value.x*value.y)

library(latex2exp)

ggplot(subset(df, q <= p), aes(p, q, z = value)) + 
  facet_wrap("n") +
  geom_contour_filled() +
  coord_fixed() +
  geom_hline(yintercept = 0.1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.1, color = "red", linetype = "dashed") +
  labs(fill = "Rejection probability") +
  xlab(unname(TeX("p $ = P(Y \\leq q_{0.1}) $"))) +
  ylab(unname(TeX("q $ = P(Y < q_{0.1}) $")))





qnbinom(p=0.1, size=5, mu=10)
pnbinom(q=4, size=5, mu=10)
pnbinom(q=3, size=5, mu=10)


?pnbinom
