library(tidyverse)

test_coverage_old <- function(df, verbose=FALSE){
  
  results <- df %>%
    mutate(I = truth <= value) %>%
    group_by(quantile) %>%
    summarize(n1 = sum(I), 
              n0 = sum(1 - I), 
              p = mean(I), 
              
              L_uc_0 = (1 - unique(quantile))^n0 * unique(quantile)^n1,
              L_uc_1 = (1 - p)^n0 * p^n1,
              LR_uc = -2*log(L_uc_0/L_uc_1),
              p_uc = 1-pchisq(LR_uc, 1)
    )
  
  if(!verbose){
    results <- results %>%
      select(c(quantile, p_uc))
  }
  
  return(results)
}

# plot_coverage <- function(df){
#   df_temp <- df %>%
#     mutate(I_l = truth < value,
#            I_u = truth <= value) %>%
#     group_by(quantile) %>%
#     summarize(l1 = sum(I_l), 
#               u1 = sum(I_u),
#               l = mean(I_l),
#               u = mean(I_u),
#               p.value_l = binom.test(x=l1, n=n(), p=unique(quantile), alternative="greater")$p.value,
#               p.value_u = binom.test(x=u1, n=n(), p=unique(quantile), alternative="less")$p.value)
#   
#   ggplot(df_temp) +
#     geom_errorbar(aes(x=quantile, ymin=l, ymax=u), width=0.05,
#                   data=df_temp, colour="black") +
#     # geom_text(aes(label = round(p.value_l, 2), x=quantile, y = l), vjust = 1.5) +
#     # geom_text(aes(label = round(p.value_u, 2), x=quantile, y = u), vjust = -.5) +
#     geom_text(data=subset(df_temp, l > quantile), aes(label = signif(p.value_l, digits=2), x=quantile, y = u), vjust = -.5) +
#     geom_text(data=subset(df_temp, u < quantile), aes(label = signif(p.value_u, digits=2), x=quantile, y = l), vjust = 1.5) +
#     geom_segment(aes(x=0,xend=1,y=0,yend=1), linetype="dashed", colour="grey70")+
#     scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
#                        labels = function(x) ifelse(x == 0, "0", x)) +
#     scale_y_continuous(labels = function(y) ifelse(y == 0, "0", y)) +
#     xlab('Quantile') +
#     ylab('')
# }


test_coverage <- function(df, correction="none", verbose=FALSE){
  results <- df %>%
    mutate(I_l = truth < value,
           I_u = truth <= value) %>%
    group_by(quantile) %>%
    summarize(l1 = sum(I_l), 
              u1 = sum(I_u),
              l = mean(I_l),
              u = mean(I_u),
              p.value_l = binom.test(x=l1, n=n(), p=unique(quantile), alternative="greater")$p.value,
              p.value_u = binom.test(x=u1, n=n(), p=unique(quantile), alternative="less")$p.value)
  
  if (correction == 'bonferroni'){
    n_tests = 2*length(unique(results$quantile))
    results$p.value_l = pmin(results$p.value_l * n_tests, 1)
    results$p.value_u = pmin(results$p.value_u * n_tests, 1)
  }
  
  if (correction == 'bonferroni_ind'){
    results$p.value_l = pmin(results$p.value_l * 2, 1)
    results$p.value_u = pmin(results$p.value_u * 2, 1)
  }
  
  if(!verbose){
    results <- results %>%
      select(c(quantile, p.value_l, p.value_u))
  }
  
  return(results)
}


plot_coverage <- function(df, correction="none"){
  df_temp <- test_coverage(df, correction, verbose=TRUE)
  
  ggplot(df_temp) +
    geom_errorbar(aes(x=quantile, ymin=l, ymax=u), width=0.05,
                  data=df_temp, colour="black") +
    geom_text(data=subset(df_temp, l > quantile), aes(label = signif(p.value_l, digits=2), x=quantile, y = u), vjust = -.5) +
    geom_text(data=subset(df_temp, u < quantile), aes(label = signif(p.value_u, digits=2), x=quantile, y = l), vjust = 1.5) +
    geom_segment(aes(x=0,xend=1,y=0,yend=1), linetype="dashed", colour="grey70")+
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = function(x) ifelse(x == 0, "0", x)) +
    scale_y_continuous(labels = function(y) ifelse(y == 0, "0", y)) +
    xlab('Quantile') +
    ylab('')
}

plot_coverage(df)
plot_coverage(df, "bonferroni")
plot_coverage(df, "bonferroni_ind")

test_coverage(df)
test_coverage(df, "bonferroni")


### EXAMPLE

alphas <- round(seq(0.1, 0.9, 0.1), 3)

# Negative binomial
size=5
mu=10
sample <- rnbinom(300, size=size, mu=mu)
F <- qnbinom(p=alphas, size=size, mu=mu)
F <- data.frame(quantile=alphas, value=F)
df <- bind_rows(replicate(length(sample), F, simplify = FALSE))
df$truth <- rep(sample, each=length(alphas))
df$index <- rep(1:length(sample), each=length(alphas))

df$truth <- df$truth -1
res <- test_coverage(df)

test_coverage_old(df)

plot_coverage(df)
plot_coverage(df, "bonferroni_ind")


ps = c(res$p.value_l, res$p.value_u)
ps_adj = p.adjust(ps, 'holm')

ps = pmin(res$p.value_l, res$p.value_u)
ps_adj = p.adjust(ps, 'holm')


# Normal distribution
sample <- rnorm(3000)
F <- qnorm(p=alphas)
hist(sample)

F <- data.frame(quantile=alphas, value=F)
df <- bind_rows(replicate(length(sample), F, simplify = FALSE))
df$truth <- rep(sample, each=length(alphas))
df$index <- rep(1:length(sample), each=length(alphas))

res <- test_coverage(df)
plot_coverage(df)


sample_norm <- function(n=100, alphas=round(seq(0.1, 0.9, 0.1), 3), bias=0){
  sample <- rnorm(n)
  F <- qnorm(p=alphas)
  F <- data.frame(quantile=alphas, value=F)
  df <- bind_rows(replicate(length(sample), F, simplify = FALSE))
  df$truth <- rep(sample, each=length(alphas))
  df$truth <- df$truth + bias
  df$index <- rep(1:length(sample), each=length(alphas))
  return(df)
}

sample_negbin <- function(n=100, size=5, mu=10, alphas=round(seq(0.1, 0.9, 0.1), 3)){
  sample <- rnbinom(n, size=size, mu=mu)
  F <- qnbinom(p=alphas, size=size, mu=mu)
  F <- data.frame(quantile=alphas, value=F)
  df <- bind_rows(replicate(length(sample), F, simplify = FALSE))
  df$truth <- rep(sample, each=length(alphas))
  df$index <- rep(1:length(sample), each=length(alphas))
  return(df)
}


res <- test_coverage(sample_norm())
plot_coverage(sample_norm())

plot_coverage(sample_norm(bias=-1))


### Simulation study

# standard normal

results_df = data.frame()
for (i in 1:5000){
  res <- test_coverage(sample_norm(500, bias=0), "bonferroni_ind")
  res$index <- i
  results_df <- bind_rows(results_df, res)
}

results_df %>%
  group_by(quantile) %>%
  summarize(fps = mean(p.value_l < 0.05 | p.value_u < 0.05))

results_df %>%
  group_by(quantile) %>%
  summarize(fns = mean(p.value_l > 0.05 & p.value_u > 0.05))



hist(c(results_df$p.value_l, results_df$p.value_u))

false_positives <- results_df %>%
  group_by(quantile) %>%
  summarize(std_norm = mean(p_uc < 0.05))


# negative binomial 
test_coverage(sample_negbin(500, mu=10))
plot_coverage(sample_negbin(500, mu=10))
plot_coverage(sample_negbin(500, mu=100))

hist(rnbinom(500, size=5, mu=10))


results_nb10 = data.frame()
for (i in 1:5000){
  res <- test_coverage(sample_negbin(500, mu=10))
  res$index <- i
  results_nb10 <- bind_rows(results_nb10, res)
}

results_nb10 %>%
  group_by(quantile) %>%
  summarize(fps = mean(p.value_l < 0.05 | p.value_u < 0.05))



results_nb50 = data.frame()
for (i in 1:5000){
  res <- test_coverage(sample_negbin(500, mu=50))
  res$index <- i
  results_nb50 <- bind_rows(results_nb50, res)
}

results_nb50 %>%
  group_by(quantile) %>%
  summarize(fps = mean(p.value_l < 0.05 | p.value_u < 0.05))






hist(c(results_nb$p.value_l, results_nb$p.value_u))
hist(results_nb$p.value_l)
hist(results_nb$p.value_u)

results_nb50 %>%
  group_by(quantile) %>%
  summarize(fps = mean(mean(p.value_l < 0.05), mean(p.value_u < 0.05)))

results_nb %>%
  group_by(quantile) %>%
  summarize(fps = mean(mean(p.value_l < 0.05), mean(p.value_u < 0.05)))

results_df %>%
  group_by(index) %>%
  summarize(fp = any(p.value_l < 0.05 | p.value_u < 0.05)) %>%
  summarize(mean(fp))

results_df %>%
  group_by(quantile) %>%
  summarize(fps = mean(mean(p.value_l < 0.05), mean(p.value_u < 0.05)))

hist(c(results_df$p.value_l, results_df$p.value_u))


false_positives10 <- results_df %>%
  group_by(quantile) %>%
  summarize(false_positive = mean(p_uc < 0.05))

false_positives20 <- results_df %>%
  group_by(quantile) %>%
  summarize(false_positive = mean(p_uc < 0.05))

false_positives100 <- results_df %>%
  group_by(quantile) %>%
  summarize(false_positive = mean(p_uc < 0.05))


false_positives10_1000 <- false_positives10

fp <- left_join(left_join(false_positives10, false_positives20, by='quantile'), false_positives100, by='quantile')

fp <- fp %>%
  rename(mu10=false_positive.x, mu20=false_positive.y, mu100=false_positive)

write.csv(fp, "data/examples/false_positives.csv", row.names=FALSE)

fp2 <- left_join(false_positives, fp)
write.csv(fp2, "data/examples/false_positives_all.csv", row.names=FALSE)

print(xtable(fp2, digits=c(0, 2, 3, 3, 3, 3)), include.rownames=FALSE)



### NEW SIMULATION STUDY

r = 10000
n = 500
ns = c(50, 100, 250, 500, 1000)

fp_all = data.frame()

for (n in ns){
  results_df = data.frame()
  for (i in 1:r){
    res <- test_coverage(sample_norm(n, bias=0))
    res$index <- i
    results_df <- bind_rows(results_df, res)
  }
  
  fp_temp <- results_df %>%
    group_by(quantile) %>%
    summarize(fp = mean(p.value_l < 0.05 | p.value_u < 0.05))
  
  fp_temp$n <- n
  
  fp_all <- bind_rows(fp_all, fp_temp)
}

a <- fp_all %>%
  filter(quantile==0.1)

ggplot(a, aes(x=n, y=fp)) +
  geom_line()


ggplot(fp_all, aes(x=n, y=fp)) +
  facet_wrap('quantile', scales="free_y") +
  geom_line() +
  labs(title="Standard normal")


fp_all2 = data.frame()

test_coverage(sample_negbin(500, 10))

for (n in ns){
  results_df = data.frame()
  for (i in 1:r){
    res <- test_coverage(sample_negbin(n, mu=10))
    res$index <- i
    results_df <- bind_rows(results_df, res)
  }
  
  fp_temp <- results_df %>%
    group_by(quantile) %>%
    summarize(fp = mean(p.value_l < 0.05 | p.value_u < 0.05))
  
  fp_temp$n <- n
  
  fp_all2 <- bind_rows(fp_all2, fp_temp)
}


a <- fp_all2 %>%
  filter(quantile==0.5)

ggplot(a, aes(x=n, y=fp)) +
  geom_line()

ggplot(fp_all2, aes(x=n, y=fp)) +
  facet_wrap('quantile', scales="free_y") +
  geom_line() +
  labs(title="Negative binomial, new test")


fp_all3 = data.frame()

for (n in ns){
  results_df = data.frame()
  for (i in 1:r){
    res <- test_coverage_old(sample_negbin(n, mu=10))
    res$index <- i
    results_df <- bind_rows(results_df, res)
  }
  
  fp_temp <- results_df %>%
    group_by(quantile) %>%
    summarize(fp = mean(p_uc < 0.05))
  
  fp_temp$n <- n
  
  fp_all3 <- bind_rows(fp_all3, fp_temp)
}

compute_ps <- function(sample_sizes = c(50, 100, 250, 500, 1000), n_rep = 10000){
  results_all = data.frame()
  for (n in sample_sizes){
    print(n)
    results_df = data.frame()
    for (i in 1:n_rep){
      
      # standard normal
      sample_df <- sample_norm(n, bias=0)
      
      ps <- test_coverage_old(sample_df)
      ps$index <- i
      ps$test <- 'old'
      ps$distr <- 'standard_normal'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df)
      ps$index <- i
      ps$test <- 'new'
      ps$distr <- 'standard_normal'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df, correction="bonferroni_ind")
      ps$index <- i
      ps$test <- 'new_bonf_ind'
      ps$distr <- 'standard_normal'
      results_df <- bind_rows(results_df, ps)
      
      
      # negative binomial, 10
      sample_df <- sample_negbin(n, mu=10)
      
      ps <- test_coverage_old(sample_df)
      ps$index <- i
      ps$test <- 'old'
      ps$distr <- 'negbin_10'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df)
      ps$index <- i
      ps$test <- 'new'
      ps$distr <- 'negbin_10'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df, correction="bonferroni_ind")
      ps$index <- i
      ps$test <- 'new_bonf_ind'
      ps$distr <- 'negbin_10'
      results_df <- bind_rows(results_df, ps)
      
      
      # negative binomial, 20
      sample_df <- sample_negbin(n, mu=20)
      
      ps <- test_coverage_old(sample_df)
      ps$index <- i
      ps$test <- 'old'
      ps$distr <- 'negbin_20'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df)
      ps$index <- i
      ps$test <- 'new'
      ps$distr <- 'negbin_20'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df, correction="bonferroni_ind")
      ps$index <- i
      ps$test <- 'new_bonf_ind'
      ps$distr <- 'negbin_20'
      results_df <- bind_rows(results_df, ps)
      
      
      # negative binomial, 50
      sample_df <- sample_negbin(n, mu=50)
      
      ps <- test_coverage_old(sample_df)
      ps$index <- i
      ps$test <- 'old'
      ps$distr <- 'negbin_50'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df)
      ps$index <- i
      ps$test <- 'new'
      ps$distr <- 'negbin_50'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df, correction="bonferroni_ind")
      ps$index <- i
      ps$test <- 'new_bonf_ind'
      ps$distr <- 'negbin_50'
      results_df <- bind_rows(results_df, ps)
      
      
      # negative binomial, 100
      sample_df <- sample_negbin(n, mu=100)
      
      ps <- test_coverage_old(sample_df)
      ps$index <- i
      ps$test <- 'old'
      ps$distr <- 'negbin_100'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df)
      ps$index <- i
      ps$test <- 'new'
      ps$distr <- 'negbin_100'
      results_df <- bind_rows(results_df, ps)
      
      ps <- test_coverage(sample_df, correction="bonferroni_ind")
      ps$index <- i
      ps$test <- 'new_bonf_ind'
      ps$distr <- 'negbin_100'
      results_df <- bind_rows(results_df, ps)
      
    }
    results_df$n <- n 
    results_all <- bind_rows(results_all, results_df)
  }
  return(results_all)
}

# b <- a
a <- compute_ps(sample_sizes = c(50, 100, 250, 500, 1000), n_rep = 10000)


write.csv(a, "sim_results_all.csv", row.names=FALSE)

a$p <- a$p_uc

a <- a %>%
  mutate(p = pmin(p_uc, p.value_l, p.value_u, na.rm=TRUE))

ggplot(subset(a, test=="new" & distr=="standard_normal" & n==1000), aes(x=p)) +
  facet_wrap("quantile") +
  geom_histogram(binwidth = 0.05) 

ggplot(subset(a, test=="old" & distr=="standard_normal" & n==500), aes(x=p)) +
  facet_wrap("quantile") +
  geom_histogram(binwidth = 0.1) 

ggplot(subset(a, test=="old" & distr=="negbin_10" & n==500), aes(x=p)) +
  facet_wrap("quantile") +
  geom_histogram(binwidth = 0.05) 

ggplot(subset(a, test=="new" & distr=="negbin_10" & n==1000), aes(x=p)) +
  facet_wrap("quantile") +
  geom_histogram(binwidth = 0.05) 

ggplot(subset(a, test=="new" & distr=="negbin_20" & n==1000), aes(x=p)) +
  facet_wrap("quantile") +
  geom_histogram(binwidth = 0.05) 

ggplot(subset(a, test=="new" & distr=="negbin_50" & n==1000), aes(x=p)) +
  facet_wrap("quantile") +
  geom_histogram(binwidth = 0.05) 

ggplot(subset(a, test=="new" & distr=="negbin_100" & n==250), aes(x=p)) +
  facet_wrap("quantile") +
  geom_histogram(binwidth = 0.1) 


ggplot(fp_all3, aes(x=n, y=fp)) +
  facet_wrap('quantile', scales="free_y") +
  geom_line() +
  labs(title="Negative binomial, old test")

test_coverage_old(sample_negbin(50, mu=10))


ggplot(results_df, aes(x=p_uc)) +
  facet_wrap("quantile") +
  geom_histogram(binwidth = 0.05) 

ggplot(results_df, aes(x=p_uc)) +
  facet_wrap("quantile") +
  geom_histogram(binwidth = 0.05, aes(y=10*..count../sum(..count..)))

# find critical value
C <- 0
while(pbinom(C, 500, 0.1) < 0.05) C <- C+1
C <- C - 1

get_cr <- function(sample_size, p, alpha){
  C <- 0
  while(pbinom(C, sample_size, p) < alpha) C <- C+1
  C - 1
}

get_cr <- function(sample_size, p, alpha, alternative='less'){
  if(alternative == 'less'){
    C <- 0
    while(pbinom(C, sample_size, p) < alpha) C <- C+1
    C - 1
  } else if(alternative == 'greater'){
    C <- sample_size
    while (1 - pbinom(C-1, sample_size, p) < alpha) C <- C -1
    C + 1
  }

}

get_cr(100, 0.1, 0.05, 'less')
get_cr(100, 0.1, 0.05, 'greater')

1- pbinom(15, 100, 0.1)


C <- 500
while (1 - pbinom(C-1, 500, 0.1) < 0.05) C <- C -1
C <- C + 1


get_cr(50, 0.1, 0.05)
get_cr(100, 0.1, 0.05)
get_cr(500, 0.1, 0.05)
get_cr(1000, 0.1, 0.05)
pbinom(84, 1000, 0.1)

df <- data.frame(p=seq(0.1, 0.14, 0.005))
df$y50 <- sapply(df$p, function(x){pbinom(get_cr(50, 0.1, 0.05), 50, x)})
df$y100 <- sapply(df$p, function(x){pbinom(get_cr(100, 0.1, 0.05), 100, x)})
df$y500 <- sapply(df$p, function(x){pbinom(get_cr(500, 0.1, 0.05), 500, x)})
df$y1000 <- sapply(df$p, function(x){pbinom(get_cr(1000, 0.1, 0.05), 1000, x)})

df <- pivot_longer(df, !p)

ggplot(df, aes(x=p, y=value)) + 
  geom_line(aes(color=name)) +
  labs(y="P_x(T <= Crit. value)")

pbinom(C, 500, 0.1)
pbinom(C, 500, 0.11)



df <- data.frame(p=seq(0.05, 0.14, 0.005))
df$y50 <- sapply(df$p, function(x){1 - pbinom(get_cr(50, 0.1, 0.05, "greater")-1, 50, x)})
df$y100 <- sapply(df$p, function(x){1 - pbinom(get_cr(100, 0.1, 0.05, "greater")-1, 100, x)})
df$y500 <- sapply(df$p, function(x){1 - pbinom(get_cr(500, 0.1, 0.05, "greater")-1, 500, x)})
df$y1000 <- sapply(df$p, function(x){1 - pbinom(get_cr(1000, 0.1, 0.05, "greater")-1, 1000, x)})

df <- pivot_longer(df, !p)

ggplot(df, aes(x=p, y=value)) + 
  geom_line(aes(color=name)) +
  labs(y="P_x(T >= Crit. value)")


### POWER FUNCTIONS

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

df2 <- data.frame(q=seq(0, 1, 0.005))
df2$y50 <- sapply(df2$q, function(x){1 - pbinom(get_cr(50, 0.1, 0.025, "greater")-1, 50, x)})
df2$y100 <- sapply(df2$q, function(x){1 - pbinom(get_cr(100, 0.1, 0.025, "greater")-1, 100, x)})
df2$y500 <- sapply(df2$q, function(x){1 - pbinom(get_cr(500, 0.1, 0.025, "greater")-1, 500, x)})
df2$y1000 <- sapply(df2$q, function(x){1 - pbinom(get_cr(1000, 0.1, 0.025, "greater")-1, 1000, x)})

df2 <- pivot_longer(df2, !q, names_to='n')
df2$n <- factor(str_sub(df2$n, start=2), levels=c("50", "100", "500", "1000"))

ggplot(df2, aes(x=q, y=value)) + 
  geom_line(aes(color=n), size=1.05) +
  labs(y="g(p)") +
  geom_hline(yintercept=0.05, linetype='dashed', size=1.05)

# combined

q = 0.1

df <- data.frame(p=seq(0, 1, 0.005))
df$y50 <- sapply(df$p, function(x){pbinom(get_cr(50, q, 0.025), 50, x)})
df$y100 <- sapply(df$p, function(x){pbinom(get_cr(100, q, 0.025), 100, x)})
df$y500 <- sapply(df$p, function(x){pbinom(get_cr(500, q, 0.025), 500, x)})
df$y1000 <- sapply(df$p, function(x){pbinom(get_cr(1000, q, 0.025), 1000, x)})

df$y50_2 <- sapply(df$p, function(x){1 - pbinom(get_cr(50, q, 0.025, "greater")-1, 50, x)})
df$y100_2 <- sapply(df$p, function(x){1 - pbinom(get_cr(100, q, 0.025, "greater")-1, 100, x)})
df$y500_2 <- sapply(df$p, function(x){1 - pbinom(get_cr(500, q, 0.025, "greater")-1, 500, x)})
df$y1000_2 <- sapply(df$p, function(x){1 - pbinom(get_cr(1000, q, 0.025, "greater")-1, 1000, x)})

df$y50 <- df$y50 + df$y50_2
df$y100 <- df$y100 + df$y100_2
df$y500 <- df$y500 + df$y500_2
df$y1000 <- df$y1000 + df$y1000_2

df <- df %>%
  select(-c(y50_2, y100_2, y500_2, y1000_2))

df <- pivot_longer(df, !p, names_to='n')

df$n <- factor(str_sub(df$n, start=2), levels=c("50", "100", "500", "1000"))

ggplot(df, aes(x=p, y=value)) + 
  geom_line(aes(color=n), size=1.05) +
  labs(y="g(p)") +
  geom_hline(yintercept=0.05, linetype='dashed', size=1.05)


### Plot densities

size=5
mu=10
sample <- rnbinom(300, size=size, mu=mu)
hist(sample)
max(sample)
dnbinom(1:100, size=size, mu=mu)

x=1:150
df <- data.frame(x = x,
                 mu10 = dnbinom(x, size=5, mu=10),
                 mu20 = dnbinom(x, size=5, mu=20),
                 mu100 = dnbinom(x, size=5, mu=100))

df <- pivot_longer(df, cols=c(mu10, mu20, mu100), names_to='mu')

ggplot(df, aes(x=x, y=value, fill=mu)) + 
  geom_bar(stat='identity', position = 'dodge') +
  scale_fill_discrete(limits=c("mu10", "mu20", "mu100"),  labels = c("10", "20", "100"), name="mu") +
  labs(y='Probability') +
  theme_gray(base_size=20)

ggsave('plots/examples/neg_bin_ex.png', width=36, height=14, dpi=500, unit='cm', device='png')
