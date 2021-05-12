library(dplyr)
library(ggplot2)

test_coverage <- function(df, verbose=FALSE){
  
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


df_temp <- df %>%
  select(quantile | value | truth)

df_temp$l <- df_temp$truth < floor(df_temp$value)
df_temp$u <- df_temp$truth <= floor(df_temp$value)

df_temp <- df_temp %>%
  group_by(quantile) %>%
  summarize(l = mean(l), u=mean(u))



results <- df %>%
  mutate(I_l = truth < floor(value),
         I_u = truth <= floor(value)) %>%
  group_by(quantile) %>%
  summarize(l1 = sum(I_l), 
            l0 = sum(1 - I_l), 
            c_l = mean(I_l), 
            
            u1 = sum(I_u), 
            u0 = sum(1 - I_u), 
            c_u = mean(I_u)
  )



test_coverage <- function(df, correction="none"){
  results <- df %>%
    mutate(I_l = truth < (value),
           I_u = truth <= (value)) %>%
    group_by(quantile) %>%
    summarize(l1 = sum(I_l), 
              u1 = sum(I_u),
              p.value_l = binom.test(x=l1, n=n(), p=unique(quantile), alternative="greater")$p.value,
              p.value_u = binom.test(x=u1, n=n(), p=unique(quantile), alternative="less")$p.value)
  
  if (correction == 'bonferroni'){
    results$p.value_l = pmin(results$p.value_l * 2*length(unique(results$quantile)), 1)
    results$p.value_u = pmin(results$p.value_u * 2*length(unique(results$quantile)), 1)
  }
  return(results)
}

test_coverage(df)

results <- df %>%
  mutate(I_l = truth < floor(value),
         I_u = truth <= floor(value)) %>%
  group_by(quantile) %>%
  summarize(l1 = sum(I_l), 
            u1 = sum(I_u),
            p.value_l = binom.test(x=l1, n=n(), p=unique(quantile), alternative="greater")$p.value,
            p.value_u = binom.test(x=u1, n=n(), p=unique(quantile), alternative="less")$p.value)

mutate(CI = list(enframe(binom.test(Defaults, Count, PD/100, alternative =  "two.sided", conf.level = 0.90)$conf.int))) %>%

binom.test(x = 80, n = 160, p = 3/5, alternative = 'less')$p.value

plot_coverage <- function(df){
  df_temp <- df %>%
    mutate(I_l = truth < value,
           I_u = truth <= value) %>%
    group_by(quantile) %>%
    summarize(l1 = sum(I_l), 
              u1 = sum(I_u),
              l = mean(I_l),
              u = mean(I_u),
              p.value_l = binom.test(x=l1, n=n(), p=unique(quantile), alternative="greater")$p.value,
              p.value_u = binom.test(x=u1, n=n(), p=unique(quantile), alternative="less")$p.value)
  
  ggplot(df_temp) +
    geom_errorbar(aes(x=quantile, ymin=l, ymax=u), width=0.05,
                  data=df_temp, colour="black") +
    # geom_text(aes(label = round(p.value_l, 2), x=quantile, y = l), vjust = 1.5) +
    # geom_text(aes(label = round(p.value_u, 2), x=quantile, y = u), vjust = -.5) +
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


### EXAMPLE

alphas <- round(seq(0.1, 0.9, 0.1), 3)

# Negative binomial
size=5
mu=10
sample <- rnbinom(300, size=size, mu=mu)
F <- qnbinom(p=alphas, size=size, mu=mu)
hist(sample)

F <- data.frame(quantile=alphas, value=F)
df <- bind_rows(replicate(length(sample), F, simplify = FALSE))
df$truth <- rep(sample, each=length(alphas))
df$index <- rep(1:length(sample), each=length(alphas))

df$truth <- df$truth -1
res <- test_coverage(df)
plot_coverage(df)

res <- test_coverage(df, verbose=FALSE)

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


sample_norm <- function(n=100, alphas=round(seq(0.1, 0.9, 0.1), 3)){
  sample <- rnorm(n)
  F <- qnorm(p=alphas)
  F <- data.frame(quantile=alphas, value=F)
  df <- bind_rows(replicate(length(sample), F, simplify = FALSE))
  df$truth <- rep(sample, each=length(alphas))
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

results_df = data.frame()
for (i in 1:10000){
  res <- test_coverage(sample_norm(5000))
  res$index <- i
  results_df <- bind_rows(results_df, res)
}
# mean(results_df$p_cc < 0.05)

results_df %>%
  group_by(quantile) %>%
  summarize(fps = mean(mean(p.value_l < 0.05), mean(p.value_u < 0.05)))

hist(c(results_df$p.value_l, results_df$p.value_u))

false_positives <- results_df %>%
  group_by(quantile) %>%
  summarize(std_norm = mean(p_uc < 0.05))


test_coverage(sample_negbin(500, mu=10))
plot_coverage(sample_negbin(500, mu=10))


results_df = data.frame()
for (i in 1:1000){
  res <- test_coverage(sample_negbin(500, mu=10), correction='bonferroni')
  res$index <- i
  results_df <- bind_rows(results_df, res)
}
# mean(results_df$p_cc < 0.05)


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
