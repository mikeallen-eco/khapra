# Sys.setenv(JAGS_HOME="C:\\Program Files\\JAGS\\JAGS-4.2.0")

library(tidyverse)
library(scales)
library(patchwork)

# load in bayesian model data
khapra_PCR <- readRDS("output/khapra_PCR_jags.rds")


# get samples for n. tech reps 1-6
det_df = map_df(1:12, ~tibble(ntr = .x, p_det = khapra_PCR$sims.list$p.PCRrep.low[, .x]))

det_df %>% 
  ggplot(aes(group = ntr, y = p_det)) + 
  geom_boxplot()


# p absnece
p_abs = function(S, Pprior){
  A = 1 - ( ((1-S)*Pprior) / ( ((1-S)*Pprior) + 1 - Pprior ) )
  return(A)
}

## test
abs_df  = det_df %>% mutate(p_abs = p_abs(p_det, 0.5))
abs_df %>% 
  ggplot(aes(group = ntr, y = p_abs)) + 
  geom_boxplot()


# run for i iterations for Pprior; rbeta(1e6, 0.25, 4)
temp = rbeta(1e6, 0.25, 4); plot(density(temp)); quantile(temp, c(0.025, 0.1, 0.5, 0.9, 0.975))
i = 100
abs_df = map_df(1:i, function(x){
  rvals = rbeta(i, 0.25, 4)
  det_df %>% 
    mutate(p_abs = p_abs(p_det, rvals[x])) %>% 
    mutate(i = x) %>% 
    select(i, ntr, p_abs)
  }  
)


## process results
## get summary stats
dfs1 = 
  abs_df %>% 
  group_by(ntr) %>% 
  summarise(q2.5 = quantile(p_abs, 0.025),
            q10 = quantile(p_abs, 0.1),
            med = median(p_abs),
            q90 = quantile(p_abs, 0.9),
            q97.5 = quantile(p_abs, 0.975)) %>% 
  ungroup()


dfs1 %>%
    ggplot() +
    geom_abline(slope = 0, intercept = 0.95, 
                color = "darkgray", size=1, linetype = 2) +
    geom_errorbar(aes(x = ntr, ymin = q2.5, ymax = q97.5),
                  width = 0, size = .5) +
    geom_errorbar(aes(x = ntr, ymin = q10, ymax = q90),
                  width = 0, size = 1) +
    geom_line(aes(x = ntr, y = med), size = 1) +
    geom_point(aes(x = ntr, y = med), 
               size = 4, color = "black") +
    theme_bw() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15)) +
    labs(y = "Probability of absence given no detections\n~ 1 beetle density", 
         x = "Number of qPCR technical replicates") +
    scale_x_continuous(breaks=pretty_breaks(6), limits = c(1, 6)) + 
    theme(panel.grid.minor.x = element_blank()) # 


ggsave(plot = last_plot(), filename = "plots/prob_abs.png",
       width = 6.5, height = 5, dpi = 300)

# Sensitivity Analysis Under Different Ppriors

## define Pprior & plot

### 0. Beta(0.25, 4)
pp0 = dbeta(seq(0,1,0.01), 0.25, 4) %>% 
  tibble() %>% mutate(x = seq(0,1,0.01)) %>% 
  ggplot(aes(x=x, y = .)) +
  geom_line() + 
  scale_x_continuous(limits = c(0,1), breaks = scales::pretty_breaks(5)) + 
  labs(x = "probability", y = "density", title = "Beta(0.25, 4)")
pp0

### 1. Uniform(0.01, 0.1)
pp1 = dunif(seq(0,1,0.01), 0.01, 0.1) %>% 
  tibble() %>% mutate(x = seq(0,1,0.01)) %>% 
  ggplot(aes(x=x, y=.)) +
  # geom_density(adjust = 2) + 
  geom_line() + 
  scale_x_continuous(limits = c(0,1), breaks = scales::pretty_breaks(5)) +
  labs(x = "probability", y = "density", title = "Uniform(0.01, 0.10)")
pp1

### 2. Beta(1, 5)
pp2 = dbeta(seq(0,1,0.01), 1, 5) %>% 
  tibble() %>% mutate(x = seq(0,1,0.01)) %>% 
  ggplot(aes(x=x, y = .)) +
  geom_line() +
  scale_x_continuous(limits = c(0,1), breaks = scales::pretty_breaks(5)) +
  labs(x = "probability", y = "density", title = "Beta(1, 5)")
pp2


### 2. Beta(1, 3)
pp3 = dbeta(seq(0,1,0.01), 1, 3) %>% 
  tibble() %>% mutate(x = seq(0,1,0.01)) %>% 
  ggplot(aes(x=x, y = .)) +
  geom_line() +
  scale_x_continuous(limits = c(0,1), breaks = scales::pretty_breaks(5)) +
  labs(x = "probability", y = "density", title = "Beta(1, 3)")
pp3



### To one plot
pp0 + pp1+ pp2 + pp3 

ggsave(plot = last_plot(), filename = "plots/priors_tested.png",
       width = 6.5, height = 5, dpi = 300)


## Figure out what the sensitivity analysis will be on?
## N TRs needed to clear median, 80th, 95th of P(A) > 0.95


### Run analysis for each Pprior
i = 100


# Define function to run analysis
SA = function(i, Pprior) {
  
  ## get raw results from each Bayesian sample
  df_raw = 
    map_df(1:i, function(x){
    rvals = Pprior # rbeta(i, 1, 0.2)
    det_df %>% 
      mutate(p_abs = p_abs(p_det, rvals[x])) %>% 
      mutate(i = x) %>% 
      select(i, ntr, p_abs)
    })  
  
  ## get summary stats
  dfs = 
    df_raw %>% 
    group_by(ntr) %>% 
    summarise(q2.5 = quantile(p_abs, 0.025),
              q10 = quantile(p_abs, 0.1),
              med = median(p_abs),
              q90 = quantile(p_abs, 0.9),
              q97.5 = quantile(p_abs, 0.975)) %>% 
    ungroup() 
  
  return(dfs)
}


## run
sa1 = SA(i, rbeta(i, 0.25, 4)) %>% mutate(Pprior = "Beta(0.25, 4)")
sa2 = SA(i, runif(i, 0.01, 0.1)) %>% mutate(Pprior = "Uniform(0.01, 0.10)")
sa3 = SA(i, rbeta(i, 1, 5)) %>% mutate(Pprior = "Beta(1, 5)")
sa4 = SA(i, rbeta(i, 1, 3)) %>% mutate(Pprior = "Beta(1, 3)")



# Combine all
sa_all = sa1 %>% bind_rows(sa2) %>% bind_rows(sa3) %>% 
  bind_rows(sa4)



# Get difference to reference
sa_min = 
  sa_all %>% 
  select(-q90, -q97.5) %>% 
  mutate(Pprior = 
           forcats::fct_relevel(Pprior, 
                   c("Beta(0.25, 4)", "Uniform(0.01, 0.10)", 
                     "Beta(1, 5)", "Beta(1, 3)")
                     )) %>% 
  pivot_longer(q2.5:med) %>% 
  mutate(over95 = value > 0.95) %>% 
  filter(over95) %>% 
  group_by(name, Pprior) %>% 
  summarise(min_trs = min(ntr)) %>% 
  ungroup() %>% 
  complete(name, Pprior)
  
sa_diff = 
  sa_min %>% 
  group_by(name) %>% 
  mutate(dif = min_trs - first(min_trs, na_rm = FALSE)) %>% 
  ungroup() %>% 
  mutate(pos = dif > 0)



# Plot
sa_diff %>%
  filter(Pprior != "Beta(0.25, 4)") %>% 
  ggplot(aes(x = name, y = dif, fill = pos)) + 
  geom_bar(stat = "identity", color = "black", size = 0.5) + 
  scale_y_continuous(breaks = scales::pretty_breaks(6)) + 
  scale_x_discrete(#guide = guide_axis(n.dodge = 2),
                   labels = c("q2.5" = "2.5th\nQuantile", 
                              "q10" =  "10th\nQuantile", 
                              "med" =  "Median")) + 
  scale_fill_manual(values = c("#eed5d2", "#b0e0e6")) +
  guides(fill = "none") + 
  labs(x = "P(absence | no detections) threshold", 
       y ="Difference in the number of qPCR techincal replicates needed\ncompared to Beta(0.25, 4)", 
       subtitle = ""
  )+ 
  facet_grid(~Pprior #, 
             # labeller = labeller(p = c("p2.5" = "Lower 5th", 
             #                           "p10" =  "Lower 20th", 
             #                           "med" =  "0.50")
                                 # )
  ) + 
  theme_bw() + 
  theme() # axis.text.x = element_text(angle = 30, vjust=0.75)


ggsave(plot = last_plot(), filename = "plots/prior_sens_analysis.png",
       width = 6.5, height = 5, dpi = 300)






# library(patchwork)
# p1 + p2
# 
# 
# 
# # p absnece
# p_abs = function(S, Pprior){
#   A = 1 - ( ((1-S)*Pprior) / ( ((1-S)*Pprior) + 1 - Pprior ) )
#   return(A)
# }
# 
# 
# vals = expand.grid(seq(0,1,0.1), seq(0,1,0.1))
# df = map_df(1:nrow(vals), 
#             ~tibble(S = vals$Var1[.x], Pprior = vals$Var2[.x], A = p_abs(vals$Var1[.x], vals$Var2[.x])) )
# 
# df %>% 
#   ggplot(aes(x = S, y = A, color = Pprior, group = Pprior)) + geom_line() + 
#   geom_abline(aes(slope = 1, intercept = 0))
# 
# df %>% 
#   ggplot(aes(x = S, fill = A, y = Pprior)) + 
#   geom_tile() + 
#   scale_fill_distiller()
