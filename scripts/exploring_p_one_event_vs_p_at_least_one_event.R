

library(tidyverse)
library(directlabels)
library(patchwork)


p_event = function(n_events, p_one_event){
  return(1 - (1-p_one_event)^n_events)
}

p_event(1:100, 0.01)

p_event_multi = function(num_events, p_one_event){
  temp = list(); d=1
  for(p in p_one_event){
    temp[[d]] = tibble(n_events = 1:num_events, 
                       p_one_event = rep(p, num_events),
                       p_at_least_one = p_event(1:num_events, p))
    d = d + 1
  }
  return(bind_rows(temp))
}


df = p_event_multi(100, seq(0.001,0.25,0.001))

p1  = 
  df %>% 
  ggplot(aes(y = n_events, x = p_one_event, fill = p_at_least_one)) +
  geom_tile() + 
  scale_fill_fermenter(breaks=seq(0,0.8,0.1)) + 
  labs(x = "P(success of any given event)", 
       fill = "P(at leat one\nevent successful)", y = "N events")
p1


p2 = 
  df %>% 
  filter(n_events %in% c(1:10,seq(20, 100, 10))) %>%  
  ggplot(aes(x = p_one_event, y = p_at_least_one, 
             group = n_events, color = n_events)) + 
  geom_line() + 
  # geom_abline(aes(slope = 1, intercept = 0), size = 2, linetype = "dashed") + 
  scale_y_continuous(breaks = scales::pretty_breaks(4)) + 
  # scale_x_continuous(trans="log10") +
  guides(color = "none") + 
  directlabels::geom_dl(aes(label = n_events), method = "last.qp") + 
  labs(x = "P(success of any given event)", 
       y = "P(at leat one event successful)", color = "N\nevents")
p2

p3 = 
  df %>% 
  filter(n_events %in% c(1:10,seq(20, 100, 10))) %>%  
  ggplot(aes(x = p_one_event, y = p_at_least_one, 
             group = n_events, color = n_events)) + 
  geom_line() + 
  # geom_abline(aes(slope = 1, intercept = 0), size = 2, linetype = "dashed") + 
  scale_y_continuous(breaks = scales::pretty_breaks(4)) + 
  scale_x_continuous(trans="log10") + 
  # guides(color = "none") + 
  directlabels::geom_dl(aes(label = n_events), method = "last.qp") + 
  labs(x = "P(success of any given event)", 
       y = "P(at leat one event successful)", color = "N\nevents")
p3

(p2 + p3)

