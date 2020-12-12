library(data.table)
library(dplyr)
library(ggplot2)
covid <- fread("D:/Github/coronavirus/csv/coronavirus.csv")
summary(covid)
#Format data
covid[,date := as.Date(date)][, province := factor(province)][, country := factor(country)][, type := factor(type)]
#Find the country with the most cases
covid[type == 'confirmed', sum(cases), by = 'country'][order(-V1)]
#Summary the new cases during the last 24h
covid[date == max(date), sum(cases), by = c('country', 'type')] %>% 
  dcast.data.table(country ~ type) %>% 
  .[order(-confirmed)]
#Stacked line plot (area chat) the totle cases worldwide
covid_world <- covid[, sum(cases), by = c('type', 'date')] %>% 
  dcast.data.table(date ~ type, value.var = "V1") %>% 
  .[, active := confirmed - death - recovered] %>% 
  .[, active_total := cumsum(active)] %>% 
  .[, recovered_total := cumsum(recovered)] %>% 
  .[, death_total := cumsum(death)]
plot <- covid_world[, c(1,7,8,6)] %>% 
  melt.data.table(id = 1, measure.vars = c(2,3,4), variable.name = 'type', value.name = 'cases', value.factor = T) %>% 
  ggplot(aes(x = date, y = cases/10^6, fill = type)) +
  geom_area() +
  labs(fill = "Type", x = "Date", y = "Cases, million") +
  scale_fill_manual(labels = c('Recovered', 'Death', 'Active'), values = c('#af8dc3','#f7f7f7','#7fbf7b'))
plot  
#Treemap plot
library(treemapify)
covid[type == 'confirmed', sum(cases), by = 'country'][order(-V1)] %>% 
  ggplot(aes(area = V1, label = country)) +
  geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre", grow = T)
  theme(legend.position = 'none')
  
