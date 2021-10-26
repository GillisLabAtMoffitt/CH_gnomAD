data_plot <- gnomad_decoded %>% 
  distinct(IDs, age_bin, .keep_all = TRUE) %>% 
  filter(IDs == "2-25234330-T-C" | IDs == "All individuals") %>% 
  select("IDs", "nbr_individuals", "age_bin", "sample_count",starts_with("sample_density"), sample_frequency, center)

data_plot %>% 
  ggplot(aes(x = age_bin, y= sample_count, fill= IDs))+
  geom_bar(stat = "identity", alpha= 0.5)
data_plot %>% 
  ggplot(aes(x = age_bin, y= sample_density, fill= IDs))+
  geom_bar(stat = "identity", alpha= 0.5)
data_plot %>% 
  ggplot(aes(x = age_bin, y= sample_density0, fill= IDs))+
  geom_bar(stat = "identity", alpha= 0.5)
data_plot %>% 
  ggplot(aes(x = age_bin, y= sample_density1, fill= IDs))+
  geom_bar(stat = "identity", alpha= 0.5)
data_plot %>% 
  ggplot(aes(x = age_bin, y= sample_frequency, fill= IDs))+
  geom_bar(stat = "identity", alpha= 0.5)

data_plot %>% 
  ggplot(aes(x = center, y= sample_frequency, color= IDs))+
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess", span = 0.6)

p1 <- data_plot %>% 
  ggplot(aes(x = center, y= sample_frequency, color= IDs))+
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess", span = 0.6)
ggplot_build(p1)

p1_1 <- layer_data(p1, 1)
p1_1

t.test(p1_1 %>% filter(group == 1) %>% select(y), p1_1 %>% filter(group == 2) %>% select(y))





########################################################
p1 <- data_plot %>% 
  ggplot(aes(x = center, y = sample_density)) + 
  geom_bar(stat = "identity")+
  theme_minimal()
p1

p1_1 <- layer_data(p1, 1)
p1_1

p1+ geom_smooth(data=p1_1, aes(x= x, y= ymax), se = FALSE, method = "loess", span = 0.6)
p1+ geom_density(data=p1_1, aes(x= xmin, y= ymax), stat = "identity")

p3<- p1+ geom_smooth(data=p1_1, aes(x= x, y= ymax), se = FALSE, method = "loess", span = 0.6)
layer_data(p3, 1) # had the same data

ref <- tibble(age_bin = c("<30", "30-35","35-40","40-45","45-50","50-55","55-60",
                          "60-65","65-70","70-75","75-80", ">80"),
              y = c(2547, 3423, 4546, 8487, 10355, 12693, 11933, 10534, 8882, 5991, 4136, 1935)) %>% 
  mutate(age_bin = factor(age_bin, 
                          levels = c("<30", "30-35","35-40","40-45","45-50","50-55",
                                     "55-60","60-65","65-70","70-75","75-80", ">80"))) %>% 
  mutate(center = case_when(age_bin == "<30" ~ 15,
                            age_bin == "30-35" ~ 32.5,
                            age_bin == "35-40" ~ 47.5,
                            age_bin == "40-45" ~ 42.5,
                            age_bin == "45-50" ~ 57.5,
                            age_bin == "50-55" ~ 52.5,
                            age_bin == "55-60" ~ 67.5,
                            age_bin == "60-65" ~ 62.5,
                            age_bin == "65-70" ~ 77.5,
                            age_bin == "70-75" ~ 72.5,
                            age_bin == "75-80" ~ 87.5,
                            TRUE ~ 90
  )) %>% 
  mutate(sum = sum(y)) %>% 
  mutate(density = y / sum)

ggplot() + 
  geom_col(data = data_plot ,aes(x= center, y= sample_count), 
           position = "identity", 
           fill = "darkblue", 
           alpha = 0.6)+
  geom_smooth(data = p1_1 ,aes(x= x, y= ymax), 
              stat = "identity", 
              color = "darkblue", 
              alpha = 0.6)+
  theme_minimal()+
  facet_wrap(.~ IDs, ncol = 3, scales = "free_y")+
  geom_smooth(data =ref, aes(x =center, y = density), stat = "identity", color = "yellow", alpha = 0.5)


ggplot() + geom_smooth(data=p1_1, aes(x= x, y= ymax), se = FALSE, method = "loess", span = 0.6)+
  geom_smooth(data =ref, aes(x =center, y = density), color = "yellow", alpha = 0.5, se = FALSE, method = "loess", span = 0.6)
p5 <- ggplot() + geom_smooth(data=p1_1, aes(x= x, y= ymax), se = FALSE, method = "loess", span = 0.6)
p6 <- ggplot() + geom_smooth(data =ref, aes(x =center, y = density), color = "yellow", alpha = 0.5, se = FALSE, method = "loess", span = 0.6)
p5_1 <- layer_data(p5, 1)
p6_1 <- layer_data(p6, 1)
t.test(p5_1$y, p6_1$y) ################### Here









#############

  
p1 <- data_plot %>% 
  ggplot(aes(x = age_bin, y = sample_count)) + 
  geom_bar(stat = "identity")+
  theme_minimal()
p1

ggplot_build(p1)
ggplot_build(p1)$plot$data
ggplot_build(p1)$plot$layers
ggplot_build(p1)$plot$mapping
ggplot_build(p1)$plot$coordinates

p1_1 <- layer_data(p1, 1)

p1_1 %>% 
  ggplot(aes(x= xmin, y= ymax))+
  geom_smooth()+
  geom_bar(aes(p1), stat = "identity")

p1+ geom_smooth(data=p1_1, aes(x= xmin, y= ymax))
p1+ geom_density(data=p1_1, aes(x= xmin, y= ymax), stat = "identity")
###########
gnomad_decoded1 %>% 
  ggplot(aes(x= center, y= count_in_bin))+
  geom_density(stat = "identity")+
  theme_minimal()


gnomad_decoded1 %>% 
  ggplot(aes(x = age_bin, y = count_in_bin)) + 
  geom_bar(stat = "identity")+
  theme_minimal()+
  facet_wrap(.~ID)

p2 <- ggplot(gnomad_decoded1, aes(x = center, y = count_in_bin)) + stat_smooth(aes(y = count_in_bin,x = center), method = "gam",se = FALSE,formula = y ~ s(x, k = 7))
p2

ggplot_build(p2)
ggplot_build(p2)$plot$data
ggplot_build(p2)$plot$layers
ggplot_build(p2)$plot$mapping
ggplot_build(p2)$plot$coordinates
layer_data(p2, 1)

p2_build = ggplot_build(p2)
p2_fill <- data_frame(
  x = p2_build$data[[1]]$x,
  y = p2_build$data[[1]]$y#,
  # group = factor(p2_build$data[[1]]$group, levels = c(1,2), labels = c("apples","bananas"))
)

ggplot(gnomad_decoded1, aes(left, count_in_bin))+
  stat_smooth(aes(y = count_in_bin,x = left), method = "gam",se = FALSE, formula = y ~ s(x, k = 7))+
  # geom_area(data = p2_fill[p2_fill$group == "apples", ], 
  #           aes(x=x, y=y), fill =  "#F8766D", alpha = 0.2, inherit.aes = F)+
  geom_area(data = p2_fill, #[p2_fill$group == "bananas", ], 
            aes(x=x, y=y), fill = "#00BFC4", alpha = 0.2, inherit.aes = F)+
  theme_classic()












#################################################################
############# ____________Infering data____________ #############
#################################################################
# https://www.r-bloggers.com/2019/08/inferring-a-continuous-distribution-from-binned-data-by-ellis2013nz/

gnomad_decoded1 %>% 
  filter(count_in_bin != 0) %>% 
  group_by(ID, age_bin) %>%
  summarise(freq = n())

volumes_bin <- dplyr::select(gnomad_decoded1, left, right) %>% 
  as.data.frame()

library(fitdistrplus)
fitted_distribution_gamma <- fitdistcens(volumes_bin, "gamma")
# overall fit:
summary(fitted_distribution_gamma)
# estimated mean
fitted_distribution_gamma$estimate["shape"] / fitted_distribution_gamma$estimate["rate"]

ggplot(volumes) +
  geom_density(aes(x = volume)) +
  stat_function(fun = dgamma, args = fitted_distribution_gamma$estimate, colour = "blue") +
  annotate("text", x = 700, y = 0.0027, label = "Blue line shows modelled distribution; black is density of the actual data.")

# bin_based_mean (358 - very wrong)
gnomad_decoded1 %>%
  mutate(mid = (left + replace_na(right, 2000)) / 2) %>%
  summarise(crude_mean = mean(mid)) %>%
  pull(crude_mean)
#################################################################



#################################################################
############# ____________creating bins____________ #############
#################################################################
library(ggplot2)
qplot(gnomad_decoded1, data=cbind(age_bin,count_in_bin), weight=count_in_bin, geom="histogram")

gnomad_decoded1 %>% qplot(aes(x= age_bin, y= count_in_bin, weight=count_in_bin), geom="histogram")

gnomad_decoded1 %>% 
  ggplot(aes(x = age_bin, y = count_in_bin, fill = ID)) + 
  geom_bar(stat = "identity")+
  theme_minimal()+
  facet_wrap(.~ID)

gnomad_decoded1$binRange <- str_extract(gnomad_decoded1$age_bin, "[0-9]-*[0-9]+")

# split the data using the , into to columns:
# one for the start-point and one for the end-point
library(splitstackshape)
df <- cSplit(gnomad_decoded1, "binRange")

# plot it, you actually dont need the second column
ggplot(df, aes(x = binRange_1, y = count_in_bin, width = 0.025)) +
  geom_bar(stat = "identity"#, breaks=seq(0,0.125, by=0.025)
  )
#################################################################

