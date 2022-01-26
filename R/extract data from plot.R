data_plot <- gnomad_decoded4 %>% 
  distinct(IDs, age_bin, .keep_all = TRUE) %>% 
  filter(IDs == "2-25328714-C-G" | IDs == "All individuals") %>% 
  select("IDs", "nbr_individuals", "age_bin", "sample_count",starts_with("sample_density"),
         AN_samples_count, AN,
         sample_frequency, age_bin, variant_frequency)

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
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess", span = 0.45)+
  ylim(0, max(data_plot$sample_frequency))
ggplot_build(p1)

p1_1 <- layer_data(p1, 1)
p1_1

t.test(p1_1 %>% filter(group == 1) %>% select(y), p1_1 %>% filter(group == 2) %>% select(y))

###################### Using sampling
x <- seq(15, 90,5)
y1 <- rnorm(length(x), 0.05, 0.02)
y2 <- rnorm(length(x), 0.06, 0.02)
n1 <- sample(x, 50, prob = y1, replace = TRUE)
n2 <- sample(x, 50, prob = y2, replace = TRUE)
plot(density(n1))
t.test(n1, n2)
y3 <- rexp(length(x),0.2)
plot(density(y3))
n3 <- sample(x, 50, prob = y3, replace = TRUE)
kurtosis(n3) # WRONG

kurtosis(y3)

y4 <- rnorm(length(x), mean = 0, sd = 4)
n1 <- sample(x, 50, prob = y4, replace = TRUE)

x <- seq(-4, 4, length=100)
y <- dunif(x, min = -3, max = 3)
kurtosis(y)
plot(x, y, type = 'l')
plot(density(y))

success <- 0:20
a <- dbinom(success, size=20, prob=.3)
plot(success, dbinom(success, size=20, prob=.3),type='h')
kurtosis(a)

plot(density(y4))
kurtosis(y4)
# The values for asymmetry and kurtosis between -2 and +2 are considered 
# acceptable in order to prove normal univariate distribution (George & Mallery, 2010). 
# (2010) and Bryne (2010) argued that data is considered to be normal if 
# skewness is between ‐2 to +2 and kurtosis is between ‐7 to +7.
# Kurtosis is a measure of the combined sizes of the two tails. ... 
# If the kurtosis is greater than 3, then the dataset has heavier tails than a 
# normal distribution (more in the tails). If the kurtosis is less than 3, then 
# the dataset has lighter tails than a normal distribution (less in the tails).

p <- data_plot %>% 
  ggplot(aes(x = center, y= sample_frequency, color= IDs))+
  # ggplot(aes(x = age_bin, y= sample_frequency, color= IDs))+
  # geom_bar(stat = "identity", aes(alpha= IDs))
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
              # method.args=list(family=quasibinomial)#,
              span = 0.6
  )+
  ylim(0, max(data_plot$sample_frequency))
p_ <- layer_data(p, 1) %>% 
  mutate(colour =  as.factor(colour)) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))

ref <- p_ %>% filter(group == 2)
ref <- sample(ref$x, 180, prob = ref$y, replace = TRUE)

variant <- p_ %>% filter(group == 1)
variant <- sample(variant$x, 180, prob = variant$y, replace = TRUE)

t.test(ref, variant)


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

