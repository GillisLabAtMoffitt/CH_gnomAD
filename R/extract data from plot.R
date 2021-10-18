
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

