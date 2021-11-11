# spline


spline(data_plot$age_bin, data_plot$sample_frequency) %>% as.data.frame() %>% 
  ggplot(aes(x=x, y=y))+
  geom_line()


b <- data_plot %>% filter(IDs == "All individuals")
b <- spline(b$age_bin, b$sample_frequency) %>% as.data.frame() %>% mutate(IDs = "All individuals")
a <- data_plot %>% filter(IDs == "2-25328714-C-G")
a <- spline(a$age_bin, a$sample_frequency) %>% as.data.frame() %>% mutate(IDs = "2-25328714-C-G") %>% 
  bind_rows(b)
# a %>% 
#   ggplot(aes(x=x, y=y, color = IDs))+
#   geom_line()
# 
# p <- a %>% 
#   ggplot(aes(x=x, y=y, color = IDs))+
#   geom_line()+
#   ylim(0, max(a$y))

a %>% 
  ggplot(aes(x=x, y=y, color = IDs))+
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
              span = 0.45)+
  ylim(0, max(a$y))

p <- a %>% 
  ggplot(aes(x=x, y=y, color = IDs))+
  geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
              span = 0.45)+
  ylim(0, max(a$y))

p_ <- layer_data(p, 1)
p_

t.test(p_ %>% filter(group == 1) %>% select(y), p_ %>% filter(group == 2) %>% select(y))$p.value







