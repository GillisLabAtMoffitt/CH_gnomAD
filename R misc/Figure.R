# Load data
library(tidyverse)
library(gtsummary)
theme_gtsummary_compact()
theme_set(theme_classic(base_size = 10))

# List pf variants----
library(tidyverse)
CH_var <- 
  read.csv(paste0(here::here(), "/Pop freq data for CH.csv")) %>% 
  mutate(`Is a CH variant` = "Yes")
SM_var <- 
  read.csv(paste0(here::here(), "/Pop freq data for CM.csv")) %>% 
  mutate(`SM variant IDs` = IDs)

a <- left_join(SM_var, CH_var, by = "IDs") %>% 
  mutate(`Is a CH variant` = case_when(
    `Is a CH variant` == "Yes"      ~ "Yes",
    TRUE                            ~ "No"
  )) %>% 
  select(`SM variant IDs`, `Is a CH variant`)

write_csv(a, "List of gnomAD identified SM and CH IDs.csv")


# Figure 4 A-D----
library(tidyverse)
library(ggupset)
theme_set(theme_classic(base_size = 12))



# Figure 4 B-D
CH_ <- 
  read.delim(paste0(here::here(), "/PAFs of CH variants with zero values.txt"), sep = " ")

fig4_B <- # CH / genes
  CH_ %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(!is.na(SYMBOL)) %>%
  select(IDs, SYMBOL, 
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  pivot_longer(cols = -c(IDs, SYMBOL)) %>% 
  filter(value > 0) %>% 
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         ancestry = case_when(
           name == "nfe"           ~ "European (non-Finnish)",
           name == "afr"           ~ "African/African American",
           name == "amr"           ~ "Latino/Admixed American",
           name == "eas"           ~ "East Asian",
           name == "sas"           ~ "South Asian",
           name == "asj"           ~ "Ashkenazi Jewish",
           name == "oth"           ~ "Other",
           name == "fin"           ~ "European (Finnish)",
         ))
b <- fig4_B %>% 
  distinct(SYMBOL, ancestry) %>% 
  group_by(SYMBOL) %>%
  summarize(ancestrys = list(ancestry))

# pdf("Figure 4 B upset CH-gene.pdf",
#     width = 10,
#     height = 5)
b %>% 
  ggplot(aes(x = ancestrys)) +
  geom_bar() +
  scale_x_upset()
# dev.off()

c <-
  b %>% 
  ggplot(aes(x = ancestrys)) +
  geom_bar(fill = c("black", "#00BFC4", "#00A9FF", 
                    "#C77CFF", "#7CAE00", "#F8766D", 
                    "#00BE67", "#CD9600", "#FF61CC")
    )+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(intersections = list(c("African/African American"),
                                     c("Ashkenazi Jewish"),
                                     c("East Asian"),
                                     c("European (Finnish)"),
                                     c("European (non-Finnish)"),
                                     c("Latino/Admixed American"),
                                     c("South Asian"),
                                     c("Other"),
                                     c("Other",
                                       "South Asian",
                                       "Latino/Admixed American",
                                       "European (non-Finnish)",
                                       "European (Finnish)",
                                       "East Asian",
                                       "Ashkenazi Jewish",
                                       "African/African American")
  ), sets =  c("African/African American",
               "Ashkenazi Jewish",
               "East Asian",
               "European (Finnish)",
               "European (non-Finnish)",
               "Latino/Admixed American",
               "South Asian",
               "Other"))+
  scale_y_continuous(breaks = NULL, lim = c(0, 1050), name = "")+
  labs(x = NULL)

  # scale_color_manual(values= c("African/African American" = "red", "Non-Action" = "black", "not observed" = "lightgrey")) +
c

ggsave("Figure 4 B upset CH-gene2.pdf", plot = c,
       width = 6,
       height = 5,
       dpi = 600, limitsize = FALSE)

# Figure 4 D upset CH-variants

b <- fig4_B %>% 
  distinct(IDs, ancestry) %>% 
  group_by(IDs) %>%
  summarize(ancestrys = list(ancestry))

c <-
  b %>% 
  ggplot(aes(x = ancestrys)) +
  geom_bar(fill = c("#00BFC4", "#00A9FF", 
                    "#C77CFF", "#7CAE00", "#F8766D",
                    "#00BE67", "#FF61CC", "#CD9600", "black")
  )+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(intersections = list(c("African/African American"),
                                     c("Ashkenazi Jewish"),
                                     c("East Asian"),
                                     c("European (Finnish)"),
                                     c("European (non-Finnish)"),
                                     c("Latino/Admixed American"),
                                     c("South Asian"),
                                     c("Other"),
                                     c("Other",
                                       "South Asian",
                                       "Latino/Admixed American",
                                       "European (non-Finnish)",
                                       "European (Finnish)",
                                       "East Asian",
                                       "Ashkenazi Jewish",
                                       "African/African American")
  ), sets =  c("African/African American",
               "Ashkenazi Jewish",
               "East Asian",
               "European (Finnish)",
               "European (non-Finnish)",
               "Latino/Admixed American",
               "South Asian",
               "Other"))+
  scale_y_continuous(breaks = NULL, lim = c(0, 22000), name = "")+
  labs(x = NULL)

# scale_color_manual(values= c("African/African American" = "red", "Non-Action" = "black", "not observed" = "lightgrey")) +
c

ggsave("Figure 4 D upset CH-variants2.pdf", plot = c,
       width = 6,
       height = 5,
       dpi = 600, limitsize = FALSE)








# Figure 4 A-C SM----
SM_ <- 
  read.delim(paste0(here::here(), "/PAFs of SM variants with zero values.txt"), sep = " ")

fig4_B <- # SM / genes
  SM_ %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  filter(!is.na(SYMBOL)) %>%
  select(IDs, SYMBOL, 
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  pivot_longer(cols = -c(IDs, SYMBOL)) %>% 
  filter(value > 0) %>% 
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         ancestry = case_when(
           name == "nfe"           ~ "European (non-Finnish)",
           name == "afr"           ~ "African/African American",
           name == "amr"           ~ "Latino/Admixed American",
           name == "eas"           ~ "East Asian",
           name == "sas"           ~ "South Asian",
           name == "asj"           ~ "Ashkenazi Jewish",
           name == "oth"           ~ "Other",
           name == "fin"           ~ "European (Finnish)",
         ))
b <- fig4_B %>% 
  distinct(SYMBOL, ancestry) %>% 
  group_by(SYMBOL) %>%
  summarize(ancestrys = list(ancestry))

c <-
  b %>% 
  ggplot(aes(x = ancestrys)) +
  geom_bar(fill = c("black", "#00BFC4", "#C77CFF", 
                    "#00A9FF", "#7CAE00", "#F8766D", 
                    "#00BE67", "#CD9600", "#C77CFF")
  )+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(intersections = list(c("African/African American"),
                                     c("Ashkenazi Jewish"),
                                     c("East Asian"),
                                     c("European (Finnish)"),
                                     c("European (non-Finnish)"),
                                     c("Latino/Admixed American"),
                                     c("South Asian"),
                                     c("Other"),
                                     c("Other",
                                       "South Asian",
                                       "Latino/Admixed American",
                                       "European (non-Finnish)",
                                       "European (Finnish)",
                                       "East Asian",
                                       "Ashkenazi Jewish",
                                       "African/African American")
  ), sets =  c("African/African American",
               "Ashkenazi Jewish",
               "East Asian",
               "European (Finnish)",
               "European (non-Finnish)",
               "Latino/Admixed American",
               "South Asian",
               "Other"))+
  scale_y_continuous(breaks = NULL, lim = c(0, 9800), name = "") +
  labs(x = NULL)

c

ggsave("Figure 4 A upset SM-gene.pdf", plot = c,
       width = 6,
       height = 5,
       dpi = 600, limitsize = FALSE)

# Figure 4 D upset SM-variants
b <- fig4_B %>% 
  distinct(IDs, ancestry) %>% 
  group_by(IDs) %>%
  summarize(ancestrys = list(ancestry))

c <-
  b %>% 
  ggplot(aes(x = ancestrys)) +
  geom_bar(fill = c("#00BFC4", "#C77CFF", "#00A9FF", 
                    "#7CAE00", "#F8766D", "black", 
                    "#00BE67", "#FF61CC", "#CD9600")
  )+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(intersections = list(c("African/African American"),
                                     c("Ashkenazi Jewish"),
                                     c("East Asian"),
                                     c("European (Finnish)"),
                                     c("European (non-Finnish)"),
                                     c("Latino/Admixed American"),
                                     c("South Asian"),
                                     c("Other"),
                                     c("Other",
                                       "South Asian",
                                       "Latino/Admixed American",
                                       "European (non-Finnish)",
                                       "European (Finnish)",
                                       "East Asian",
                                       "Ashkenazi Jewish",
                                       "African/African American")
  ), sets =  c("African/African American",
               "Ashkenazi Jewish",
               "East Asian",
               "European (Finnish)",
               "European (non-Finnish)",
               "Latino/Admixed American",
               "South Asian",
               "Other"))+
  scale_y_continuous(breaks = NULL, lim = c(0, 95000), name = "")+
  labs(x = NULL)

c

ggsave("Figure 4 C upset SM-variants.pdf", plot = c,
       width = 6,
       height = 5,
       dpi = 600, limitsize = FALSE)

# Using upsetR----
library(UpSetR)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")
movies

b <- fig4_B %>% 
  distinct(IDs, ancestry) #%>% 
  # group_by(IDs) %>%
  # summarize(ancestrys = list(ancestry))
library(fastDummies)
dummy_cols(b$ancestry)
b

dum <- dummy_cols(.data = b, select_columns = "ancestry", remove_selected_columns = TRUE, omit_colname_prefix = TRUE) %>% 
  group_by(IDs) %>% 
  summarise(across(everything(), sum)) %>% ungroup() %>% 
  as.data.frame()

# UpSetR::upset(dum, nintersects = 10, 
#               empty.intersections = "on", order.by = "freq",
#               nsets = 10, sets =  c("African/African American",
#                                                          "Ashkenazi Jewish",
#                                                          "East Asian",
#                                                          "ancestry_European (Finnish)",
#                                                          "ancestry_European (non-Finnish)",
#                                                          "Latino/Admixed American",
#                                                          "South Asian",
#                                                          "Other")
#               )


input <- c(
  "African/African American" = 104,
  "Ashkenazi Jewish" = 9,
  "East Asian" = 126,
  "European (Finnish)" = 33,
  "European (non-Finnish)" = 588,
  "Latino/Admixed American" = 197,
  "South Asian" = 136,
  "Other" = 7,
  "African/African American&Ashkenazi Jewish&East Asian&European (Finnish)&European (non-Finnish)&Latino/Admixed American&South Asian&Other" = 953)
# p <- 
pdf("UpsetR CHgene.pdf", onefile=FALSE,
    width = 7,
    height = 4)
  upset(fromExpression(input), 
      nintersects = 40, 
      nsets = 10, 
      order.by = c("freq", "degree"),# decreasing = c(TRUE,TRUE),
      # decreasing = T, 
      # mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1,
      set_size.show = TRUE,
      sets =  c("African/African American",
                "Ashkenazi Jewish",
                "East Asian",
                "European (Finnish)",
                "European (non-Finnish)",
                "Latino/Admixed American",
                "South Asian",
                "Other"),
      queries = list(list(query = intersects, params = list("African/African American",
                                                            "Ashkenazi Jewish",
                                                            "East Asian",
                                                            "European (Finnish)",
                                                            "European (non-Finnish)",
                                                            "Latino/Admixed American",
                                                            "South Asian",
                                                            "Other"), color = "blue", active = T))
)
dev.off()

ups <- upset(fromExpression(input), 
      nintersects = 40, 
      nsets = 10, 
      order.by = c("freq", "degree"),# decreasing = c(TRUE,TRUE),
      # decreasing = T, 
      # mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1,
      set_size.show = TRUE,
      sets =  c("African/African American",
                "Ashkenazi Jewish",
                "East Asian",
                "European (Finnish)",
                "European (non-Finnish)",
                "Latino/Admixed American",
                "South Asian",
                "Other"),
      queries = list(list(query = intersects, params = list("African/African American",
                                                            "Ashkenazi Jewish",
                                                            "East Asian",
                                                            "European (Finnish)",
                                                            "European (non-Finnish)",
                                                            "Latino/Admixed American",
                                                            "South Asian",
                                                            "Other"), color = "blue", active = T))
)
library(grid)
library(gridExtra)
skip_set_size_plot <- function(ups) {
  main <- ups$Main_bar
  ## find panel grob
  panel_idx <- grep("panel", main$layout$name, fixed = TRUE)
  ## find text grob
  text_idx <- which(
    vapply(main$grobs[[panel_idx]]$children, 
           \(x) inherits(x, "text"), 
           logical(1)))
  tG <- main$grobs[[panel_idx]]$children[[text_idx]]
  # tG$label <- paste0(tG$label, " (",
  #                    scales::label_percent(0.1)(as.numeric(tG$label) / 
  #                                                 sum(as.numeric(tG$label))),
  #                    ")")
  main$grobs[[panel_idx]]$children[[text_idx]] <- tG
  grid::grid.newpage()
  # grid::grid.draw(gridExtra::arrangeGrob(main, ups$Matrix, heights = ups$mb.ratio, widths = c(2,0.2)))
  ggpubr::ggarrange(plotlist = list(main, ups$Matrix), ncol = 1, align = "v", vjust = 1, hjust = 1)
}
pdf("UpsetR CHgene without size plot.pdf", onefile=FALSE,
    width = 7,
    height = 4)
skip_set_size_plot(ups)
dev.off()

# ggsave("UpsetR CHgene.pdf", plot = p,
#        width = 7,
#        height = 4,
#        dpi = 600, limitsize = FALSE)

# Figure 5 for table 3 ACs number-----
library(tidyverse)
theme_set(theme_classic(base_size = 10))
dat <- tibble::tribble(
  ~Symbol, ~`African/African.American`, ~Ashkenazi.Jewish, ~East.Asian, ~`European.(Finnish)`, ~`European.(non-Finnish)`, ~`Latino/Admixed.American`, ~Other, ~South.Asian, ~All.Populations, ~`African/African.American`, ~Ashkenazi.Jewish, ~East.Asian, ~`European.(Finnish)`, ~`European.(non-Finnish)`, ~`Latino/Admixed.American`, ~Other, ~South.Asian, ~All.populations,
  "DNMT3A",                         39L,               13L,         24L,                   71L,                      161L,                        22L,    12L,          23L,             365L,                         26L,               11L,         21L,                   18L,                      101L,                         8L,     8L,          21L,             214L,
  "CUX1",                          3L,                5L,          4L,                   11L,                       45L,                        21L,     1L,          32L,             122L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               0L,
  "SF3B1",                         29L,                5L,          2L,                   16L,                       44L,                         7L,     4L,           7L,             114L,                          1L,                3L,          1L,                    2L,                       29L,                         6L,     1L,           3L,              46L,
  "NF1",                         41L,                2L,         17L,                    5L,                       34L,                         4L,     0L,           0L,             103L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               0L,
  "TP53",                         11L,                5L,          4L,                    9L,                       58L,                         3L,     1L,           7L,              98L,                          0L,                1L,          0L,                    0L,                        6L,                         0L,     0L,           0L,               7L,
  "JAK2",                          5L,                6L,          3L,                   12L,                       45L,                         6L,     1L,          10L,              88L,                          5L,                6L,          3L,                   12L,                       45L,                         6L,     1L,          10L,              88L,
  "PDS5B",                          5L,                3L,          0L,                    1L,                       44L,                         4L,     0L,           1L,              58L,                          4L,                3L,          0L,                    1L,                       30L,                         3L,     0L,           1L,              42L,
  "MFSD11",                          5L,                2L,          2L,                    0L,                       37L,                         1L,     1L,           2L,              50L,                          5L,                2L,          2L,                    0L,                       37L,                         1L,     1L,           2L,              50L,
  "EZH2",                          9L,                0L,         11L,                    7L,                       22L,                         1L,     0L,           0L,              50L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               0L,
  "ASXL1",                          2L,                2L,          5L,                    4L,                       26L,                         1L,     1L,           0L,              41L,                          2L,                2L,          2L,                    2L,                       16L,                         1L,     0L,           0L,              25L,
  "GNB1",                          0L,                1L,          2L,                    3L,                       11L,                         5L,     1L,           2L,              25L,                          0L,                1L,          2L,                    3L,                       11L,                         5L,     1L,           2L,              25L,
  "KMT2D",                          3L,                0L,          0L,                    0L,                       11L,                         8L,     1L,           1L,              24L,                          2L,                0L,          0L,                    1L,                        6L,                         1L,     0L,           0L,              10L,
  "WT1",                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,          22L,              23L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               0L,
  "TET2",                          0L,                0L,          0L,                    1L,                        4L,                         2L,     0L,          11L,              17L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               0L,
  "PRPF40B",                          0L,                3L,          0L,                    2L,                        2L,                         2L,     1L,           2L,              12L,                          0L,                3L,          0L,                    2L,                        2L,                         2L,     1L,           2L,              12L,
  "NOTCH1",                          0L,                0L,          0L,                    3L,                        6L,                         0L,     0L,           0L,               9L,                          0L,                0L,          0L,                    1L,                        4L,                         0L,     0L,           0L,               5L,
  "PRPF8",                          0L,                0L,          1L,                    0L,                        6L,                         1L,     0L,           0L,               8L,                          0L,                0L,          1L,                    0L,                        6L,                         1L,     0L,           0L,               8L,
  "IDH2",                          1L,                0L,          1L,                    0L,                        6L,                         0L,     0L,           0L,               8L,                          1L,                0L,          1L,                    0L,                        6L,                         0L,     0L,           0L,               8L,
  "GNAS",                          0L,                1L,          0L,                    0L,                        2L,                         1L,     0L,           0L,               4L,                          0L,                1L,          0L,                    0L,                        2L,                         1L,     0L,           0L,               4L,
  "CSDE1",                          0L,                0L,          0L,                    1L,                        1L,                         0L,     0L,           0L,               2L,                          0L,                0L,          0L,                    1L,                        1L,                         0L,     0L,           0L,               2L,
  "KMT2A",                          2L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               2L,                          2L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               2L,
  "SETD2",                          0L,                0L,          0L,                    0L,                        1L,                         0L,     0L,           0L,               1L,                          0L,                0L,          0L,                    0L,                        1L,                         0L,     0L,           0L,               1L,
  "CREBBP",                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               1L,                          0L,                0L,          0L,                    1L,                        0L,                         0L,     0L,           0L,               1L,
  "SMC3",                          0L,                0L,          0L,                    1L,                        1L,                         0L,     0L,           0L,               1L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               0L,
  "STAT3",                          0L,                0L,          0L,                    0L,                        1L,                         0L,     0L,           0L,               1L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               0L,
  "CBL",                          0L,                0L,          0L,                    0L,                        0L,                         1L,     0L,           0L,               1L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               0L,
  "RUNX1",                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           1L,               1L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,     0L,           0L,               0L
)


dat <- dat %>% 
  `colnames<-`(str_replace_all(colnames(.), "\\.", " ")) 
dat1 <- dat %>% 
  `colnames<-`(str_c(colnames(.), c("",rep(" SM in M-CHIP genes", 9), rep(" M-CHIP", 9)))) %>% 
  select(-starts_with("All")) %>% 
  pivot_longer(cols = -Symbol) %>% 
  mutate(data = str_extract(name, "SM in M-CHIP genes|M-CHIP")) %>% 
  mutate(name = str_remove(name, " SM in M-CHIP genes| M-CHIP")) %>% 
  mutate(name = factor(name, levels = c(
    "African/African American",
    "Ashkenazi Jewish",
    "East Asian",
    "European (Finnish)",
    "European (non-Finnish)",
    "Latino/Admixed American",
    "South Asian",
    "Other"
  )))

dat1 %>% 
  mutate(data = factor(data, levels = c("SM in M-CHIP genes", "M-CHIP"))) %>% 
  group_by(data, Symbol) %>% 
  mutate(sum = case_when(
    data == "SM in M-CHIP genes"      ~ sum(value)
  )) %>% 
  arrange(data, desc(sum)) %>% 
  # ggplot(aes(x= fct_reorder(Symbol, value), y= value, fill = name))+
  ggplot(aes(x= fct_reorder(Symbol, sum), y= value, fill = name))+
  geom_bar(stat = "identity")+
  coord_flip()+
  facet_wrap(. ~ data)+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "italic"))

ggsave("Figure with table 3 ACs.pdf",
       width = 7,
       height = 4,
       dpi = 600)

# L-CHIP

dat <- tibble::tribble(
            ~Symbol, ~`African/African.American`, ~Ashkenazi.Jewish, ~East.Asian, ~`European.(Finnish)`, ~`European.(non-Finnish)`, ~`Latino/Admixed.American`, ~South.Asian, ~Other, ~All.Populations, ~`African/African.American`, ~Ashkenazi.Jewish, ~East.Asian, ~`European.(Finnish)`, ~`European.(non-Finnish)`, ~`Latino/Admixed.American`, ~South.Asian, ~Other, ~All.populations,
             "TSC2",                         21L,                0L,          0L,                    1L,                        1L,                         8L,           0L,     1L,              32L,                          0L,                0L,          0L,                    1L,                        1L,                         0L,           0L,     1L,               3L,
          "SMARCA4",                          4L,                1L,          3L,                    1L,                       15L,                         0L,           0L,     0L,              24L,                          4L,                1L,          3L,                    1L,                       14L,                         0L,           0L,     0L,              23L,
           "ARID5B",                          4L,                0L,          3L,                    2L,                        6L,                         0L,           0L,     0L,              15L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
            "KMT2D",                          3L,                0L,          0L,                    1L,                       11L,                         8L,           1L,     0L,              24L,                          2L,                0L,          0L,                    1L,                        6L,                         1L,           0L,     0L,              10L,
            "PTPRD",                          4L,                0L,          0L,                    0L,                        6L,                         5L,           0L,     0L,              15L,                          0L,                0L,          0L,                    0L,                        2L,                         0L,           0L,     0L,               2L,
           "NOTCH1",                          0L,                0L,          0L,                    3L,                        6L,                         1L,           1L,     0L,              11L,                          0L,                0L,          0L,                    1L,                        4L,                         0L,           1L,     0L,               6L,
           "ARID1A",                          0L,                0L,          0L,                    0L,                        3L,                         3L,           3L,     0L,               9L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
            "CHEK2",                          0L,                0L,          0L,                    0L,                        4L,                         1L,           2L,     1L,               8L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
            "KMT2C",                          0L,                0L,          1L,                    1L,                        4L,                         0L,           0L,     1L,               7L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
             "POT1",                          1L,                0L,          0L,                    0L,                        5L,                         0L,           0L,     0L,               6L,                          1L,                0L,          0L,                    0L,                        5L,                         0L,           0L,     0L,               6L,
           "ARID1B",                          0L,                1L,          0L,                    0L,                        2L,                         1L,           0L,     0L,               4L,                          0L,                1L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               1L,
             "DTX1",                          0L,                0L,          0L,                    0L,                        2L,                         0L,           1L,     0L,               3L,                          0L,                0L,          0L,                    0L,                        1L,                         0L,           1L,     0L,               2L,
             "KLF2",                          0L,                1L,          0L,                    0L,                        0L,                         0L,           2L,     0L,               3L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
           "CDKN2A",                          0L,                0L,          0L,                    0L,                        3L,                         2L,           0L,     0L,               5L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
             "PAX5",                          0L,                0L,          0L,                    0L,                        1L,                         1L,           1L,     0L,               3L,                          0L,                0L,          0L,                    0L,                        0L,                         1L,           0L,     0L,               1L,
             "SPEN",                          0L,                0L,          0L,                    0L,                        2L,                         0L,           0L,     0L,               2L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
             "BTG1",                          0L,                0L,          1L,                    0L,                        0L,                         0L,           1L,     0L,               2L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
              "MGA",                          0L,                0L,          0L,                    0L,                        2L,                         0L,           0L,     0L,               2L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
            "CIITA",                          0L,                0L,          0L,                    1L,                        1L,                         0L,           0L,     0L,               2L,                          0L,                0L,          0L,                    1L,                        1L,                         0L,           0L,     0L,               2L,
         "TNFRSF14",                          2L,                0L,          0L,                    0L,                        2L,                         0L,           1L,     0L,               5L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
            "STAT3",                          0L,                0L,          0L,                    0L,                        1L,                         0L,           0L,     0L,               1L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
            "FBXW7",                          0L,                0L,          0L,                    0L,                        0L,                         0L,           1L,     0L,               1L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
            "CCND3",                          0L,                0L,          0L,                    0L,                        0L,                         1L,           0L,     0L,               1L,                          0L,                0L,          0L,                    0L,                        0L,                         1L,           0L,     0L,               1L,
            "PRDM1",                          0L,                0L,          0L,                    0L,                        0L,                         0L,           1L,     0L,               1L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
          "TNFAIP3",                          0L,                0L,          0L,                    1L,                        0L,                         0L,           0L,     0L,               1L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,
             "FAT1",                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L,                          0L,                0L,          0L,                    0L,                        0L,                         0L,           0L,     0L,               0L
         )


dat <- dat %>% 
  `colnames<-`(str_replace_all(colnames(.), "\\.", " ")) 
dat1 <- dat %>% 
  `colnames<-`(str_c(colnames(.), c("",rep(" SM in L-CHIP genes", 9), rep(" L-CHIP", 9)))) %>% 
  select(-starts_with("All")) %>% 
  pivot_longer(cols = -Symbol) %>% 
  mutate(data = str_extract(name, "SM in L-CHIP genes|L-CHIP")) %>% 
  mutate(name = str_remove(name, " SM in L-CHIP genes| L-CHIP")) %>% 
  mutate(name = factor(name, levels = c(
    "African/African American",
    "Ashkenazi Jewish",
    "East Asian",
    "European (Finnish)",
    "European (non-Finnish)",
    "Latino/Admixed American",
    "South Asian",
    "Other"
  )))

dat1 %>% 
  mutate(data = factor(data, levels = c("SM in L-CHIP genes", "L-CHIP"))) %>% 
  group_by(data, Symbol) %>% 
  mutate(sum = case_when(
    data == "SM in L-CHIP genes"      ~ sum(value)
  )) %>% 
  arrange(data, desc(sum)) %>% 
  # ggplot(aes(x= fct_reorder(Symbol, value), y= value, fill = name))+
  ggplot(aes(x= fct_reorder(Symbol, sum), y= value, fill = name))+
  geom_bar(stat = "identity")+
  coord_flip()+
  facet_wrap(. ~ data)+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "italic"))

ggsave("Figure B with table 3 ACs.pdf",
       width = 7,
       height = 4,
       dpi = 600)





# Figure 4 C bins----
library(tidyverse)
library(ggridges)
theme_set(theme_classic(base_size = 10))

CH_ <- 
  read.delim(paste0(here::here(), "/PAFs of CH variants with zero values.txt"), sep = " ")

fig4_F <-
  CH_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, 
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  # filter(nc_freq_allele_count_nfe > 0 & nc_freq_allele_count_afr > 0 &
  #          nc_freq_allele_count_eas > 0 & nc_freq_allele_count_amr > 0 &
  #          nc_freq_allele_count_sas > 0 & nc_freq_allele_count_asj > 0 &
  #          nc_freq_allele_count_oth > 0 & nc_freq_allele_count_fin > 0) %>% 
  filter(if_all(where(is.numeric), \(x)(`>`(x, 0)))) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = case_when(
           name == "nfe"           ~ "European (non-Finnish)",
           name == "afr"           ~ "African/African American",
           name == "amr"           ~ "Latino/Admixed American",
           name == "eas"           ~ "East Asian",
           name == "sas"           ~ "South Asian",
           name == "asj"           ~ "Ashkenazi Jewish",
           name == "oth"           ~ "Other",
           name == "fin"           ~ "European (Finnish)",
         ), name = factor(name, levels = c("African/African American",
                                           "Ashkenazi Jewish",
                                           "East Asian",
                                           "European (Finnish)",
                                           "European (non-Finnish)",
                                           "Latino/Admixed American",
                                           "South Asian",
                                           "Other"))) %>% 
  filter(!is.na(value)) %>% 
  select(name, value) %>% 
  ggplot(aes(x=value, y=name, fill=name, color=name)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"
                                 ))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"
                                  ))+
  scale_x_continuous(limits = c(-0.00001,0.00015), 
                     labels = function(x) format(x, scientific = TRUE))+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")


  

legend <- ggpubr::get_legend(fig4_F)
# Convert to a ggplot and print
ggpubr::as_ggplot(legend)
ggsave("Figure4_legend_bottom.pdf",
       width = 11,
       height = 5,
       dpi = 600)

fig4_F+
  theme(legend.position = "none"
  )

ggsave("Figure 4 PAFs common CH variants.pdf",
       width = 5,
       height = 5,
       dpi = 600)


SM_ <- 
  read.delim(paste0(here::here(), "/PAFs of SM variants with zero values.txt"), sep = " ")
SM_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, 
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  filter(if_all(where(is.numeric), \(x)(`>`(x, 0)))) %>% 
  
  pivot_longer(cols = -c(IDs, SYMBOL)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = case_when(
           name == "nfe"           ~ "European (non-Finnish)",
           name == "afr"           ~ "African/African American",
           name == "amr"           ~ "Latino/Admixed American",
           name == "eas"           ~ "East Asian",
           name == "sas"           ~ "South Asian",
           name == "asj"           ~ "Ashkenazi Jewish",
           name == "oth"           ~ "Other",
           name == "fin"           ~ "European (Finnish)",
         ), name = factor(name, levels = c("Other",
                                           "South Asian",
                                           "Latino/Admixed American",
                                           "European (non-Finnish)",
                                           "European (Finnish)",
                                           "East Asian",
                                           "Ashkenazi Jewish",
                                           "African/African American"))) %>% 
  filter(!is.na(value)) %>% 
  select(name, value) %>% 
  # mutate(mutations_in_BickWHO = factor(mutations_in_BickWHO, levels= c("Yes", "No"))) %>% 
  
  ggplot(aes(x=value, y=name, fill=name, color=name)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins = 100, scale = 0.95, draw_baseline = FALSE)+
  # ggtitle("Popultaion Allele Frequency Distribution in CH Variants")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  scale_x_continuous(limits = c(-0.00001,0.00015), 
                     labels = function(x) format(x, scientific = TRUE))+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("Figure 4 PAFs SM common variants.pdf",
       width = 5,
       height = 5,
       dpi = 600)

# Fig4 new C - ALL Variant present all pop + scaling----
library(tidyverse)
library(ggridges)
theme_set(theme_classic(base_size = 15))
SM_ <- 
  read.delim(paste0(here::here(), "/PAFs of SM variants with zero values.txt"), sep = " ")

# 1.SM----
scaled_data <-
  SM_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, 
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  # filter to variants present in 2 pop
  # mutate(filter_2pop = rowSums(select(.,nc_freq_allele_count_afr:nc_freq_allele_count_sas) > 0, na.rm = TRUE)) %>% 
  # filter(filter_2pop >= 2) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_")) %>% 
  mutate(sample_size = case_when(
    name == "nfe"           ~ 51377,
    name == "afr"           ~ 7451,
    name == "amr"           ~ 17130,
    name == "eas"           ~ 8846,
    name == "sas"           ~ 15263,
    name == "asj"           ~ 4786,
    name == "oth"           ~ 2810,
    name == "fin"           ~ 10816,
  )) %>% 
  mutate(
    name = case_when(
      name == "nfe"           ~ "European (non-Finnish)",
      name == "afr"           ~ "African/African American",
      name == "amr"           ~ "Latino/Admixed American",
      name == "eas"           ~ "East Asian",
      name == "sas"           ~ "South Asian",
      name == "asj"           ~ "Ashkenazi Jewish",
      name == "oth"           ~ "Other",
      name == "fin"           ~ "European (Finnish)",
    ), name = factor(name, levels = c("Other",
                                      "South Asian",
                                      "Latino/Admixed American",
                                      "European (non-Finnish)",
                                      "European (Finnish)",
                                      "East Asian",
                                      "Ashkenazi Jewish",
                                      "African/African American"
    ))) %>% 
  filter(!is.na(value)) %>% 
  # re-scale
  mutate(value1 = value * sample_size) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value_scaled)) %>% 
  ungroup()

legend <- scaled_data %>% 
  select(name, value_scaled) %>% 
  filter(value_scaled > 0) %>% 
  ggplot(aes(x=value_scaled, y=name, fill=name, color=name, height = after_stat(density))) +
  stat_density_ridges(alpha=0.3, stat="density", scale = 1.3, 
                      bandwidth = 2e-05)+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  scale_x_continuous(limits = c(0,0.00034),
                     breaks = c(6.55e-05, 0.00013, 0.00019, 0.00026, 0.00032)
  )+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.title = element_blank(),
        legend.position = "bottom")

legend <- ggpubr::get_legend(legend)
# Convert to a ggplot and print
ggpubr::as_ggplot(legend)
ggsave("Figure4_legend_bottom.pdf",
       width = 11,
       height = 5,
       dpi = 600)


scaled_data %>% 
  select(name, value_scaled) %>% 
  filter(value_scaled > 0) %>% 
  ggplot(aes(x=value_scaled, y=name, fill=name, color=name, height = after_stat(density))) +
  stat_density_ridges(alpha=0.3, stat="density", scale = 1.2, 
                      bandwidth = 2e-05,
                      quantile_lines = TRUE, quantiles = 2, # add meadian
                      # quantile_lines = TRUE, quantile_fun=function(value_scaled,...)mean(value_scaled), # add mean
                      vline_color = "black", vline_width = 0.25)+
  # geom_vline(xintercept = 6.55e-05, linetype = 2)+
  # geom_vline(xintercept = 0.00013, linetype = 2)+
  # geom_vline(xintercept = 0.00019, linetype = 2)+
  # geom_vline(xintercept = 0.00026, linetype = 2)+
  # geom_vline(xintercept = 0.00032, linetype = 2)+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  scale_x_continuous(limits = c(0,0.00034),#c(2e-05,0.00013), # c(-0.00001,0.00015),
                     breaks = c(6.55e-05, 0.00013, 0.00019, 0.00026, 0.00032)
                     # labels = function(x) format(x, scientific = TRUE)
  )+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ggsave("Figure 4 PAFs SM common variants_02042025.pdf",
       width = 5,
       height = 5,
       dpi = 600)

# axis break
library(ggbreak)
scaled_data %>% 
  select(name, value_scaled) %>% 
  filter(value_scaled > 0) %>% 
  ggplot(aes(x=value_scaled, y=name, fill=name, color=name, height = after_stat(density))) +
  stat_density_ridges(alpha=0.3, stat="density", scale = 1.2, 
                      bandwidth = 2e-05,
                      quantile_lines = TRUE, quantiles = 2, # add meadian
                      # quantile_lines = TRUE, quantile_fun=function(value_scaled,...)mean(value_scaled), # add mean
                      vline_color = "black", vline_width = 0.25)+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  # scale_x_continuous(limits = c(0.00034,7)
  # )+
  expand_limits(y = c(0, 6.732228))+
  scale_x_break(c(0.02, 6.731888),# scales = 6.732228
                ticklabels=c(6.731888)
  )+
  # scale_x_continuous(limits = c(0,0.00034),#c(2e-05,0.00013), # c(-0.00001,0.00015),
  #                    breaks = c(6.55e-05, 0.00013, 0.00019, 0.00026, 0.00032)
                     # labels = function(x) format(x, scientific = TRUE)
  # )+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

# scaled_data %>% 
#   select(name, value_scaled, dat) %>% 
#   filter(value_scaled > 0) %>% 
#   ggplot(aes(x=value_scaled, y=dat, fill=name, color=name, height = after_stat(density))) +
#   stat_density_ridges(fill = "transparent", #alpha=0.2, 
#                       # stat="density", 
#                       scale = 1.8, bandwidth = 2e-05,
#                       quantile_lines = TRUE, quantiles = 2#, # add meadian
#                       # quantile_lines = TRUE, quantile_fun=function(value_scaled,...)mean(value_scaled), # add mean
#                       # vline_color = "black"
#   )+
#   # geom_vline(xintercept = 6.55e-05, linetype = 2)+
#   # geom_vline(xintercept = 0.00013, linetype = 2)+
#   # geom_vline(xintercept = 0.00019, linetype = 2)+
#   # geom_vline(xintercept = 0.00026, linetype = 2)+
#   # geom_vline(xintercept = 0.00032, linetype = 2)+
#   labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
#   # scale_fill_discrete(limits = c("European (non-Finnish)",
#   #                                "Latino/Admixed American",
#   #                                "South Asian",
#   #                                "European (Finnish)",
#   #                                "East Asian",
#   #                                "African/African American",
#   #                                "Ashkenazi Jewish"))+#,
#   #                                # "Other"))+
#   scale_color_discrete(limits = c("European (non-Finnish)",
#                                   "Latino/Admixed American",
#                                   "South Asian",
#                                   "European (Finnish)",
#                                   "East Asian",
#                                   "African/African American",
#                                   "Ashkenazi Jewish"))+#,
#   # "Other"))+
#   scale_x_continuous(limits = c(0,0.00034),#c(2e-05,0.00013), # c(-0.00001,0.00015),
#                      breaks = c(6.55e-05, 0.00013, 0.00019, 0.00026, 0.00032)
#                      # labels = function(x) format(x, scientific = TRUE)
#   )+
#   scale_y_discrete(expand = c(0, 0.2))+
#   theme(legend.title = element_blank(),
#         # axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         # legend.position = "none"
#   )



CH_ <- 
  read.delim(paste0(here::here(), "/PAFs of CH variants with zero values.txt"), sep = " ")

# 2.CH----
scaled_data <-
  CH_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, 
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  # mutate(filter_2pop = rowSums(select(.,nc_freq_allele_count_afr:nc_freq_allele_count_sas) > 0, na.rm = TRUE)) %>% 
  # filter(filter_2pop >= 2) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_")) %>% 
  mutate(sample_size = case_when(
    name == "nfe"           ~ 51377,
    name == "afr"           ~ 7451,
    name == "amr"           ~ 17130,
    name == "eas"           ~ 8846,
    name == "sas"           ~ 15263,
    name == "asj"           ~ 4786,
    name == "oth"           ~ 2810,
    name == "fin"           ~ 10816,
  )) %>% 
  mutate(
    name = case_when(
      name == "nfe"           ~ "European (non-Finnish)",
      name == "afr"           ~ "African/African American",
      name == "amr"           ~ "Latino/Admixed American",
      name == "eas"           ~ "East Asian",
      name == "sas"           ~ "South Asian",
      name == "asj"           ~ "Ashkenazi Jewish",
      name == "oth"           ~ "Other",
      name == "fin"           ~ "European (Finnish)",
    ), name = factor(name, levels = c("Other",
                                      "South Asian",
                                      "Latino/Admixed American",
                                      "European (non-Finnish)",
                                      "European (Finnish)",
                                      "East Asian",
                                      "Ashkenazi Jewish",
                                      "African/African American"
    ))) %>% 
  filter(!is.na(value)) %>% 
  # re-scale
  mutate(value1 = value * sample_size) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value_scaled)) %>% 
  ungroup()


scaled_data %>% 
  select(name, value_scaled) %>% 
  filter(value_scaled > 0) %>% 
  ggplot(aes(x=value_scaled, y=name, fill=name, color=name, height = after_stat(density))) +
  stat_density_ridges(alpha=0.3, stat="density", scale = 1.3, 
                      bandwidth = 2e-05,
                      quantile_lines = TRUE, quantiles = 2, # add meadian
                      # quantile_lines = TRUE, quantile_fun=function(value_scaled,...)mean(value_scaled), # add mean
                      vline_color = "black", vline_width = 0.25)+
  # geom_vline(xintercept = 6.55e-05, linetype = 2)+
  # geom_vline(xintercept = 0.00013, linetype = 2)+
  # geom_vline(xintercept = 0.00019, linetype = 2)+
  # geom_vline(xintercept = 0.00026, linetype = 2)+
  # geom_vline(xintercept = 0.00032, linetype = 2)+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  scale_x_continuous(limits = c(0,0.00034),#c(2e-05,0.00013), # c(-0.00001,0.00015),
                     breaks = c(6.55e-05, 0.00013, 0.00019, 0.00026, 0.00032)
                     # labels = function(x) format(x, scientific = TRUE)
  )+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ggsave("Figure 4 PAFs CH common variants_02042025.pdf",
       width = 5,
       height = 5,
       dpi = 600)


























# Table 3 A PAF CH by race ----
library(tidyverse)
library(gtsummary)
CH_ <- 
  read.delim(paste0(here::here(), "/PAFs of CH variants with zero values.txt"), sep = " ")

scaled_ch <-
  CH_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP, 
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  # mutate(filter_2pop = rowSums(select(.,nc_freq_allele_count_afr:nc_freq_allele_count_sas) > 0, na.rm = TRUE)) %>% 
  # filter(filter_2pop >= 2) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_")) %>% 
  mutate(sample_size = case_when(
    name == "nfe"           ~ 51377,
    name == "afr"           ~ 7451,
    name == "amr"           ~ 17130,
    name == "eas"           ~ 8846,
    name == "sas"           ~ 15263,
    name == "asj"           ~ 4786,
    name == "oth"           ~ 2810,
    name == "fin"           ~ 10816,
  )) %>% 
  mutate(
    name = case_when(
      name == "nfe"           ~ "European (non-Finnish)",
      name == "afr"           ~ "African/African American",
      name == "amr"           ~ "Latino/Admixed American",
      name == "eas"           ~ "East Asian",
      name == "sas"           ~ "South Asian",
      name == "asj"           ~ "Ashkenazi Jewish",
      name == "oth"           ~ "Other",
      name == "fin"           ~ "European (Finnish)",
    ), name = factor(name, levels = c("African/African American",
                                      "Ashkenazi Jewish",
                                      "East Asian",
                                      "European (Finnish)",
                                      "European (non-Finnish)",
                                      "Latino/Admixed American",
                                      "South Asian",
                                      "Other"
    ))) %>% 
  filter(!is.na(value)) %>% 
  # re-scale
  mutate(value1 = value * sample_size) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value_scaled)) %>% 
  ungroup()

CH_ %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  # select(IDs, nc_freq_allele_count) %>% 
  mutate(value1 = nc_freq_allele_count * 118479) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  # mutate(value_scaled = (nc_freq_allele_count) * 2) %>% 
  
  # filter(mutations_in_BickWHO == "Yes") %>% 
  # filter(L_CHIP == "Yes") %>% 
  
  filter(value_scaled > 0) %>% 
  select(value_scaled) %>% 
  tbl_summary(
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  )

paf_table_ch <- scaled_ch %>% 
  filter(value_scaled > 0)
  
paf_table_ch %>% 
  # select(IDs, nc_freq_allele_count_afr,
  #        nc_freq_allele_count_amr,
  #        nc_freq_allele_count_eas,
  #        nc_freq_allele_count_sas,
  #        nc_freq_allele_count_nfe) %>%
  # pivot_longer(cols = -IDs) %>%
  # mutate(value = as.numeric(value)) %>%
  # mutate(name = str_remove(name, "nc_freq_allele_count_"),
  #        name = factor(name, levels = c(
  #          "nfe_male", "nfe_female", "nfe",
  #          "eas_male", "eas_female", "eas", 
  #          "sas_male", "sas_female", "sas", 
  #          "amr_male", "amr_female", "amr", 
  #          "afr_male", "afr_female", "afr")
  #        )) %>% 
  # mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
  #   is.na(Sex)             ~ "overall",
  #   TRUE                   ~ Sex
  # )) %>%
  # mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
  #   is.na(Ethnicity)       ~ "overall",
  #   TRUE                   ~ Ethnicity
  # )) %>%
  # select(name, value,
  #        Ethnicity, Sex) %>% 
  # mutate(Ethnicity = case_when(
  #   Ethnicity == "nfe"           ~ "White",
  #   Ethnicity == "afr"           ~ "Black",
  #   Ethnicity == "amr"           ~ "Hispanic",
  #   Ethnicity == "eas"           ~ "East Asian",
  #   Ethnicity == "sas"           ~ "South Asian"
  # ), Ethnicity = factor(Ethnicity, levels = c("White",
  #                                             "South Asian",
  #                                             "Hispanic",
  #                                             "East Asian",
  #                                             "Black"))) %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**CH Pop AF race/eth include all variants population without zero in each pop - median**"
  ) %>% as_kable()

# M-CH
paf_table_ch <- scaled_ch %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  filter(value_scaled > 0)

paf_table_ch %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**M-CH Pop AF race/eth include all variants population without zero in each pop - median**"
  ) %>% as_kable()

# L-CH
paf_table_ch <- scaled_ch %>% 
  filter(L_CHIP == "Yes") %>% 
  filter(value_scaled > 0)

paf_table_ch %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**L-CH Pop AF race/eth include all variants population without zero in each pop - median**"
  ) %>% as_kable()


# Table 3A CH by sex----
scaled_ch <-
  CH_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP, 
         nc_freq_allele_count_male, nc_freq_allele_count_female
  ) %>%
  # mutate(filter_2pop = rowSums(select(.,nc_freq_allele_count_afr:nc_freq_allele_count_sas) > 0, na.rm = TRUE)) %>% 
  # filter(filter_2pop >= 2) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_")) %>% 
  mutate(sample_size = case_when(
    name == "male"             ~ 64629,
    name == "female"           ~ 53850
  )) %>% 
  mutate(
    name = factor(name, levels = c("female",
                                   "male"
    ))) %>% 
  filter(!is.na(value)) %>% 
  # re-scale
  mutate(value1 = value * sample_size) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value_scaled)) %>% 
  ungroup()

paf_table_ch <- scaled_ch %>% 
  filter(value_scaled > 0)

paf_table_ch %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**CH Pop AF sex include all variants population without zero in each pop - median**"
  ) %>% as_kable()

# M-CH
paf_table_ch <- scaled_ch %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  filter(value_scaled > 0)

paf_table_ch %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**M-CH Pop AF sex include all variants population without zero in each pop - median**"
  ) %>% as_kable()

# L-CH
paf_table_ch <- scaled_ch %>% 
  filter(L_CHIP == "Yes") %>% 
  filter(value_scaled > 0)

paf_table_ch %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**L-CH Pop AF sex include all variants population without zero in each pop - median**"
  ) %>% as_kable()



# Table 3 A PAF SM by race ----
library(tidyverse)
library(gtsummary)
SM_ <- 
  read.delim(paste0(here::here(), "/PAFs of SM variants with zero values.txt"), sep = " ")

SM_ %>% 
  distinct(IDs, .keep_all = TRUE) %>% 
  # select(IDs, nc_freq_allele_count) %>% 
  mutate(value1 = nc_freq_allele_count * 118479) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  # mutate(value_scaled = (nc_freq_allele_count) * 2) %>% 
  
  # filter(mutations_in_BickWHO == "Yes") %>%
  # filter(L_CHIP == "Yes") %>%
  
  filter(value_scaled > 0) %>% 
  select(value_scaled) %>% 
  tbl_summary(
    digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  )

scaled_SM <-
  SM_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP, 
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  # mutate(filter_2pop = rowSums(select(.,nc_freq_allele_count_afr:nc_freq_allele_count_sas) > 0, na.rm = TRUE)) %>% 
  # filter(filter_2pop >= 2) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_")) %>% 
  mutate(sample_size = case_when(
    name == "nfe"           ~ 51377,
    name == "afr"           ~ 7451,
    name == "amr"           ~ 17130,
    name == "eas"           ~ 8846,
    name == "sas"           ~ 15263,
    name == "asj"           ~ 4786,
    name == "oth"           ~ 2810,
    name == "fin"           ~ 10816,
  )) %>% 
  mutate(
    name = case_when(
      name == "nfe"           ~ "European (non-Finnish)",
      name == "afr"           ~ "African/African American",
      name == "amr"           ~ "Latino/Admixed American",
      name == "eas"           ~ "East Asian",
      name == "sas"           ~ "South Asian",
      name == "asj"           ~ "Ashkenazi Jewish",
      name == "oth"           ~ "Other",
      name == "fin"           ~ "European (Finnish)",
    ), name = factor(name, levels = c("African/African American",
                                      "Ashkenazi Jewish",
                                      "East Asian",
                                      "European (Finnish)",
                                      "European (non-Finnish)",
                                      "Latino/Admixed American",
                                      "South Asian",
                                      "Other"
    ))) %>% 
  filter(!is.na(value)) %>% 
  # re-scale
  mutate(value1 = value * sample_size) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value_scaled)) %>% 
  ungroup()

paf_table_SM <- scaled_SM %>% 
  filter(value_scaled > 0)

paf_table_SM %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM Pop AF race/eth include all variants population without zero in each pop - median**"
  ) %>% as_kable()

# M-CH
paf_table_SM <- scaled_SM %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  filter(value_scaled > 0)

paf_table_SM %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM in M-CH Pop AF race/eth include all variants population without zero in each pop - median**"
  ) %>% as_kable()

# L-CH
paf_table_SM <- scaled_SM %>% 
  filter(L_CHIP == "Yes") %>% 
  filter(value_scaled > 0)

paf_table_SM %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM in L-CH Pop AF race/eth include all variants population without zero in each pop - median**"
  ) %>% as_kable()


# Table 3A SM by sex----
scaled_SM <-
  SM_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP, 
         nc_freq_allele_count_male, nc_freq_allele_count_female
  ) %>%
  # mutate(filter_2pop = rowSums(select(.,nc_freq_allele_count_afr:nc_freq_allele_count_sas) > 0, na.rm = TRUE)) %>% 
  # filter(filter_2pop >= 2) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_")) %>% 
  mutate(sample_size = case_when(
    name == "male"             ~ 64629,
    name == "female"           ~ 53850
  )) %>% 
  mutate(
    name = factor(name, levels = c("female",
                                   "male"
    ))) %>% 
  filter(!is.na(value)) %>% 
  # re-scale
  mutate(value1 = value * sample_size) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value_scaled)) %>% 
  ungroup()

paf_table_SM <- scaled_SM %>% 
  filter(value_scaled > 0)

paf_table_SM %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM Pop AF sex include all variants population without zero in each pop - median**"
  ) %>% as_kable()

# M-CH
paf_table_SM <- scaled_SM %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  filter(value_scaled > 0)

paf_table_SM %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM in M-CH Pop AF sex include all variants population without zero in each pop - median**"
  ) %>% as_kable()

# L-CH
paf_table_SM <- scaled_SM %>% 
  filter(L_CHIP == "Yes") %>% 
  filter(value_scaled > 0)

paf_table_SM %>% 
  select(`Population Allele Frequency` = value_scaled, name) %>% 
  tbl_summary(by= name,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM in L-CH Pop AF sex include all variants population without zero in each pop - median**"
  ) %>% as_kable()










# FigS2 -Pie chart for nucleotide types in CH----
nucleotide_data <- 
  read_rds(paste0(here::here(), "/CH_ref_alt.rds"))
nucleotide_plot <- nucleotide_data %>% 
  mutate(alt_num = str_length(REF), .before = 2) %>% 
  mutate(ref_num = str_length(ALT), .before = 2) %>% 
  mutate(nuclotide_type = case_when(
    ref_num == 1 &
      alt_num == 1           ~ "single nucleotide",
    ref_num == 1 &
      alt_num < 3            ~ "small insertion",
    ref_num == 1 &
      alt_num >= 3            ~ "large insertion",
    ref_num < 3 &
      alt_num == 1            ~ "small deletion",
    ref_num >= 3 &
      alt_num == 1            ~ "large deletion",
    TRUE                      ~ str_c(REF, ALT, sep = ", ")
  ))
table(nucleotide_plot$nuclotide_type)

pdf("nucleotide type in CH.pdf")
nucleotide_plot %>% 
  group_by(nuclotide_type) %>% 
  mutate(ypos = n()) %>% 
  mutate(ypos = ypos - 0.5*ypos) %>% 
  ungroup() %>% 
  ggplot(aes(x= "", fill= nuclotide_type))+
  geom_bar(width=1, color= "white") +
  scale_fill_manual(values = c("#a6cee3", "#ff7f00",
                               "#993399",
                               "#1f78b4", "#fdbf6f"
                               ), name= "")+
  coord_polar("y", start=1)+
  theme_void()+
  theme(legend.position = "none")
dev.off()

legend <- nucleotide_plot %>% 
  group_by(nuclotide_type) %>% 
  mutate(ypos = n()) %>% 
  mutate(ypos = ypos - 0.5*ypos) %>% 
  ungroup() %>% 
  ggplot(aes(x= "", fill= nuclotide_type))+
  geom_bar(width=1, color= "white") +
  scale_fill_manual(values = c("#a6cee3", "#ff7f00",
                               "#993399",
                               "#1f78b4", "#fdbf6f"
  ), name= "")
legend <- ggpubr::get_legend(legend)
# Convert to a ggplot and print
pdf("nucleotide_legend.pdf")
ggpubr::as_ggplot(legend)
dev.off()

nucleotide_plot %>% 
  group_by(nuclotide_type) %>% 
  mutate(ypos = n()) %>% 
  distinct(nuclotide_type, ypos) %>% 
  ungroup() %>% 
  mutate(tot = sum(ypos)) %>% 
  mutate(perc = ypos / tot *100)
# A tibble: 5 × 4
# nuclotide_type     ypos   tot   perc
# <chr>             <int> <int>  <dbl>
#   1 single nucleotide 85083 89361 95.2  
# 2 small insertion    1453 89361  1.63 
# 3 small deletion      731 89361  0.818
# 4 large insertion    1322 89361  1.48 
# 5 large deletion      772 89361  0.864

# Pie chart for nucleotide types in gnomAD
nucleotide_data <- 
  read_rds(paste0(here::here(), "/clonal_mosaicism_ref_alt.rds"))
nucleotide_plot <- nucleotide_data %>% 
  mutate(alt_num = str_length(REF), .before = 2) %>% 
  mutate(ref_num = str_length(ALT), .before = 2) %>% 
  mutate(nuclotide_type = case_when(
    ref_num == 1 &
      alt_num == 1           ~ "single nucleotide",
    ref_num == 1 &
      alt_num < 3            ~ "small insertion",
    ref_num == 1 &
      alt_num >= 3            ~ "large insertion",
    ref_num < 3 &
      alt_num == 1            ~ "small deletion",
    ref_num >= 3 &
      alt_num == 1            ~ "large deletion",
    TRUE                      ~ str_c(REF, ALT, sep = ", ")
  ))
table(nucleotide_plot$nuclotide_type)

pdf("nucleotide type in SM.pdf")
nucleotide_plot %>% 
  group_by(nuclotide_type) %>% 
  mutate(ypos = n()) %>% 
  mutate(ypos = ypos - 0.5*ypos) %>% 
  ungroup() %>% 
  ggplot(aes(x= "", fill= nuclotide_type))+
  geom_bar(width=1, color= "white") +
  scale_fill_manual(values = c("#a6cee3", "#ff7f00",
                               "#993399",
                               "#1f78b4", "#fdbf6f"
  ), name= "")+
  coord_polar("y", start=1)+
  theme_void()+
  theme(legend.position = "none")
dev.off()

nucleotide_plot %>% 
  group_by(nuclotide_type) %>% 
  mutate(ypos = n()) %>% 
  distinct(nuclotide_type, ypos) %>% 
  ungroup() %>% 
  mutate(tot = sum(ypos)) %>% 
  mutate(perc = ypos / tot *100)
# A tibble: 5 × 4
# nuclotide_type      ypos    tot   perc
# <chr>              <int>  <int>  <dbl>
#   1 single nucleotide 473487 503703 94.0  
# 2 small deletion      3917 503703  0.778
# 3 small insertion     7945 503703  1.58 
# 4 large insertion    11408 503703  2.26 
# 5 large deletion      6946 503703  1.38 

# Pie chart for nucleotide types in gnomAD
nucleotide_data <- 
  read_rds(paste0(here::here(), "/gnomad_variants_ref_alt.rds"))
nucleotide_plot <- nucleotide_data %>% 
  mutate(alt_num = str_length(REF), .before = 2) %>% 
  mutate(ref_num = str_length(ALT), .before = 2) %>% 
  mutate(nuclotide_type = case_when(
    ref_num == 1 &
      alt_num == 1           ~ "single nucleotide",
    ref_num == 1 &
      alt_num < 3            ~ "small insertion",
    ref_num == 1 &
      alt_num >= 3            ~ "large insertion",
    ref_num < 3 &
      alt_num == 1            ~ "small deletion",
    ref_num >= 3 &
      alt_num == 1            ~ "large deletion",
    TRUE                      ~ str_c(REF, ALT, sep = ", ")
  ))
table(nucleotide_plot$nuclotide_type)

pdf("nucleotide type in gnomAD.pdf")
nucleotide_plot %>% 
  group_by(nuclotide_type) %>% 
  mutate(ypos = n()) %>% 
  mutate(ypos = ypos - 0.5*ypos) %>% 
  ungroup() %>% 
  ggplot(aes(x= "", fill= nuclotide_type))+
  geom_bar(width=1, color= "white") +
  scale_fill_manual(values = c("#a6cee3", "#ff7f00",
                               "#993399",
                               "#1f78b4", "#fdbf6f"
  ), name= "")+
  coord_polar("y", start=1)+
  theme_void()+
  theme(legend.position = "none")
dev.off()

nucleotide_plot %>% 
  group_by(nuclotide_type) %>% 
  mutate(ypos = n()) %>% 
  distinct(nuclotide_type, ypos) %>% 
  ungroup() %>% 
  mutate(tot = sum(ypos)) %>% 
  mutate(perc = ypos / tot *100)

# A tibble: 5 × 4
# nuclotide_type       ypos     tot  perc
# <chr>               <int>   <int> <dbl>
#   1 single nucleotide 3637547 3934541 92.5 
# 2 large deletion      57609 3934541  1.46
# 3 large insertion    109304 3934541  2.78
# 4 small insertion     74543 3934541  1.89
# 5 small deletion      55538 3934541  1.41


# SM
pop_freq <- 
  read.csv(paste0(here::here(), "/Pop freq data for CM.csv"))




# Table PAF by race SM ----
pop_freq <- 
  read.csv(paste0(here::here(), "/Pop freq data for CM.csv"))

paf_table_ch <- pop_freq %>% distinct(IDs, .keep_all = TRUE) %>% 
  filter(nc_freq_allele_count_nfe > 0 &
           nc_freq_allele_count_afr > 0 &
           nc_freq_allele_count_eas > 0 &
           nc_freq_allele_count_amr > 0 &
           nc_freq_allele_count_sas > 0)

paf_table_ch %>% 
  select(IDs, nc_freq_allele_count_afr,
         nc_freq_allele_count_amr,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_sas,
         nc_freq_allele_count_nfe) %>%
  pivot_longer(cols = -IDs) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_"),
         name = factor(name, levels = c(
           "nfe_male", "nfe_female", "nfe",
           "eas_male", "eas_female", "eas", 
           "sas_male", "sas_female", "sas", 
           "amr_male", "amr_female", "amr", 
           "afr_male", "afr_female", "afr")
         )) %>% 
  mutate(Sex = str_extract(name, "female|male"), Sex = case_when(
    is.na(Sex)             ~ "overall",
    TRUE                   ~ Sex
  )) %>%
  mutate(Ethnicity = str_extract(name, "afr|amr|nfe_nwe|nfe|eas|sas"), Ethnicity = case_when(
    is.na(Ethnicity)       ~ "overall",
    TRUE                   ~ Ethnicity
  )) %>%
  select(name, value,
         Ethnicity, Sex) %>% 
  mutate(Ethnicity = case_when(
    Ethnicity == "nfe"           ~ "White",
    Ethnicity == "afr"           ~ "Black",
    Ethnicity == "amr"           ~ "Hispanic",
    Ethnicity == "eas"           ~ "East Asian",
    Ethnicity == "sas"           ~ "South Asian"
  ), Ethnicity = factor(Ethnicity, levels = c("White",
                                              "South Asian",
                                              "Hispanic",
                                              "East Asian",
                                              "Black"))) %>% 
  select(`Population Allele Frequency` = value, Ethnicity) %>% 
  tbl_summary(by= Ethnicity,
              digits = all_continuous() ~ function(x) format(x, digits = 3, scientific = TRUE)
  ) %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**SM Pop AF race/eth include non-zero for all population - median**"
  ) %>% as_kable()







# VAF cutoff----

vaf_cutpoint <-
  readRDS("~/Documents/GitHub/Gillis/CH_gnomAD/population allele frequency data for VAF cut-point.rds")

# vaf_cutpoint %>%
#   filter(data == "Variants considered") %>%
#   ggplot(aes(x=nc_freq_allele_count)) +
#   geom_density(color = "#F8766D")+
#   labs(x= "Popultaion Allele Frequency", y= "Frequency Density")+
# 
#   scale_x_continuous(limits = c(-0.00001,0.0005),
#                      labels = function(x) format(x, scientific = TRUE))
# 
# vaf_cutpoint %>%
#   filter(data == "Clonal hematopoiesis (CH) variants") %>%
#   ggplot(aes(x=nc_freq_allele_count)) +
#   geom_density(color = "#00BFC4")+
#   labs(x= "Popultaion Allele Frequency", y= "Frequency Density")+
# 
#   scale_x_continuous(limits = c(-0.00001,0.0005),
#                      labels = function(x) format(x, scientific = TRUE))
# 
# 
# vaf_cutpoint %>%
#   filter(data == "Variants considered") %>%
#   ggplot(aes(x=nc_freq_allele_count)) +
#   geom_density(color = "#00BFC4")+
#   geom_density(data = vaf_cutpoint %>%
#                  filter(data == "Clonal hematopoiesis (CH) variants"),
#                aes(x=nc_freq_allele_count), color = "#F8766D")+
#   labs(x= "Popultaion Allele Frequency", y= "Frequency Density")+
# 
#   scale_x_continuous(limits = c(-0.00001,0.0005),
#                      labels = function(x) format(x, scientific = TRUE))

vaf_cutpoint %>%
  filter(data == "Variants considered") %>%
  ggplot(aes(x=nc_freq_allele_count)) +
  geom_density(fill = "#F8766D", alpha = 0.5)+
  geom_density(data = vaf_cutpoint %>%
                 filter(data == "Clonal hematopoiesis (CH) variants"),
               aes(x=nc_freq_allele_count), fill = "#00BFC4", alpha = 0.5)+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density")+

  scale_x_continuous(limits = c(-0.00001,0.0002),# labels = function(x) format(x, scientific = TRUE)
                     breaks = c(seq(0, 0.0002, by= 0.00005), 0.000025)
  )
  # scale_x_continuous(limits = c(-0.00001,0.0005),
  #                    labels = function(x) format(x, scientific = TRUE))



# vaf_cutpoint %>%
#   mutate(nc_freq_allele_count = round(nc_freq_allele_count, 2)) %>%
#   filter(!is.na(nc_freq_allele_count)) %>%
#   group_by(nc_freq_allele_count, data) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
# 
#   ggplot(aes(x= nc_freq_allele_count, y=count, fill= data
#   ))+
#   geom_area(position = "identity") +
#   scale_x_continuous(limits = c(0,0.02), labels = function(x) format(x, scientific = TRUE)
#   )+
#   labs(x= "Popultaion Allele Frequency", y= "Variant (#)")+
#   coord_cartesian(ylim = c(0, 90000))+
#   theme(legend.position = "bottom")
# 
# vaf_cutpoint %>%
#   mutate(nc_freq_allele_count = round(nc_freq_allele_count, 2)) %>%
#   filter(!is.na(nc_freq_allele_count)) %>%
#   group_by(nc_freq_allele_count, data) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
# 
#   ggplot(aes(x= nc_freq_allele_count, y=count, fill= data
#   ))+
#   geom_area(position = "identity") +
#   scale_x_continuous(limits = c(0,0.02), labels = function(x) format(x, scientific = TRUE)
#   )+
#   scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
#   labs(x= "Popultaion Allele Frequency", y= "Variant (#)")+
#   theme(legend.position = "bottom")+
#   facet_zoom(ylim = c(0, 90000))

ch_in_known_bickwho <- 
  read.csv(paste0(here::here(), "/known Bick+WHO variants in CH_variants.csv")) %>% 
  mutate(is_ch_variants = "Yes")

vaf_cutpoint <- vaf_cutpoint %>% 
  left_join(ch_in_known_bickwho %>% 
              select(IDs, is_ch_variants), by = "IDs")

vaf_cutpoint %>%
  mutate(is_ch_variants = case_when(
    is_ch_variants == "Yes"        ~ "M-CHIP in gnomAD",
    TRUE                           ~ "No M-CHIP in gnomAD"
  )) %>% 
  # filter(data == "Variants considered") %>%
  ggplot(aes(x=nc_freq_allele_count, fill = is_ch_variants)) +
  geom_density( alpha = 0.5)+
  # geom_density(data = vaf_cutpoint %>%
  #                filter(data == "Clonal hematopoiesis (CH) variants"),
  #              aes(x=nc_freq_allele_count), color = "#F8766D")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density", fill = "")+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  scale_x_continuous(limits = c(-0.00001,0.0002),# labels = function(x) format(x, scientific = TRUE)
                     breaks = c(seq(0, 0.0002, by= 0.00005), 0.000025)
  )
  scale_x_continuous(limits = c(-0.00001,0.0005),
                     labels = function(x) format(x, scientific = TRUE))

# vaf_cutpoint %>%
#   mutate(nc_freq_allele_count = round(nc_freq_allele_count, 2)) %>%
#   filter(!is.na(nc_freq_allele_count)) %>%
#   group_by(nc_freq_allele_count, data) %>%
#   summarise(count = n()) %>%
#   ungroup() %>% 
#   ggplot(aes(x= nc_freq_allele_count, fill= data
#   ))+
#   # ggplot(aes(x=nc_freq_allele_count)) +
#   geom_density(alpha = 0.5)+
#   labs(x= "Popultaion Allele Frequency", y= "Frequency Density")+
#   
#   scale_x_continuous(limits = c(-0.00001,1),# labels = function(x) format(x, scientific = TRUE)
#                      breaks = c(seq(0, 1, by= 0.05), 0.01, 0.02)
#   )+
#   # scale_y_discrete(expand = c(0, 0.2))+
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())+
#   facet_zoom(xlim = c(-0.00001, 0.2))

# Age

# age_data <- readRDS("~/Documents/GitHub/Gillis/CH_gnomAD/age.rds")
# 
# age_data <- age_data %>% filter(variant_in_cosmic == "Yes") %>% 
#   arrange(desc(data), Primary.site)
# 
# 
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>%
#   # filter(variant_in_cosmic == "Yes") %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   nrow()
# 
# age_data %>% 
#   filter(data == "Variants considered") %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   nrow()
# 
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Age) %>% 
#   tbl_summary()
# 
# # Table S2
# age_data %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Age, data) %>% 
#   tbl_summary(by = data,
#               sort = list(everything() ~ "frequency")) %>% 
#   add_overall() %>% 
#   add_p()
# 
# dat <- age_data %>% 
#   filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
#   distinct(Sample.name, .keep_all = TRUE) 
# 
# dat %>% 
#   select(Age, data) %>% 
#   tbl_summary(by = data,
#               sort = list(everything() ~ "frequency")) %>% 
#   add_overall() %>% 
#   add_p()
# 
# b <- age_data %>% 
#   arrange(desc(data)) %>% 
#   filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   filter(!Sample.name %in% c(dat$Sample.name))
# b %>% 
#   select(Age, data) %>% 
#   tbl_summary(by = data,
#               sort = list(everything() ~ "frequency")) %>% 
#   add_overall() %>% 
#   add_p()
# 
# # Table 4
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants" &
#            ) %>% 
#   arrange(Primary.site) %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Primary.site) %>% 
#   tbl_summary(sort = list(everything() ~ "frequency"))
# 
# # Table S3
# age_data %>% 
#   arrange(desc(data), Primary.site) %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Primary.site, data) %>% 
#   tbl_summary(by = data,
#               sort = list(everything() ~ "frequency"),
#               percent = "row") %>% 
#   add_overall() %>% as_kable()
# 
# # Table S4
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>% 
#   arrange(desc(mutations_in_BickWHO)) %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Age, mutations_in_BickWHO) %>% 
#   tbl_summary(by = mutations_in_BickWHO) %>% 
#   add_overall() %>% 
#   add_p()
# 
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>% 
#   filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
#   arrange(desc(mutations_in_BickWHO)) %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Age, mutations_in_BickWHO) %>% 
#   tbl_summary(by = mutations_in_BickWHO) %>% 
#   add_overall() %>% 
#   add_p()
# 
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>% 
#   filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
#   arrange(desc(mutations_in_BickWHO)) %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Age, mutations_in_BickWHO) %>% 
#   tbl_summary(by = mutations_in_BickWHO) %>% 
#   add_overall() %>% 
#   add_p()
# 
# 
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>% 
#   arrange(desc(L_CHIP)) %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Age, L_CHIP) %>% 
#   tbl_summary(by = L_CHIP) %>% 
#   add_overall() %>% 
#   add_p()
# 
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>% 
#   filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
#   arrange(desc(L_CHIP)) %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Age, L_CHIP) %>% 
#   tbl_summary(by = L_CHIP) %>% 
#   add_overall() %>% 
#   add_p()
# 
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>% 
#   filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
#   arrange(desc(L_CHIP)) %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Age, L_CHIP) %>% 
#   tbl_summary(by = L_CHIP) %>% 
#   add_overall() %>% 
#   add_p()


# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>% 
#   filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Primary.site) %>% 
#   tbl_summary(sort = list(everything() ~ "frequency"))
# 
# age_data %>% 
#   filter(data == "Clonal hematopoiesis (CH) variants") %>% 
#   filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
#   arrange(desc(L_CHIP)) %>% 
#   distinct(Sample.name, .keep_all = TRUE) %>% 
#   select(Sample.name, L_CHIP) %>% 
#   tbl_summary(sort = list(everything() ~ "frequency"))


# Table S4----
cosmic_primary <- 
  read_rds(paste0(here::here(), "/New_Fig5B_percent_Prevalence_Cancers_inCosmic.rds"))

a <- cosmic_primary %>% 
  filter(data == "CH") %>% 
  filter(sum_individuals_by_primary > 500) %>% 
  mutate(Primary.histology = coalesce(Primary.histology, Primary.site)) %>% 
  group_by(Primary.site) %>%
  mutate(perc_sum = sum(perc)) %>% 
  arrange(desc(perc_sum))

# Table 4 mutations ----
library(tidyverse)
library(gtsummary)

CH_mut <-
  read.delim(paste0(here::here(), "/COSMIC all cancers and mutations info in CH variants.txt"), sep = " ")
CH_mut2 <- 
  read.delim(paste0(here::here(), "/new COSMIC all cancers and mutations info in CH variants.txt"), sep = " ")

a <- CH_mut2 %>% distinct(Sample.name, Primary.site, Primary.histology)

CH_mut2 %>% 
  arrange(Sample.name, Primary.site) %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  select(Primary.site) %>% 
  tbl_summary(
    sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% 
  modify_header(
    label = "**Cosmic Patients/Cell Line Characteristics**"
  ) %>% 
  as_kable()
  
CH_mut2 %>% 
  arrange(Sample.name, Primary.site) %>% 
  distinct(Sample.name, .keep_all = TRUE) %>%
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  select(Primary.histology) %>% 
  tbl_summary(
    sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% 
  modify_header(
    label = "**Cosmic Patients/Cell Line Characteristics**"
  ) %>% 
  as_kable()

CH_mut2 %>% 
  select(Mutation.Description,
         FATHMM.prediction, Mutation.somatic.status) %>% 
  mutate(Mutation.somatic.status = case_when(
    str_detect(Mutation.somatic.status, "somatic")        ~ "reported or confirmed somatic",
    TRUE                                                  ~ Mutation.somatic.status
  )) %>% 
  mutate(Mutation.Description = case_when(
    Mutation.Description == ""           ~ NA_character_,
    Mutation.Description == "Unknown"    ~ NA_character_,
    TRUE                                 ~ Mutation.Description
  )) %>% 
  mutate(FATHMM.prediction = case_when(
    FATHMM.prediction == ""              ~ NA_character_,
    FATHMM.prediction == "Unknown"       ~ NA_character_,
    TRUE                                 ~ FATHMM.prediction
  )) %>% 
  tbl_summary(
    sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% 
  modify_header(
    label = "**Cosmic Patients/Cell Line Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Samples, careful if 1 sample shows multiple gene-mutations**") %>% 
  as_kable()

# Table S6 ----
library(tidyverse)
library(gtsummary)
cosmic_inGnom <- readRDS("~/Documents/GitHub/Gillis/CH_gnomAD/gnomad_variants.rds")
CH_mut2 <- 
  read.delim(paste0(here::here(), "/new COSMIC all cancers and mutations info in CH variants.txt"), sep = " ")

unique_CH_in_cosmic <- CH_mut2 %>% 
  arrange(Sample.name, Primary.site) %>% 
  distinct(Sample.name, .keep_all = TRUE)

bind_rows(unique_CH_in_cosmic %>% 
            mutate(data = "CH"),
          cosmic_inGnom %>% 
            arrange(Sample.name, Primary.site) %>% 
            filter(!Sample.name %in% c(unique_CH_in_cosmic$Sample.name)) %>%
            distinct(Sample.name, .keep_all = TRUE) %>% 
            mutate(data = "No CH")) %>% 
  distinct(Sample.name, .keep_all = TRUE) %>% 
  
  select(Age, data) %>%
  tbl_summary(by = data,
    sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% 
  add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients/Cell Line Characteristics**"
  ) %>% 
  as_kable()

bind_rows(unique_CH_in_cosmic %>% 
            filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
            mutate(data = "CH"),
          cosmic_inGnom %>% 
            arrange(Sample.name, Primary.site) %>% 
            filter(!Sample.name %in% c(unique_CH_in_cosmic$Sample.name)) %>%
            distinct(Sample.name, .keep_all = TRUE) %>% 
            filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
            mutate(data = "No CH")) %>% 
  distinct(Sample.name, .keep_all = TRUE) %>% 
  
  select(Age, data) %>%
  tbl_summary(by = data,
              sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% 
  add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients/Cell Line Characteristics**"
  ) %>% 
  as_kable()

bind_rows(unique_CH_in_cosmic %>% 
            filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
            mutate(data = "CH"),
          cosmic_inGnom %>% 
            arrange(Sample.name, Primary.site) %>% 
            filter(!Sample.name %in% c(unique_CH_in_cosmic$Sample.name)) %>%
            distinct(Sample.name, .keep_all = TRUE) %>% 
            filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
            mutate(data = "No CH")) %>% 
  distinct(Sample.name, .keep_all = TRUE) %>% 
  
  select(Age, data) %>%
  tbl_summary(by = data,
              sort = everything() ~ "frequency") %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_p() %>% 
  add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients/Cell Line Characteristics**"
  ) %>% 
  as_kable()

# Table S7----
library(tidyverse)
library(gtsummary)
cosmic_inGnom <- readRDS("~/Documents/GitHub/Gillis/CH_gnomAD/gnomad_variants.rds")

unique_CH_in_cosmic <- CH_mut2 %>% 
  arrange(Sample.name, Primary.site) %>% 
  distinct(Sample.name, .keep_all = TRUE)

bind_rows(unique_CH_in_cosmic %>% 
            mutate(data = "CH"),
          cosmic_inGnom %>% 
            arrange(Sample.name, Primary.site) %>% 
            distinct(Sample.name, .keep_all = TRUE) %>% 
            mutate(data = "No CH")) %>% 
  distinct(Sample.name, .keep_all = TRUE) %>% 
  
  select(Primary.site, data) %>%
  tbl_summary(by=data,
              sort = everything() ~ "frequency",
              percent = "row") %>%
  bold_labels() %>% add_stat_label() %>% add_overall() %>%
  add_p() %>% bold_p(t=.05) %>%
  modify_header(
    label = "**Prevalence of CH by cancer row percent**"
  ) %>%
  as_kable()


# Table S8----
CH_mut2 %>% distinct(IDs) %>% nrow()

CH_mut2 %>% 
  arrange(Sample.name, Primary.site, desc(mutations_in_BickWHO)) %>% 
  # arrange(desc(mutations_in_BickWHO)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In Bick+WHO**") %>% 
  modify_header(
    label = "**Age of COSMIC patients with CH variants in Bick & WHO**"
  ) %>% 
  as_kable()

CH_mut2 %>% 
  arrange(Sample.name, Primary.site, desc(mutations_in_BickWHO)) %>% 
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(mutations_in_BickWHO)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In Bick+WHO**") %>% 
  modify_header(
    label = "**Age of COSMIC heme patients with CH variants in Bick & WHO**"
  ) %>% 
  as_kable()

CH_mut2 %>% 
  arrange(Sample.name, Primary.site, desc(mutations_in_BickWHO)) %>% 
  filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
  arrange(desc(mutations_in_BickWHO)) %>%
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, mutations_in_BickWHO) %>% 
  tbl_summary(by= mutations_in_BickWHO#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In Bick+WHO**") %>% 
  modify_header(
    label = "**Age of COSMIC non-heme patients with CH variants in Bick & WHO**"
  ) %>% 
  as_kable()

CH_mut2 %>% 
  arrange(Sample.name, Primary.site, desc(L_CHIP)) %>% 
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, L_CHIP) %>% 
  tbl_summary(by= L_CHIP#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In L_CHIP**") %>% 
  modify_header(
    label = "**Age of COSMIC patients with CH variants L_CHIP**"
  ) %>% 
  as_kable()

CH_mut2 %>% 
  arrange(Sample.name, Primary.site, desc(L_CHIP)) %>% 
  filter(Primary.site == "haematopoietic_and_lymphoid_tissue") %>% 
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, L_CHIP) %>% 
  tbl_summary(by= L_CHIP#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In L_CHIP**") %>% 
  modify_header(
    label = "**Age of COSMIC heme patients with CH variants L_CHIP**"
  ) %>% 
  as_kable()

CH_mut2 %>% 
  arrange(Sample.name, Primary.site, desc(L_CHIP)) %>% 
  filter(Primary.site != "haematopoietic_and_lymphoid_tissue") %>% 
  distinct(Sample.name, .keep_all = TRUE) %>% 
  select(Age, L_CHIP) %>% 
  tbl_summary(by= L_CHIP#,
              # statistic = list(all_continuous() ~ "{mean} ({sd})"),
              # digits = all_continuous() ~ 2
  ) %>% 
  bold_labels() %>% add_stat_label() %>%
  add_p() %>% bold_p(t=.05) %>% add_overall() %>% 
  modify_header(
    label = "**Cosmic Patients Characteristics**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**CH In L_CHIP**") %>% 
  modify_header(
    label = "**Age of COSMIC non-heme patients with CH variants L_CHIP**"
  ) %>% 
  as_kable()





































# Figure S5----
library(tidyverse)
theme_set(theme_classic(base_size = 14))
SM_ <- 
  read.delim(paste0(here::here(), "/PAFs of SM variants with zero values.txt"), sep = " ")

# find order of variants by overall PAF
order <- SM_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  mutate(considered_variants = case_when(
    IDs == "2-25234373-C-T"                     ~ "DNMT3A R882H",
    IDs == "2-25234374-G-A"                     ~ "DNMT3A R882C",
    TRUE                                        ~ considered_variants
  )) %>% 
  filter(!is.na(considered_variants)) %>%
  filter(considered_variants != "DNMT3A W860R") %>% 
  select(considered_variants, nc_freq_allele_count) %>% 
  arrange(desc(nc_freq_allele_count))
order$considered_variants


scaled_data <-
  SM_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, considered_variants,
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  mutate(considered_variants = case_when(
    IDs == "2-25234373-C-T"                     ~ "DNMT3A R882H",
    IDs == "2-25234374-G-A"                     ~ "DNMT3A R882C",
    TRUE                                        ~ considered_variants
  )) %>% 
  # mutate(filter_2pop = rowSums(select(.,nc_freq_allele_count_afr:nc_freq_allele_count_sas) > 0, na.rm = TRUE)) %>% # We don't apply this here
  # filter(filter_2pop >= 2) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL, considered_variants)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_")) %>% 
  mutate(sample_size = case_when(
    name == "nfe"           ~ 51377,
    name == "afr"           ~ 7451,
    name == "amr"           ~ 17130,
    name == "eas"           ~ 8846,
    name == "sas"           ~ 15263,
    name == "asj"           ~ 4786,
    name == "oth"           ~ 2810,
    name == "fin"           ~ 10816,
  )) %>% 
  mutate(
    name = case_when(
      name == "nfe"           ~ "European (non-Finnish)",
      name == "afr"           ~ "African/African American",
      name == "amr"           ~ "Latino/Admixed American",
      name == "eas"           ~ "East Asian",
      name == "sas"           ~ "South Asian",
      name == "asj"           ~ "Ashkenazi Jewish",
      name == "oth"           ~ "Other",
      name == "fin"           ~ "European (Finnish)",
    ), name = factor(name, levels = c("Other",
                                      "South Asian",
                                      "Latino/Admixed American",
                                      "European (non-Finnish)",
                                      "European (Finnish)",
                                      "East Asian",
                                      "Ashkenazi Jewish",
                                      "African/African American"
    ))) %>% 
  filter(!is.na(value)) %>% 
  # re-scale
  mutate(value1 = value * sample_size) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value_scaled)) %>% 
  ungroup() %>% 
  mutate(considered_variants = factor(considered_variants, levels = c(order$considered_variants)))


scaled_data %>%
  filter(!is.na(considered_variants)) %>%
  filter(considered_variants != "DNMT3A W860R") %>%
  ggplot(aes(x=name, y= value_scaled, color= name, shape = value_scaled == 0))+
  geom_point(size= 3)+
  labs(#title = "Variant Allele Frequency in Non-Cancer population", 
    y= "Popultaion Allele Frequency", x= NULL)+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  scale_shape_manual(values=c(16, 1), labels = c("PAF > 0", "PAF = 0"), )+
  theme(legend.title=element_blank(),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1,
                                   hjust=1),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "italic"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))+
  facet_wrap(. ~ considered_variants, ncol = 2)+ 
  guides(color = guide_legend(nrow = 3),
         shape = guide_legend(nrow = 2))+
  coord_flip()
ggsave("Fig S5 PAF SM considered variants_02102025_shape.pdf",
       width = 6.5,
       height = 10,
       dpi = 600)


# Figure 6 B-F----
library(tidyverse)
library(ggridges)
theme_set(theme_classic(base_size = 15))
SM_ <- 
  read.delim(paste0(here::here(), "/PAFs of SM variants with zero values.txt"), sep = " ")

# 1.SM----
scaled_data <-
  SM_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP,
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  # mutate(filter_2pop = rowSums(select(.,nc_freq_allele_count_afr:nc_freq_allele_count_sas) > 0, na.rm = TRUE)) %>% 
  # filter(filter_2pop >= 2) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_")) %>% 
  mutate(sample_size = case_when(
    name == "nfe"           ~ 51377,
    name == "afr"           ~ 7451,
    name == "amr"           ~ 17130,
    name == "eas"           ~ 8846,
    name == "sas"           ~ 15263,
    name == "asj"           ~ 4786,
    name == "oth"           ~ 2810,
    name == "fin"           ~ 10816,
  )) %>% 
  mutate(
    name = case_when(
      name == "nfe"           ~ "European (non-Finnish)",
      name == "afr"           ~ "African/African American",
      name == "amr"           ~ "Latino/Admixed American",
      name == "eas"           ~ "East Asian",
      name == "sas"           ~ "South Asian",
      name == "asj"           ~ "Ashkenazi Jewish",
      name == "oth"           ~ "Other",
      name == "fin"           ~ "European (Finnish)",
    ), name = factor(name, levels = c("Other",
                                      "South Asian",
                                      "Latino/Admixed American",
                                      "European (non-Finnish)",
                                      "European (Finnish)",
                                      "East Asian",
                                      "Ashkenazi Jewish",
                                      "African/African American"
    ))) %>% 
  filter(!is.na(value)) %>% 
  # re-scale
  mutate(value1 = value * sample_size) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value_scaled)) %>% 
  ungroup()

scaled_data %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  select(name, value_scaled) %>% 
  filter(value_scaled > 0) %>% 
  ggplot(aes(x=value_scaled, y=name, fill=name, color=name, height = after_stat(density))) +
  stat_density_ridges(alpha=0.3, stat="density", scale = 1.3, 
                      bandwidth = 2e-05,
                      quantile_lines = TRUE, quantiles = 2, # add meadian
                      vline_color = "black", vline_width = 0.25)+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  scale_x_continuous(limits = c(0,0.00034),#c(2e-05,0.00013), # c(-0.00001,0.00015),
                     breaks = c(6.55e-05, 0.00013, 0.00019, 0.00026, 0.00032)
                     # labels = function(x) format(x, scientific = TRUE)
  )+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ggsave("Figure 7 PAFs SM MCHIP_02042025.pdf",
       width = 5,
       height = 5,
       dpi = 600)

# In L-CHIP, there is pop without any variants. To avaoid blank spot in the plot
# I create a data with very low value to get a base line for these pop
dat <- tibble(name = c("Ashkenazi Jewish", "Other"), 
              # x = c(0,0,0), 
              value_scaled = c(0.0003269032, 0.0003269032))
scaled_data %>% 
  filter(L_CHIP == "Yes") %>% 
  select(name, value_scaled) %>% 
  filter(value_scaled > 0) %>% 
  ggplot(aes(x=value_scaled, y=name, fill=name, color=name, height = after_stat(density))) +
  stat_density_ridges(alpha=0.3, stat="density", scale = 1.3, 
                      bandwidth = 2e-05,
                      quantile_lines = TRUE, quantiles = 2, # add meadian
                      vline_color = "black", vline_width = 0.25)+
  # geom_hline(yintercept = 2, color = "#CD9600")+
  # geom_hline(yintercept = 3, color = "#7CAE00")+
  # geom_hline(yintercept = 8, color = "#FF61CC")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  scale_x_continuous(limits = c(0,0.00034),#c(2e-05,0.00013), # c(-0.00001,0.00015),
                     breaks = c(6.55e-05, 0.00013, 0.00019, 0.00026, 0.00032)
                     # labels = function(x) format(x, scientific = TRUE)
  )+
  scale_y_discrete(expand = c(0, 0.2))+
  geom_segment(data = dat, 
               aes(y = name, x = value_scaled-0.0003213779, xend = value_scaled), 
               show.legend = F, inherit.aes = F, color = c("#CD9600", "#FF61CC"))+
  # geom_segment(data = dat, 
  #              aes(xend = name, y = x, 
  #                  yend = xend), show.legend = F)
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ggsave("Figure 7 PAFs SM LCHIP_02042025.pdf",
       width = 5,
       height = 5,
       dpi = 600)


CH_ <- 
  read.delim(paste0(here::here(), "/PAFs of CH variants with zero values.txt"), sep = " ")

# 2.CH----
scaled_data <-
  CH_ %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  select(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP,
         nc_freq_allele_count_afr, nc_freq_allele_count_asj,
         nc_freq_allele_count_eas,
         nc_freq_allele_count_fin, nc_freq_allele_count_nfe,
         nc_freq_allele_count_amr, nc_freq_allele_count_oth,
         nc_freq_allele_count_sas
  ) %>%
  # mutate(filter_2pop = rowSums(select(.,nc_freq_allele_count_afr:nc_freq_allele_count_sas) > 0, na.rm = TRUE)) %>% 
  # filter(filter_2pop >= 2) %>% 
  pivot_longer(cols = -c(IDs, SYMBOL, mutations_in_BickWHO, L_CHIP)) %>%
  mutate(name = str_remove(name, "nc_freq_allele_count_")) %>% 
  mutate(sample_size = case_when(
    name == "nfe"           ~ 51377,
    name == "afr"           ~ 7451,
    name == "amr"           ~ 17130,
    name == "eas"           ~ 8846,
    name == "sas"           ~ 15263,
    name == "asj"           ~ 4786,
    name == "oth"           ~ 2810,
    name == "fin"           ~ 10816,
  )) %>% 
  mutate(
    name = case_when(
      name == "nfe"           ~ "European (non-Finnish)",
      name == "afr"           ~ "African/African American",
      name == "amr"           ~ "Latino/Admixed American",
      name == "eas"           ~ "East Asian",
      name == "sas"           ~ "South Asian",
      name == "asj"           ~ "Ashkenazi Jewish",
      name == "oth"           ~ "Other",
      name == "fin"           ~ "European (Finnish)",
    ), name = factor(name, levels = c("Other",
                                      "South Asian",
                                      "Latino/Admixed American",
                                      "European (non-Finnish)",
                                      "European (Finnish)",
                                      "East Asian",
                                      "Ashkenazi Jewish",
                                      "African/African American"
    ))) %>% 
  filter(!is.na(value)) %>% 
  # re-scale
  mutate(value1 = value * sample_size) %>% 
  mutate(value_scaled = (value1 / 15263) * 2) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value_scaled)) %>% 
  ungroup()

scaled_data %>% 
  filter(mutations_in_BickWHO == "Yes") %>% 
  select(name, value_scaled) %>% 
  filter(value_scaled > 0) %>% 
  ggplot(aes(x=value_scaled, y=name, fill=name, color=name, height = after_stat(density))) +
  stat_density_ridges(alpha=0.3, stat="density", scale = 1.3, 
                      bandwidth = 2e-05,
                      quantile_lines = TRUE, quantiles = 2, # add meadian
                      vline_color = "black", vline_width = 0.25)+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  scale_x_continuous(limits = c(0,0.00034),#c(2e-05,0.00013), # c(-0.00001,0.00015),
                     breaks = c(6.55e-05, 0.00013, 0.00019, 0.00026, 0.00032)
                     # labels = function(x) format(x, scientific = TRUE)
  )+
  scale_y_discrete(expand = c(0, 0.2))+
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ggsave("Figure 7 PAFs CH MCHIP_02042025.pdf",
       width = 5,
       height = 5,
       dpi = 600)

dat <- tibble(name = c("Ashkenazi Jewish","East Asian", 
                       "South Asian",
                       "Other"), 
              value_scaled = c(rep(0.0003269032, 4))) %>% 
  mutate(name = factor(name, levels = c("Other",
                                        "South Asian",
                                        "Latino/Admixed American",
                                        "European (non-Finnish)",
                                        "European (Finnish)",
                                        "East Asian",
                                        "Ashkenazi Jewish",
                                        "African/African American"
  )))

scaled_data1 <- scaled_data #%>% 
  # mutate(value_scaled = case_when(
  #   IDs == "16-2084567-TC-T" &
  #     name == "Ashkenazi Jewish"       ~ 6.552531e-05,
  #   IDs == "16-2084567-TC-T" &
  #     name == "Latino/Admixed American"       ~ 6.552531e-05,
  #   TRUE                              ~ value_scaled
  # ))
  # [IDs = "1-939036-G-A", name = "Ashkenazi Jewish", "value_scaled"] <- 6.552531e-05

scaled_data1 %>% 
  filter(L_CHIP == "Yes") %>% 
  select(name, value_scaled) %>% 
  filter(value_scaled > 0) %>% 
  ggplot(aes(x=value_scaled, y=name, fill=name, color=name, height = after_stat(density))) +
  stat_density_ridges(alpha=0.3, stat="density", scale = 1.3, 
                      bandwidth = 2e-05,
                      quantile_lines = TRUE, quantiles = 2, # add meadian
                      vline_color = "black", vline_width = 0.25)+
  # geom_hline(yintercept = 2, color = "#CD9600")+
  # geom_hline(yintercept = 3, color = "#7CAE00")+
  # geom_hline(yintercept = 8, color = "#FF61CC")+
  labs(x= "Popultaion Allele Frequency", y= "Frequency Density by Population")+
  scale_fill_discrete(limits = c("African/African American",
                                 "Ashkenazi Jewish",
                                 "East Asian",
                                 "European (Finnish)",
                                 "European (non-Finnish)",
                                 "Latino/Admixed American",
                                 "South Asian",
                                 "Other"))+
  scale_color_discrete(limits = c("African/African American",
                                  "Ashkenazi Jewish",
                                  "East Asian",
                                  "European (Finnish)",
                                  "European (non-Finnish)",
                                  "Latino/Admixed American",
                                  "South Asian",
                                  "Other"))+
  scale_x_continuous(limits = c(0,0.00034),#c(2e-05,0.00013), # c(-0.00001,0.00015),
                     breaks = c(6.55e-05, 0.00013, 0.00019, 0.00026, 0.00032)
                     # labels = function(x) format(x, scientific = TRUE)
  )+
  scale_y_discrete(expand = c(0, 0.2))+
  geom_segment(data = dat,
               aes(y = name, x = value_scaled-0.0003213779, xend = value_scaled),
               show.legend = F, inherit.aes = F, color = c("#CD9600", "#7CAE00",
                                                           "#C77CFF", "#FF61CC"))+
  theme(legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ggsave("Figure 7 PAFs CH LCHIP_02042025.pdf",
       width = 5,
       height = 5,
       dpi = 600)


# Table S3----
library(tidyverse)
library(gtsummary)
table_s3 <- 
  readxl::read_xlsx(paste0(here::here(), "/Table S3_List of gnomAD identified SM and CH with annotation_01.06.2025.xlsx"))
table(table_s3$Gene.refGene)

table_s3 %>% 
  select(Gene.refGene) %>% 
  tbl_summary(sort = everything() ~ "frequency")


# Age rate ----
library(tidyverse)
library(ggbreak)
theme_set(theme_classic(base_size = 15))
theme_set(theme_classic())
CH_ <- 
  read_rds(paste0(here::here(), "/CH_age_dat.rds"))



CH_ %>% 
  select(-c(IDs : is_SM)) %>% 
  mutate(nrow = n()) %>% 
  group_by(nrow) %>% 
  summarise_at(c("<30",
                 "30-35","35-40","40-45","45-50",
                 "50-55","55-60","60-65","65-70",
                 "70-75","75-80", ">80"), ~ sum(.)) %>% 
  ungroup() %>% 
  pivot_longer(cols = -nrow) %>% 
  mutate(name = factor(name, levels = c("<30",
                                        "30-35","35-40","40-45","45-50",
                                        "50-55","55-60","60-65","65-70",
                                        "70-75","75-80", ">80"))) %>% 
  mutate(All = case_when(
    name == "<30" ~ 2547, 
    name == "30-35" ~ 3423,
    name == "35-40" ~ 4546,
    name == "40-45" ~ 8487,
    name == "45-50" ~ 10355,
    name == "50-55" ~ 12693,
    name == "55-60" ~ 11933,
    name == "60-65" ~ 10534,
    name == "65-70" ~ 8882,
    name == "70-75" ~ 5991,
    name == "75-80" ~ 4136,
    name == ">80" ~ 1935
  )) %>% 
  mutate(rate = value / All * 100) %>% 
  ggplot(aes(x = name, y = rate))+
  geom_bar(stat = "identity", fill= "#5DC863FF")
  #   labs(#title = "CH variants", 
  #        x= "Age",
  #        y= "Individuals (%)")+
  #   ylim(0, 20.1)

  #   "IDs"="All individuals", "ID"="All individuals of any genotype bin",
  #         "<30" = 2547, "30-35" = 3423,
  #         "35-40" = 4546,"40-45" = 8487,
  #         "45-50" = 10355,"50-55" = 12693,
  #         "55-60" = 11933,"60-65" = 10534,
  #         "65-70" = 8882,"70-75" = 5991,
  #         "75-80" = 4136,
  #         ">80" = 1935, .before = 1
  # ) %>%
  
  
  # 
  # mutate(sum = sum(value),
  #        perc = (value / sum) * 100) %>%
  
  
  
  
CH_ <- 
  read_rds(paste0(here::here(), "/age_CH_variants.rds"))
  
CH_M <- 
  read_rds(paste0(here::here(), "/age_CH_MCHIP_variants.rds"))
  
CH_L <- 
  read_rds(paste0(here::here(), "/age_CH_LCHIP_variants.rds"))

SM_ <- 
  read_rds(paste0(here::here(), "/age_SM_variants.rds"))

SM_M <- 
  read_rds(paste0(here::here(), "/age_SM_MCHIP_variants.rds"))

SM_L <- 
  read_rds(paste0(here::here(), "/age_SM_LCHIP_variants.rds"))

theme_set(theme_grey())


CH_ %>% 
  mutate(dat = "All CH") %>% 
  bind_rows(., CH_M %>% 
              mutate(dat = "M-CHIP")) %>% 
  bind_rows(., CH_L %>% 
              mutate(dat = "L-CHIP")) %>% 
  mutate(dat = factor(dat, levels = c("All CH",
                                      "M-CHIP",
                                      "L-CHIP"))) %>% 
  
  mutate(Percent = (value / (nrow * All))     ) %>% 
  
  ggplot(aes(x = factor(name), y = Percent, shape = dat, group = dat))+
  geom_point(color = "#5DC863FF", size =2)+ 
  geom_line(color = "#5DC863FF")+
  labs(x = "Age (years)", y = "Percent (%)")+
  expand_limits(y = c(0, 0.002))+
  scale_y_continuous(breaks=c(0, 0.0005, 0.0010, 0.0015, 0.0020), 
                     labels = c("0", "0.0005", "0.0010", "0.0015", "0.0020"))+
  scale_y_break(breaks = c(0.0005, 0.0010),
                scales = "fixed", expand = TRUE,
                ticklabels = c(0.00025, 0.0010, 0.0015, 0.0020)
  ) +
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

ggsave("New age figure CH.pdf",
       width = 5,
       height = 5,
       dpi = 600)


SM_ %>% 
  mutate(dat = "All SM") %>% 
  bind_rows(., SM_M %>% 
              mutate(dat = "SM in M-CHIP")) %>% 
  bind_rows(., SM_L %>% 
              mutate(dat = "SM in L-CHIP")) %>% 
  mutate(dat = factor(dat, levels = c("All SM",
                                      "SM in M-CHIP",
                                      "SM in L-CHIP"))) %>% 
  mutate(Percent = (value / (nrow * All))) %>% 
  ggplot(aes(x = factor(name), y = Percent, shape = dat, group = dat))+
  geom_point(color = "#21908CFF", size =2)+ 
  geom_line(color = "#21908CFF")+
  labs(x = "Age (years)", y = "Percent (%)")+
  expand_limits(y = c(0, 0.002))+
  scale_y_continuous(breaks=c(0, 0.0005, 0.0010, 0.0015, 0.0020), 
                     labels = c("0", "0.0005", "0.0010", "0.0015", "0.0020"))+
  scale_y_break(breaks = c(0.0005, 0.0010),
                scales = "fixed", expand = TRUE,
                ticklabels = c(0.00025, 0.0010, 0.0015, 0.0020)
  ) +
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

ggsave("New age figure SM.pdf",
       width = 5,
       height = 5,
       dpi = 600)














