# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the visualization of the data of Experiment A 
# Depicted in I) the main manuscript
# and II) the supplements



# ----- 2. Load in needed packages ---------------------------------------------
# create ggplots
library(ggplot2)
library(ggpmisc)

# to read in table with read_rds
library(tidyverse) 

# Grammar of Data Manipulation, e.g. mutate(); select(); filter(); summarise(); arrange() 
library(dplyr) # Grammar of Data Manipulation, e.g. mutate(); select(); filter(); summarise(); arrange() 

# A Forward-Pipe Operator
library(magrittr) 

# create frequency tables
library(rstatix)

# add text to ggplot - e.g., stars in heatmap
library("ggtext")

# to put different plots into one
library("patchwork")

library(cowplot)
# use labels = "AUTO" for auto addition of e.g., A, B labels


library(tidyr)
library(vegan)

#NMDS plot
library(ggrepel)
library(UpSetR)

# ----- 3. Data: Load and Edit -------------------------------------------------
## ---- 3.1. Manual feeding data -----------------------------------------------
# table with all manual feeding data
Feeding <- read_rds("in/Feeding_data.rds")

# Design Tables for graphical processing
# unite Particle and Biofilm to create additional column for both to use to visualize
 Feeding <- Feeding %>%
   unite(part_film, c(particle, biofilm), sep = "_", remove = FALSE)

# do the same for shape
 Feeding <- Feeding %>%
   unite(shape_film, c(shape, biofilm), sep = "_", remove = FALSE)

# level entries of new column
Feeding <- Feeding %>%
  transform(part_film = as.factor(part_film),
            shape_film = as.factor(shape_film))

Feeding$part_film <- factor(Feeding$part_film, 
                           levels = c("Copepods_blank", "Artemia_cysts_blank",
                                      "PE_steril", "PVC_steril", "PS_steril",
                                      "PA_steril", "PP_steril", "PET_steril",
                                      "Sand_steril",
                                      "PE_biofilm", "PVC_biofilm", "PS_biofilm",
                                      "PA_biofilm", "PP_biofilm", "PET_biofilm",
                                      "Sand_biofilm"))
levels(Feeding$part_film)

Feeding$shape_film <- factor(Feeding$shape_film, 
                            levels = c("Artemia_blank", "Copepods_blank",
                                       "fiber_steril", "fragment_steril",
                                       "Sand_steril",
                                       "fiber_biofilm", "fragment_biofilm",
                                       "Sand_biofilm"))
levels(Feeding$shape_film)

Feeding <- Feeding %>%
  mutate(incorporation = ifelse(egestion == "yes", "no", "yes"))


## ---- 3.2. Polyp activity measurements ---------------------------------------
polyp_activity <- read_rds("in/polyp_activity_graph_data.rds") %>%
  mutate(tp = case_when(tp == "1" ~ "0",
                        tp == "2" ~ "1",
                        tp == "2.5" ~ "1.2",
                        tp == "3" ~ "2",
                        tp == "3.5" ~ "2.2",
                        tp == "4" ~ "3",
                        tp == "4.5" ~ "3.2",
                        tp == "5"~ "4",
                        tp == "5.5"~ "4.2")) %>%
  transform(tp = as.numeric(tp))

## ---- 3.3. Build up theme ----------------------------------------------------
# Create labels for species
spec_labs <- c("P. verrucosa", "S. pistillata")
names(spec_labs) <- c("Pve", "Spi")

# Create labels for Response
response_labs <-  c("No reaction",
                    "Release",
                    "Egestion", "Incorporation")
names(response_labs) <- c("neg_react",
                          "neg_ingest",
                          "egest", "incorporat")

# Create labels for offered particles
offered_labs <- c("Artemia", "Copepods",
                  "PE", "PVC", "PS",
                  "PA", "PP", "PET",
                  "Sand",
                  "Biofouled PE", "Biofouled PVC", "Biofouled PS",
                  "Biofouled PA", "Biofouled PP", "Biofouled PET",
                  "Biofouled Sand")
names(offered_labs) <- c("Artemia_cysts_blank", "Copepods_blank",
                         "PE_steril", "PVC_steril", "PS_steril",
                         "PA_steril", "PP_steril", "PET_steril",
                         "Sand_steril",
                         "PE_biofilm", "PVC_biofilm", "PS_biofilm",
                         "PA_biofilm", "PP_biofilm", "PET_biofilm",
                         "Sand_biofilm")


# ----- 4. Question 1 ----------------------------------------------------------
## ---- 4.1. Edit Tables for graphical processing ------------------------------
Q1 <- Feeding %>%
  filter(biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger"),
         student %in% c ("Celine", "Verena", "Moritz", "Tim"))



## ---- 4.2. Create Plot -------------------------------------------------------
### --- 4.2.1. Figure 3 --------------------------------------------------------
Q1_Pos_Neg_2 <- Q1 %>%
  pivot_longer(cols=c(reaction, ingestion, incorporation), 
               names_to = "variable",
               values_to = "response", values_drop_na = FALSE) 

# remove "positive response" = Reaction and Ingestion, as otherwise they will be included double
Q1_Freq_2 <- Q1_Pos_Neg_2 %>%
  freq_table(spec, part_film, trigger, variable, response)

# level entries of new column
Q1_Freq_2 <- Q1_Freq_2 %>%
  transform(variable = as.factor(variable))

Q1_Freq_2$variable <- factor(Q1_Freq_2$variable, 
                           levels = c("reaction", 
                                      "ingestion",
                                      "incorporation"))
levels(Q1_Freq_2$response)

Q1_Freq_2 <- Q1_Freq_2 %>%
  mutate(fill_var = interaction(variable, response))

# Create a color scheme
color_scheme <- c("reaction.no" = "#F1CD7E", 
                  "reaction.yes" = "#CE891E",
                  "ingestion.no" = "#529DA3", 
                  "ingestion.yes" = "#164242",
                  "incorporation.no" = "#A3A152", 
                  "incorporation.yes" = "#5E5A31")

# Create specific labels for Responses for legend
response_labs_Q1 <-  c("No reaction", "Reaction",
                    "Release", "Ingestion",
                    "Egestion", "Incorporation")
names(response_labs_Q1) <- c("reaction.no", "reaction.yes",
                          "ingestion.no", "ingestion.yes",
                          "incorporation.no", "incorporation.yes")

levels(Q1_Freq_2$fill_var)

Q1_Freq_2$fill_var<- factor(Q1_Freq_2$fill_var, 
                             levels = c("reaction.no", "reaction.yes",
                                        "ingestion.no", "ingestion.yes",
                                        "incorporation.no", "incorporation.yes"))

# Create specific labels for Responses for axis title
Var_Tit <-  c("Reaction", "Ingestion", "Incorporation")
names(Var_Tit) <- c("reaction", "ingestion", "incorporation")

labels_n <- Q1_Freq_2 %>%
  filter(response == "yes") %>%
  group_by(part_film, variable, spec) %>%
  slice_min(order_by = fill_var, with_ties = FALSE) %>%   # hängt von Reihenfolge ab
  mutate(ypos = (prop)/100,       # unterster Block startet immer bei 0
         label = n)  # Zahl, die du dranschreiben willst

Q1_graph <- ggplot(Q1_Freq_2, aes(part_film, prop, fill=fill_var)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(variable ~ spec, 
             labeller = labeller(spec = spec_labs, variable = Var_Tit)) +
  # add lines to split natural food from MP and sand
  geom_vline(xintercept = 2.5) +
  geom_vline(xintercept = + 8.5) +
  labs(x = "Offered particle", y = "Response (%)", fill = "Response") +
  scale_fill_manual(values = color_scheme,
                    labels = response_labs_Q1) +
  scale_x_discrete(labels = offered_labs) +
  # convert y axis to percent
  scale_y_continuous(labels = function(x) paste0(x*100),
                     breaks = seq(0, 1, by = 0.25)) +
  # Zahl beim untersten Block hinzufügen
  geom_text(data = labels_n,
            aes(x = part_film, y = ypos+0.11, label = label),
            vjust = 1.2,    # leicht unterhalb des Blocks
            size = 3,
            inherit.aes = FALSE) +
  # design theme
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm')) #change legend key size
Q1_graph

ggsave("out/Figure3.png", plot = Q1_graph,
       scale = 1, width = 20, height = 20, units = c("cm"),
       dpi = 900, limitsize = TRUE) 




# ----- 5. Question 2 ----------------------------------------------------------
## ---- 5.1. Edit Tables for graphical processing ------------------------------
Q2 <- Feeding %>%
  filter(biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger"),
         particle %in% c ("PP", "PA"))

# create a frequency table for shape
Q2_Pos_Neg <- Q2 %>%
  mutate(reaction = case_when(reaction == "yes"~ "pos_react", 
                              reaction == "no" ~ "neg_react"),
         ingestion = case_when(ingestion == "yes" ~ "pos_ingest", # positive ingestion
                               ingestion == "no" ~ "neg_ingest"), # negative ingestion --> release
         incorporation = case_when(incorporation == "no" ~ "egest", # positive egestion
                              incorporation == "yes" ~ "incorporat")) %>% # negative egestion --> still in polyp
  pivot_longer(cols=c(reaction, ingestion, incorporation), 
               names_to = "variable",
               values_to = "response", values_drop_na = FALSE) 

Q2_Freq <- Q2_Pos_Neg %>%
  freq_table(spec, particle, shape, trigger, variable, response) %>%
  transform(variable = as.factor(variable)) %>% 
  mutate(compare = interaction(particle, shape))

Q2_Freq$response <- factor(Q2_Freq$response, c("neg_react", "pos_react",
                                               "neg_ingest", "pos_ingest",
                                               "egest", "incorporat"))
Q2_Freq$variable <- factor(Q2_Freq$variable, 
                           levels = c("reaction", 
                                      "ingestion",
                                      "incorporation"))
Q2_Freq$compare <- factor(Q2_Freq$compare, 
                           levels = c("PA.fragment", 
                                      "PA.fiber",
                                      "PP.fragment",
                                      "PP.fiber"))

# Create a color scheme
color_scheme_Q2 <- c("neg_react" = "#F1CD7E", 
                  "pos_react" = "#CE891E",
                  "neg_ingest" = "#529DA3", 
                  "pos_ingest" = "#164242",
                  "egest" = "#A3A152", 
                  "incorporat" = "#5E5A31")

# Create specific labels for Responses for legend
response_labs_Q2 <-  c("No reaction", "Reaction",
                       "Release", "Ingestion",
                       "Egestion", "Incorporation")
names(response_labs_Q2) <- c("neg_react", "pos_react",
                             "neg_ingest", "pos_ingest",
                             "egest", "incorporat")
offered_labs_Q2 <-  c("PA Fragment", "PA Fiber",
                      "PP Fragment","PP Fiber")
names(offered_labs_Q2) <- c("PA.fragment", "PA.fiber",
                            "PP.fragment","PP.fiber")

# Create specific labels for Responses for axis title
Var_Tit <-  c("Reaction", "Ingestion", "Incorporation")
names(Var_Tit) <- c("reaction", "ingestion", "incorporation")

labels_n2 <- Q2_Freq %>%
  filter(response %in% c ("pos_react", "pos_ingest", "incorporat")) %>%
  group_by(compare, variable, spec) %>%
  slice_min(order_by = particle, with_ties = FALSE) %>%   # hängt von Reihenfolge ab
  mutate(ypos = (prop)/100,       # unterster Block startet immer bei 0
         label = n)  # Zahl, die du dranschreiben willst

Q2_graph <- ggplot(Q2_Freq, aes(compare, prop, fill=response)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(variable ~ spec, 
             labeller = labeller(spec = spec_labs, variable = Var_Tit)) +
  # add lines to split natural food from MP and sand
  geom_vline(xintercept = 2.5) +
  labs(x = "Offered particle", y = "Response (%)", fill = "Response") +
  scale_fill_manual(values = color_scheme_Q2,
                    labels = response_labs_Q2) +
  scale_x_discrete(labels = offered_labs_Q2) +
  # convert y axis to percent
  scale_y_continuous(labels = function(x) paste0(x*100)) +
  # Zahl beim untersten Block hinzufügen
  geom_text(data = labels_n2,
            aes(x = compare, y = ypos+0.11, label = label),
            vjust = 1.2,    # leicht unterhalb des Blocks
            size = 3,
            inherit.aes = FALSE) +
  # design theme
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm')) #change legend key size
Q2_graph

ggsave("out/Figure4.png", plot = Q2_graph,
      scale = 1, width = 12, height = 18, units = c("cm"),
     dpi = 900, limitsize = TRUE) 




# ----- 6. Question 3 ----------------------------------------------------------
## ---- 6.1. Edit Tables for graphical processing ------------------------------
Q3 <- Feeding %>%
  filter(student %in% c ("Celine", "Verena", "Moritz", "Tim"),
         particle %in% c ("PE", "PS", "PVC", "PA", "PP", "PET")) %>%
  mutate(biotic = interaction(biofilm, trigger))

# level biotic factors
Q3$biotic<- factor(Q3$biotic, 
                           levels = c("steril.NoTrigger", 
                                      "biofilm.NoTrigger",
                                      "steril.Trigger",
                                      "biofilm.Trigger"))
# level shapes
Q3$shape<- factor(Q3$shape, 
                   levels = c("fragment", 
                              "fiber"))


Q3_Pos_Neg <- Q3 %>%
  pivot_longer(cols=c(reaction, ingestion, incorporation), 
               names_to = "variable",
               values_to = "response", values_drop_na = FALSE) 

# remove "positive response" = Reaction and Ingestion, as otherwise they will be included double
Q3_Freq <- Q3_Pos_Neg %>%
  freq_table(spec, biotic, shape, trigger, variable, response)

# level entries of response variables
Q3_Freq <- Q3_Freq %>%
  transform(variable = as.factor(variable))

Q3_Freq$variable <- factor(Q3_Freq$variable, 
                             levels = c("reaction", 
                                        "ingestion",
                                        "incorporation"))

Q3_Freq <- Q3_Freq %>%
  mutate(fill_var = interaction(variable, response))

# level responses for graph
Q3_Freq$fill_var <- factor(Q3_Freq$fill_var, 
                      levels = c("reaction.no", 
                                 "reaction.yes",
                                 "ingestion.no", 
                                 "ingestion.yes",
                                 "incorporation.no", 
                                 "incorporation.yes"))
levels(Q3_Freq$response)

# Create specific labels for Response
response_labs_Q3 <-  c("No reaction", "Reaction",
                       "Release", "Ingestion",
                       "Egestion", "Incorporation")
names(response_labs_Q3) <- c("reaction.no", "reaction.yes",
                             "ingestion.no", "ingestion.yes",
                             "incorporation.no", "incorporation.yes")

# Create a color scheme
color_scheme <- c("reaction.no" = "#F1CD7E", 
                  "reaction.yes" = "#CE891E",
                  "ingestion.no" = "#529DA3", 
                  "ingestion.yes" = "#164242",
                  "incorporation.no" = "#A3A152", 
                  "incorporation.yes" = "#5E5A31")

# Create specific labels for Responses for axis title
Var_Tit <-  c("Reaction", "Ingestion", "Incorporation")
names(Var_Tit) <- c("reaction", "ingestion", "incorporation")

# Create specific labels for Responses for axis title
shape_labs <-  c("Fragment", "Fiber")
names(shape_labs) <- c("fragment", "fiber")

# Create labels for biotic factors
biotic_labs <- c("unfouled MPs", "biofouled MPs",
                 "unfouled MPs + trigger", "biofouled MPs + trigger")
names(biotic_labs) <- c("steril.NoTrigger", "biofilm.NoTrigger",
                        "steril.Trigger", "biofilm.Trigger")

labels_n3 <- Q3_Freq %>%
  filter(response == "yes") %>%
  group_by(biotic, shape, variable, spec) %>%
  slice_min(order_by = fill_var, with_ties = FALSE) %>%   # hängt von Reihenfolge ab
  mutate(ypos = (prop)/100,       # unterster Block startet immer bei 0
         label = n)  # Zahl, die du dranschreiben willst

## ---- 6.2. Create Plot -------------------------------------------------------
Q3_Graph <- ggplot(Q3_Freq, aes(biotic, prop, fill=fill_var)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(variable + shape ~ spec, 
             labeller = labeller(spec = spec_labs, variable = Var_Tit,
                                 shape = shape_labs)) +
  labs(x = "Offered particle", y = "Response (%)", fill = "Response") +
  scale_fill_manual(values = color_scheme,
                    labels = response_labs_Q3) +
  scale_x_discrete(labels = biotic_labs) +
  # convert y axis to percent
  scale_y_continuous(labels = function(x) paste0(x*100),
                     breaks = seq(0, 1, by = 0.25)) +
  # Zahl beim untersten Block hinzufügen
  geom_text(data = labels_n3,
            aes(x = biotic, y = ypos+0.17, label = label),
            vjust = 1.2,    # leicht unterhalb des Blocks
            size = 3,
            inherit.aes = FALSE) +
  # design theme
  theme_classic() +
theme(panel.background = element_rect(color = "black"),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(face = "italic", size = 10),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, 'cm')) #change legend key size

Q3_Graph

ggsave("out/Figure5.png", plot = Q3_Graph,
       scale = 1, width = 12, height = 24, units = c("cm"),
       dpi = 900, limitsize = TRUE) 


# ----- 7. Question 4 ----------------------------------------------------------
## ---- 7.1. Create Table ------------------------------------------------------
Q4 <- read_csv2("in/z-values.csv") %>%
  transform(z_value = as.numeric(z_value),
            OR = as.numeric(OR),
            log_OR = as.numeric(log_OR),
            graph_OR = as.numeric(graph_OR),
            p_code= as.factor(p_code),
            Response = as.factor(Response)) %>%
  select(-Specifications, -n, -vs) %>%
  # add color for z-value according to significance
  mutate(Color = case_when(p_code == "***"~ "#070B0D",
                           p_code == "**"~ "#152128",
                           p_code == "*"~ "#233743",
                           p_code == "."~ "#6C97B2",
                           p_code == "-" ~ "#A1BCCE"))
  
  
# level responses in right order
Q4$Response <- factor(Q4$Response, 
                              levels = c("Reaction", "Ingestion", "Incorporation"))

# level Contrast 1 (left side) in right order
Q4$Contrast_1 <- factor(Q4$Contrast_1, 
                      levels = c("Food", "Pve", "No Trigger",
                                 "unfouled", "Fragment", "Polymer A"))

# level Contrast 2 (right side) in right order
Q4$Contrast_2 <- factor(Q4$Contrast_2, 
                      levels = c("MP", "Spi", "Trigger",
                                 "biofouled", "Fiber", "Polymer B"))

## ---- 7.2. Create Plot -------------------------------------------------------
# create labels for left and right axis
left_labs <- c("Food",
               "<i>P. verrucosa<i>",
               "No Trigger",
               "Unfouled",
               "Fragment", 
               "Polymer A")
right_labs <- c("MP", 
                "<i>S. pistillata<i>",
                "Trigger",
                "Biofouled",
                "Fiber", 
                "Polymer B")

# separate plots to handle left and right axis
Q4_react <- Q4 %>%
  filter(Response %in% c("Reaction"))

Q4_react_graph <-  
  ggplot(Q4_react, aes(x = graph_OR, y = Contrast_1, color = Color)) +  
  facet_grid(~ Response) +
  # add horizontal line for z = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#AFC6D5", alpha = 0.7) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "#AFC6D5", alpha = 0.7) +
  geom_point()+ 
  geom_linerange(aes(xmin = 0, xmax = graph_OR), linetype = "dotted") +
  scale_color_identity() +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_discrete(labels = left_labs) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
       # axis.text.y = element_text(size = 10),
        axis.text.y = element_markdown(size = 10),
        legend.text = element_text(size = 10)) +
  xlab("-log(Odds Ratio)")

Q4_react_labs <- Q4_react %>%
  mutate(label = case_when(p_code == "***"~ "***",
                           p_code == "**"~ "**",
                           p_code == "*"~ "*",
                           p_code == "."~ ".",
                           p_code == "-" ~ NA))
  
Q4_react_graph_an <- Q4_react_graph + geom_text(data = Q4_react_labs, mapping = aes(label = label),
                           nudge_y = 0.1)
Q4_react_graph_an 


Q4_ing <- Q4 %>%
  filter(Response %in% c("Ingestion"))

Q4_ing_graph <- 
ggplot(Q4_ing, aes(x = graph_OR, y = Contrast_1, color = Color)) +  
facet_grid(~ Response) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#AFC6D5", alpha = 0.7) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "#AFC6D5", alpha = 0.7) +
#  annotate("segment", x = 0, xend = 5, y = 5, yend = 5, linetype = "dotted",
 #          col = "#A1BCCE", linewidth = 0.5) +
  geom_point()+ 
  geom_linerange(aes(xmin = 0, xmax = graph_OR), linetype = "dotted") +
  scale_color_identity() +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_discrete(labels = left_labs) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        legend.text = element_text(size = 10)) +
  xlab("-log(Odds Ratio)")

Q4_ing_graph

Q4_ing_labs <- Q4_ing %>%
  mutate(label = case_when(p_code == "***"~ "***",
                           p_code == "**"~ "**",
                           p_code == "*"~ "*",
                           p_code == "."~ ".",
                           p_code == "-" ~ NA))

Q4_ing_graph_an <- Q4_ing_graph + geom_text(data = Q4_ing_labs, mapping = aes(label = label),
                                                nudge_y = 0.1)
Q4_ing_graph_an

Q4_incor <- Q4 %>%
  filter(Response %in% c("Incorporation"))

Q4_incor_graph <-  
  ggplot(Q4_incor, aes(x = graph_OR, y = Contrast_1, color = Color)) +  
  facet_grid(~ Response) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "#AFC6D5", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#AFC6D5", alpha = 0.7) +
  geom_point()+ 
  geom_linerange(aes(xmin = 0, xmax = graph_OR), linetype = "dotted") +
  scale_color_identity() +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_discrete(position = "right", labels = right_labs) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_markdown(size = 10),
        legend.text = element_text(size = 10)) +
  xlab("-log(Odds Ratio)")

  
Q4_incor_graph

Q4_incor_labs <- Q4_incor %>%
  mutate(label = case_when(p_code == "***"~ "***",
                           p_code == "**"~ "**",
                           p_code == "*"~ "*",
                           p_code == "."~ ".",
                           p_code == "-" ~ NA))

Q4_incor_graph_an <- Q4_incor_graph + geom_text(data = Q4_incor_labs, mapping = aes(label = label),
                                            nudge_y = 0.1)

Q4_graph <- Q4_react_graph_an + Q4_ing_graph_an + Q4_incor_graph_an
Q4_graph

ggsave("out/Figure6.png", plot = Q4_graph,
       scale = 1, width = 20, height = 10, units = c("cm"),
       dpi = 600, limitsize = TRUE)  


# ---- 8. Associated molecules  ------------------------------------------------
molecules <- read_csv2("in/Molecules.csv")
data_wide <- molecules %>%
  unite("Polymer_Status", Polymer, State, sep = "_") %>%
  pivot_longer(cols = starts_with("C"), 
               names_to = "molecule", 
               values_to = "measured") %>%
  pivot_wider(names_from = Polymer_Status, values_from = measured, 
              values_fill = list(measured = 0))

upset_data <- data_wide %>%
  as.data.frame()  # Stellen sicher, dass es ein DataFrame ist

set_colors <- c(
  "#1F1F1F", "#2D2DB4", "#1F1F1F",
  "#888781", "#888781", "#BA3B3B",  
  "#1F1F1F", "#CEC1A6", "#CEC1A6",  
  "#CFCA4F", "#CFCA4F", "#2D2DB4",  
  "#2D2DB4", "#BA3B3B", "#888781",  
  "#CEC1A6", "#CFCA4F", "#BA3B3B"   
)

png(paste0("out/Figure7A.png"), units="cm", width=16, height=16, res=600)
upset(upset_data, 
      nintersects = 40, 
      nsets = 18, 
      mb.ratio = c(0.25, 0.75),
      order.by = "freq", 
      text.scale = 1.2, 
      point.size = 2.5, 
      line.size = 1,
      main.bar.color = "#788AA3", 
      sets.bar.color = set_colors,
      matrix.color = "#333333")  
dev.off()



# Convert categorical variables to row names
rownames(molecules) <- paste(molecules$Polymer, molecules$State, sep = "_")

# Remove non-numeric columns (Polymer, State) and the summary column (15)
data <- molecules[, -c(1,2,15)]
# Remove constant columns (zero variance)
data <- data[, apply(data, 2, function(x) length(unique(x)) > 1)]

# Compute the dissimilarity matrix (using Jaccard index)
dist_matrix <- vegdist(data, method = "jaccard")

# Perform NMDS (Non-metric Multidimensional Scaling)
set.seed(852)
nmds_result <- metaMDS(dist_matrix, k = 2)  # k=2 for 2D plot
stress_value <- round(nmds_result$stress, 4) # get stress value
# [1] 0.0455304 - perfect

# Extract NMDS scores
nmds_scores <- as.data.frame(nmds_result$points) %>%
  mutate(Polymer = molecules$Polymer,
         State   = molecules$State)

# Compute convex hull for each Polymer group
hull_data_polymer <- nmds_scores %>%
  group_by(Polymer) %>%
  slice(chull(MDS1, MDS2))
hull_data_polymer$Polymer <- factor(hull_data_polymer$Polymer, 
                                    levels = c("PE", "PVC", "PS",
                                               "PA", "PP", "PET"))

hull_data_polymer$State <- factor(hull_data_polymer$State, 
                                  levels = c("unfouled", "biofouled", "leachate"))

# NMDS Plot with Convex Hulls by Polymer
Moleks_B <- ggplot(nmds_scores, aes(x = MDS1, y = MDS2, color = Polymer, shape = State)) +
  # Convex Hulls by Polymer
  geom_polygon(data = hull_data_polymer, aes(group = Polymer, fill = Polymer), alpha = 0.25, linetype = "dotted") +
  # NMDS Points
  geom_point(size = 3, alpha = 0.9) +
  annotate("text", x = max(nmds_scores$MDS1), y = -1.5, 
           label = paste("Stress =", stress_value), hjust = 1, size = 3) +
  labs(x = "NMDS1", y = "NMDS2", color = "Polymer", shape = "State", fill = "Polymer") +
  scale_shape_manual(name = "Condition", values=c(16, 1, 8)) +
  scale_color_manual(name = "Polymer", 
                     values = c("PE" = "#1F1F1F", "PVC" = "#888781", "PS" = "#CEC1A6", 
                                "PP" = "#CFCA4F", "PET" = "#2D2DB4","PA" = "#BA3B3B")) +
  scale_fill_manual(name = "Polymer", 
                    values = c("PE" = "#1F1F1F", "PVC" = "#888781", "PS" = "#CEC1A6", 
                               "PP" = "#CFCA4F", "PET" = "#2D2DB4","PA" = "#BA3B3B")) +
  theme_classic() +
  theme(panel.background = element_rect(color = "black", linewidth = 1),
        legend.key = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"), # Adjust key size
        legend.position = "right")
Moleks_B


ggsave("out/Figure7B.png", plot = last_plot(),
       scale = 1, width = 14, height = 12, units = c("cm"),
       dpi = 600, limitsize = TRUE)


# Impact of State
# PERMANOVA with Jaccard
set.seed(852)  # Setze eine feste Zahl, damit die Permutationen gleich bleiben
permanova_result <- adonis2(dist_matrix ~ State, data = molecules, permutations = 4999)
print(permanova_result)
# OUTPUT: Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 4999
# adonis2(formula = dist_matrix ~ State, data = molecules, permutations = 4999)
# Df SumOfSqs      R2      F Pr(>F)  
# Model     2   1.1898 0.18817 1.7384 0.0566 .
# Residual 15   5.1332 0.81183                
# Total    17   6.3230 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# Impact of Polymer
# PERMANOVA with Jaccard
set.seed(852)  # Setze eine feste Zahl, damit die Permutationen gleich bleiben
permanova_result <- adonis2(dist_matrix ~ Polymer, data = molecules, permutations = 4999)
print(permanova_result)
# OUTPUT: Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 4999
# adonis2(formula = dist_matrix ~ Polymer, data = molecules, permutations = 4999)
# Df SumOfSqs    R2      F Pr(>F)    
# Model     5   3.8697 0.612 3.7855  2e-04 ***
# Residual 12   2.4533 0.388                  
# Total    17   6.3230 1.000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# ----- 9. Question 5 ----------------------------------------------------------
## ---- 9.1. Edit Tables for graphical processing ------------------------------
Q5 <- polyp_activity %>%
  filter(treat %in% c("MP", "MP+food")) %>%
  group_by(tp) %>%
  # Calculate the mean of the "Frequency" column for each group
  summarise_at(vars(ranks),
               list(Mean_ranks = mean))

Q5_mean <- polyp_activity %>%
  filter(treat %in% c("MP", "MP+food"),
         tp != "0") %>%
  group_by(spec, ID, col, treat, tp) %>%
  get_summary_stats(ranks, type = "mean")
  

## ---- 9.2. Create Plot -------------------------------------------------------
# Create labels for species
spec_labs <- c("P. verrucosa", "S. pistillata")
names(spec_labs) <- c("Pve", "Spi")

# Create labels for treatments
treat_labs <- c("control", "MP", "MP+Food")
names(treat_labs) <- c("control", "MP", "MP+food")

# Create labels for foodnote
food_labs <- c("MP+Food fed", "all fed")
names(food_labs) <- c("food", "all")

# create color scheme
color_scheme <- c(#"#F1DEDE", # control 
                  "#5D576B", # MP
                  "#D496A7") # MP+food


Q5_Graph <- 
  Q5_mean %>% 
  ggplot(aes(x=treat, y=mean)) +
    geom_violin(alpha=0.6,
                aes(fill = treat, color = treat)) +
    geom_boxplot(width=0.1, aes(color=treat)) +
  facet_grid(. ~ spec, 
             labeller = labeller(spec = spec_labs)) +
  scale_fill_manual(values = color_scheme) +
    scale_color_manual(values = color_scheme) +
  labs(x = "Treatment",
       y = "Polypactivity",
       title = "") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none")

Q5_Graph

ggsave("out/Figure8.png", plot = Q5_Graph,
       scale = 1, width = 10, height = 10, units = c("cm"),
       dpi = 600, limitsize = TRUE)






# ----- X. Supplements ---------------------------------------------------------
## ------- Time graphs ---------------------------------------------------------
### ------ SQ1 -----------------------------------------------------------------
Time_Q1 <- Q1 %>% 
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# level particles
Time_Q1$particle <- factor(Time_Q1$particle, 
                              levels = c("Artemia_cysts", "Copepods",
                                         "PE", "PVC", "PS",
                                         "PA", "PP", "PET",
                                         "Sand"))

# Create labels for offered particles
offered_labs_Time1 <- c("Artemia", "Copepods",
                  "PE", "PVC", "PS",
                  "PA", "PP", "PET",
                  "Sand")
names(offered_labs_Time1) <- c("Artemia_cysts", "Copepods",
                         "PE", "PVC", "PS",
                         "PA", "PP", "PET",
                         "Sand")

Q1_react_dur <- Time_Q1 %>% 
  ggplot(aes(x=particle, y=reaction_dur_num)) +
  geom_violin(alpha=0.6,
              aes(fill = particle, color = particle)) +
  geom_boxplot(width=0.1, aes(color=particle)) +
  facet_grid(. ~ spec, 
             labeller = labeller(spec = spec_labs)) +
  scale_fill_manual(values = c("#CE891E", "#CE891E", "#CE891E", "#CE891E", "#CE891E", "#CE891E","#CE891E", "#CE891E", "#CE891E")) +
  scale_color_manual(values = c("#CE891E", "#CE891E", "#CE891E", "#CE891E", "#CE891E", "#CE891E","#CE891E", "#CE891E", "#CE891E")) +
  scale_y_continuous(trans='log2', breaks=c(10, 30, 60, 120, 300, 600, 1800, 3600),
                     labels=c("10 sec", "30 sec", "1 min", "2 min", "5 min", "10 min", "30 min", "1 h")) +
  scale_x_discrete(labels = offered_labs_Time1) +
  labs(x = "Offered particle",
       y = "Reaction time",
       title = "") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none")
Q1_react_dur

Q1_ing_dur <- Time_Q1 %>% 
  ggplot(aes(x=particle, y=ingestion_dur_num)) +
  geom_violin(alpha=0.6,
              aes(fill = particle, color = particle)) +
  geom_boxplot(width=0.1, aes(color=particle)) +
  facet_grid(. ~ spec, 
             labeller = labeller(spec = spec_labs)) +
  scale_fill_manual(values = c("#164242", "#164242", "#164242", "#164242", "#164242", "#164242","#164242", "#164242", "#164242")) +
  scale_color_manual(values = c("#164242", "#164242", "#164242", "#164242", "#164242", "#164242","#164242", "#164242", "#164242")) +
  scale_y_continuous(trans='log2', breaks=c(10, 30, 60, 120, 300, 600, 1800, 3600),
                     labels=c("10 sec", "30 sec", "1 min", "2 min", "5 min", "10 min", "30 min", "1 h")) +
  scale_x_discrete(labels = offered_labs_Time1) +
  labs(x = "Offered particle",
       y = "Ingestion time",
       title = "") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none")
Q1_ing_dur

Q1_Time <- Q1_react_dur / Q1_ing_dur

ggsave("out/Figure_S3.png", plot = Q1_Time,
       scale = 1, width = 16, height = 28, units = c("cm"),
       dpi = 600, limitsize = TRUE)




### ------ SQ2 -----------------------------------------------------------------
Time_Q2 <- Q2 %>% 
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# level particles
Time_Q2$shape <- factor(Time_Q2$shape, 
                           levels = c("fragment",
                                      "fiber"))

Q2_react_dur <- Time_Q2 %>% 
  ggplot(aes(x=shape, y=reaction_dur_num)) +
  geom_violin(alpha=0.6,
              aes(fill = shape, color = shape)) +
  geom_boxplot(width=0.1, aes(color=shape)) +
  facet_grid(particle ~ spec, 
             labeller = labeller(spec = spec_labs)) +
  scale_fill_manual(values = c("#CE891E", "#CE891E")) +
  scale_color_manual(values = c("#CE891E", "#CE891E")) +
  scale_y_continuous(trans='log2', breaks=c(10, 30, 60, 120, 300, 600, 1800, 3600),
                     labels=c("10 sec", "30 sec", "1 min", "2 min", "5 min", "10 min", "30 min", "1 h")) +
 # scale_x_discrete(labels = offered_labs_Time1) +
  labs(x = "Offered particle",
       y = "Reaction time",
       title = "") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none")
Q2_react_dur

Q2_ing_dur <- Time_Q2 %>% 
  ggplot(aes(x = shape, y=ingestion_dur_num)) +
  geom_violin(alpha=0.6,
              aes(fill = shape, color = shape)) +
  geom_boxplot(width=0.1, aes(color=shape)) +
  facet_grid(particle ~ spec, 
             labeller = labeller(spec = spec_labs)) +
  scale_fill_manual(values = c("#164242", "#164242")) +
  scale_color_manual(values = c("#164242", "#164242")) +
  scale_y_continuous(trans='log2', breaks=c(10, 30, 60, 120, 300, 600, 1800, 3600),
                     labels=c("10 sec", "30 sec", "1 min", "2 min", "5 min", "10 min", "30 min", "1 h")) +
 # scale_x_discrete(labels = offered_labs_Time1) +
  labs(x = "Offered particle",
       y = "Ingestion time",
       title = "") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none")
Q2_ing_dur

Q2_Time <- Q2_react_dur / Q2_ing_dur

ggsave("out/Figure_S4.png", plot = Q2_Time,
       scale = 1, width = 16, height = 28, units = c("cm"),
       dpi = 600, limitsize = TRUE)


### ------ SQ3 -----------------------------------------------------------------
Time_Q3 <- Q3 %>% 
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# Create labels for biotic factors
biotic_labs <- c("MP", "biofouled MP",
                 "MP + trigger", "biofouled MP + trigger")
names(biotic_labs) <- c("steril.NoTrigger", "biofilm.NoTrigger",
                        "steril.Trigger", "biofilm.Trigger")


Q3_react_dur <- Time_Q3 %>% 
  ggplot(aes(x=biotic, y=reaction_dur_num)) +
  geom_violin(alpha=0.6,
              aes(fill = biotic, color = biotic)) +
  geom_boxplot(width=0.1, aes(color=biotic)) +
  facet_grid(shape ~ spec, 
             labeller = labeller(spec = spec_labs)) +
  scale_fill_manual(values = c("#CE891E", "#CE891E", "#CE891E", "#CE891E")) +
  scale_color_manual(values = c("#CE891E", "#CE891E", "#CE891E", "#CE891E")) +
  scale_y_continuous(trans='log2', breaks=c(10, 30, 60, 120, 300, 600, 1800, 3600),
                     labels=c("10 sec", "30 sec", "1 min", "2 min", "5 min", "10 min", "30 min", "1 h")) +
  scale_x_discrete(labels = biotic_labs) +
  labs(x = "Offered particle",
       y = "Reaction time",
       title = "") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none")
Q3_react_dur

Q3_ing_dur <- Time_Q3 %>% 
  ggplot(aes(x=biotic, y=ingestion_dur_num)) +
  geom_violin(alpha=0.6,
              aes(fill = biotic, color = biotic)) +
  geom_boxplot(width=0.1, aes(color=biotic)) +
  facet_grid(shape ~ spec, 
             labeller = labeller(spec = spec_labs)) +
  scale_fill_manual(values = c("#164242", "#164242", "#164242", "#164242")) +
  scale_color_manual(values = c("#164242", "#164242", "#164242", "#164242")) +
  scale_y_continuous(trans='log2', breaks=c(10, 30, 60, 120, 300, 600, 1800, 3600),
                     labels=c("10 sec", "30 sec", "1 min", "2 min", "5 min", "10 min", "30 min", "1 h")) +
  scale_x_discrete(labels = biotic_labs) +
  labs(x = "Offered particle",
       y = "Ingestion time",
       title = "") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none")
Q3_ing_dur

Q3_Time <- Q3_react_dur / Q3_ing_dur

ggsave("out/Figure_S5.png", plot = Q3_Time,
       scale = 1, width = 16, height = 28, units = c("cm"),
       dpi = 600, limitsize = TRUE)







## ------- Mean tables ---------------------------------------------------------
write_csv2(Q1_Freq_2, "out/Mean_values.csv")

