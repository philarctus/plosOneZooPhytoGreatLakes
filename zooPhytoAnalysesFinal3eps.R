# Zooplankton-phytoplankton associations GLNPO this version has edibility data
# Feb 2021, updated zooplankton data analyses May 19, January 2022, July 2022
# author: katya
# R version 4.0.2 "Taking Off Again"

# Set working directory to the file location (only RStudio) 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Part 1 Merger
install_load <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}
# Required packages
install_load(c("ggplot2", "vegan", "plyr", "dplyr", "mgcv", "stringi", "viridis", "randomForest",
               "readxl", "tidyr", "reshape2", "data.table", "labdsv", "vegan", "purrr", "compare", "ggdensity"))

#########################################
# # Get phyto data from canonical files, merging 2018 vs the rest, calculate division-level, richness etc summaries
# 
# dt <- read_excel("Phyto compiled 2001-2017.xlsx", sheet = 1) # ignore warnings (poorly named columns)??
# dt <- dt[, c("STATION", "SPECCODE", "DIV", "SAMPLENUM", "SPECIES","CELLML", "BIOVOLUME")]
# dt <- subset(dt, SAMPLENUM != "SAMPLENUM") # full dataset has an aberrant entry SAMPLENUM as a value
# dt$CELLML <- ifelse(dt$SAMPLENUM == "14GH46I72", dt$CELLML/10, dt$CELLML) # an aberrant count for one sample
# 
# dt2 <- read_excel("PHYCNT18.xls", sheet = 1)
# dt2 <- dt2[, c("STATION", "SPECCODE", "DIV", "SAMPLENUM", "SPECIES","CELLML", "BIOVOLUME")]
# dt3 <- merge(dt, dt2, all = TRUE)
# 
# # new in this version, read in edibility data
# ed <- read.csv("ER_algalSppEdibility.csv")
# ed <- mutate_if(ed, is.character, as.factor)
# # ed2 <- subset(ed, SPECCODE == "CYCCOMES")
# ed <- ed[-132,] # remove one dup
# 
# dt3 <- join(dt3, ed[ c(1,6,7)], by = "SPECCODE") # join to phyto data, not ideal that it's not renamed but need for flow
# dt3 <- subset(dt3, nutrition.edibility != "l" & shape.edibility != "l") # subset to edible only
# 
# idPhyto1 <- read_excel("sample detailsUpdated.xlsx", sheet = "PHYCNT") # get ID info for phyto sample numbers to match to zoo
# idPhyto18 <- read_excel("PHYSAM18.xls", sheet = "PHYSAM") # Add sample details for 2018
# idPhyto18$YEAR <- 2018
# setnames(idPhyto18, old = c("SEASON", "DEPTH"), new = c("CRUISE", "SAMPLE TYPE"))
# idPhyto <- join(idPhyto1, idPhyto18[,c(5,4,1,19,23)], by = "SAMPLENUM", type = "full")
# 
# idPhyto$STATION <- gsub(" ", "", idPhyto$STATION, fixed = TRUE) # remove spaces in station IDs
# idPhyto$STATION <- sub("M$", "", idPhyto$STATION) # trailing M from master incompletely removed, remove M
# idPhyto$CRUISE <- recode(idPhyto$CRUISE, "SUMMER" = "Summer", "SPRING" = "Spring") # correct spelling for season
# idPhyto$id <- paste(idPhyto$YEAR, idPhyto$CRUISE, idPhyto$STATION, sep =".") # create lake*station*year ID to match to zoo
# 
# joinPhyto <- join(dt3[,-8:-9], idPhyto[,4:6], by = "SAMPLENUM") # this changed to remove edibility
# 
# # joinPhyto2 <- joinPhyto[complete.cases(joinPhyto[,8]),] #
# joinPhyto2 <- subset(joinPhyto, `SAMPLE TYPE` == "INT") # 79K obs
# joinPhyto2$DIV <- recode(joinPhyto2$DIV, "EUG" = "CHL") # merge CHL and EUG
# 
# phytoGr <- aggregate(BIOVOLUME ~ id + DIV, joinPhyto2[,c(3,7,9)],  FUN = function(x) sum = sum(x)) # group phyto by division
# phytoGr2 <- spread(phytoGr, DIV, BIOVOLUME) # wide format
# phytoGr2[is.na(phytoGr2)] <- 0
# 
# # Phyto diversity to zoo biovolume relationships
# phyto <- spread(joinPhyto2[, c(2,7,9)], SPECCODE, BIOVOLUME)
# phyto[is.na(phyto)] <- 0
# 
# phyto$totAlBiov <- rowSums(phyto[,-1], na.rm = TRUE) # calculate total algal biovolume per sample
# phyto$algalDiversityH <- diversity(phyto[,-c(1,446)]) # corrected this January 20, 2022, previously was counting totalBiov
# phyto$algalRichness <- specnumber(phyto[,-c(1,446)])
# # algalH <- cbind.data.frame(id = phyto[,1], algalDiversityH, algalRichness)
# 
# joinPhytoH <- join(phytoGr2, phyto[,c(1,446:448,2:445)], by = "id")
# # write.csv(joinPhytoH, "allPhytoSummariesEdible.csv", row.names = F) # 1978 obs, all phyto species biovolume data, DIV, diversity

#########################################
## Zoo data prep, step 1 in Steph's data
## joinZoo from Steph's modified workflow (00.2)
## NB from Steph: I updated the script and pushed it to the repo to include D20 and D100 dataframes 
## that are '*_by_condensed_SPECCODE.' These dfs have one SPECCODE/taxa, while the "*_by_SPECCODE" 
## have some SPECCODES that represent the same taxa (e.g. one for adult copepods and one for juveniles). 
## The *_by_condensed_SPECCODE df should be used for any ordination analysis.
# VisitTable_D20$Station <- sub("M$", "", VisitTable_D20$Station) # remove master designation to match phyto stations/strip last M
# VisitTable_D20$id <- as.factor(paste(VisitTable_D20$Year, VisitTable_D20$Season, VisitTable_D20$Station, sep = "."))
# joinZoo <- join(D20_by_condensed_SPECCODE[,-1], VisitTable_D20[,c(1,15)], by = "Visit_ID")
# joinZooSpread <- spread(joinZoo[,c(2,4,10)], SPECCODE, BIOMASS_ug_m3) # wide format
# joinZooSpread[is.na(joinZooSpread)] <- 0
# 
# # # compare to Old, note 36 more taxa!
# # joinZooOld <- join(D20_by_SPECCODE, VisitTable_D20[,c(1,15)], by = "Visit_ID")
# # joinZooOldSpread <- spread(joinZooOld[,c(4,6,9)], SPECCODE, BIOMASS_ug_m3) # wide format 156 taxa
# 
# # Zoo diversity and total biomass
# joinZooSpread$totBiomass <- rowSums(joinZooSpread[,-1], na.rm = TRUE) # calculate total zoo biomass per sample/m3
# joinZooSpread$zooH <- diversity(joinZooSpread[,-c(1,122)])
# joinZooSpread$zooRichness <- specnumber(joinZooSpread[,-c(1,122)])
#   
# # zoo biomass by major group
# zooGr <- aggregate(BIOMASS_ug_m3 ~ SF_group + id, joinZoo[,c(4,7,10)], FUN = function(x) sum = sum(x)) 
# zooGr2 <- spread(zooGr, SF_group, BIOMASS_ug_m3) # wide-format zoo (from 00.2)
# zooGr2[is.na(zooGr2)] <- 0
# joinZoo2 <- join(zooGr2, joinZooSpread, by = "id")
# setnames(joinZoo2, old = c("Pred. clad.", "Other clad."), new = c("predClad", "otherClad"))
# 
# # Zoo biomass by pred/prey
# joinZoo2$predRatio <- (joinZoo2$Limno+joinZoo2$predClad)/(joinZoo2$totBiomass - joinZoo2$Limno - joinZoo2$predClad)
# # write.csv(joinZoo2[,c(1:11, 132:135, 12:131)], "allZooSummaries.csv", row.names = F) # 2209 obs/all zoo spp/group/H/predator data

#########################################
## Group-level cross-correlations, richness, diversity, pred/prey, final figures
## NB: code changes for edible vs. total, pay attention
# ph <- read.csv("allPhytoSummaries.csv")
ph <- read.csv("allPhytoSummariesEdible.csv") # to re-run edibility figures/analyses, NB smry2 changes
## zo <- read.csv("allZooSummaries.csv") # old data with errors in veligers
zo <- read.csv("63Mesh_06June2023_from_v7_4_by_group.csv") # June 6 new zoo data with length issues corrected

zo$id <- as.factor(paste(zo$Year, zo$Season, zo$Station, sep = "."))
zo2 <- reshape2::dcast(zo[,c(7,11,12)], id ~ group, value.var = "BIOMASS_ug_m3")
zo2[is.na(zo2)] <- 0

zoTot <- aggregate(BIOMASS_ug_m3 ~ id, zo[,c(11,12)], FUN = function(x) sum = sum(x))
zo3 <- join(zo2, zoTot, by = "id")

zoNoVelig <- subset(zo, group != "Mollusc")
zoTotNoVelig <- aggregate(BIOMASS_ug_m3 ~ id, zoNoVelig[,c(11,12)], FUN = function(x) sum = sum(x))
setnames(zoTotNoVelig, old = "BIOMASS_ug_m3", new = "totBiomassNoVelig")
zo4 <- join(zo3, zoTotNoVelig, by = "id")

setnames(zo4, old = c("Pred. clad.", "BIOMASS_ug_m3"), new = c("predClad", "totBiomass"))
zo4$predRatio <- (zo4$Limno + zo4$predClad)/(zo4$totBiomassNoVelig - zo4$Limno - zo4$predClad)

# write.csv(zo4, "zooplanktonDiv_06June2023.csv", row.names = F)
zoph <- join(zo4, ph, by = "id", type = "inner") # 1451 obs
# write.csv(zoph, "zoph.csv", row.names = F)
smry <- zoph[,c(1:25)] # this changed from previous version
smry$season <- as.factor(stri_sub(smry$id,6,8))
smry$lake <- as.factor(stri_sub(smry$id,-4,-3))
smry$year <- as.numeric(stri_sub(smry$id,1,4))
  
smry1 <- subset(smry, totAlBiov > 100) # drop one outlier
# # x <- smry %>% slice_min(totAlBiov, n = 5)
# write.csv(smry1, "zooPhytoJoinedEdible.csv", row.names = F)

# smry2 <- smry1[,c(1,12:14,23:28)] # for all
smry2 <- smry1[,c(1,12:14,21:23,26:28)] # for edible
smry2$lake <- "All lakes"
# smry3 <- rbind(smry1[,c(1,12:14,23:28)], smry2) # All: this is to have all lakes in the same panel
smry3 <- rbind(smry1[,c(1,12:14,21:23,26:28)], smry2) # Edible: this is to have all lakes in the same panel

smry3$lake <- recode(smry3$lake, "ER" = "Erie", "HU" = "Huron", "MI" = "Michigan",
                     "ON" = "Ontario", "SU" = "Superior")

smry3$season <- recode(smry3$season, "Spr" = "spring", "Sum" = "summer")
# smry3$difVelig <- (smry3$totBiomass - smry3$totBiomassNoVelig)/smry3$totBiomass *100
# pdf("fig1_phytoZooBiov6facetNoVelig.pdf", width = 5.5, height = 4)
ggplot(smry3, aes(y = log10(totBiomassNoVelig), x = log10(totAlBiov), col = season)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.1) +
  facet_wrap(~ lake) +
  labs(x = expression(paste("Edible phytoplankton biovolume, ", µm^3, "/L, log")),
  # labs(x = expression(paste("Phytoplankton biovolume, ", µm^3, "/L, log")), 
                  y = expression(paste("Zooplankton biomass, ", µg/m^3, " log"))) +
  scale_colour_manual(values = c("green", "blue")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.position = c(0.92,0.63),
                     legend.title=element_blank(),
                     legend.text = element_text(size = 8),
                     legend.key.size = unit(0.8, 'cm'),
                     legend.key.height = unit(0.4, 'cm'),
                     legend.key.width = unit(0.4, 'cm'),
                     legend.background = element_blank())
# ggsave("Fig1new.eps", device = "pdf", width = 6, height = 4)
# ggsave("FigS3.eps", device = "pdf", width = 6, height = 4)
# dev.off()

smryX <- subset(smry1, season == "Sum")
fit <- lm(log10(totBiomassNoVelig) ~ algalDiversityH, smryX)
summary(fit)

# pdf("fig3_phytoHZooBiov6facetJune15.pdf", width = 5.5, height = 4)
ggplot(smry3, aes(y = log10(totBiomassNoVelig), x = algalDiversityH, col = season)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.1) +
  facet_wrap(~ lake) +
  labs(x = expression(paste("Edible phytoplankton Shannon diversity")), 
       y = expression(paste("Zooplankton biomass, ", µg/m^3, " log"))) +
  # xlab("Edible phytoplankton Shannon diversity") + 
  scale_colour_manual(values = c("green", "blue")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.position = c(0.92,0.63),
                     legend.title=element_blank(),
                     legend.text = element_text(size = 8),
                     legend.key.size = unit(0.8, 'cm'),
                     legend.key.height = unit(0.4, 'cm'),
                     legend.key.width = unit(0.4, 'cm'),
                     legend.background = element_blank())
# ggsave("Fig3new.eps", device = "pdf", width = 6, height = 4)
ggsave("FigS4new.eps", device = "pdf", width = 6, height = 4)
# dev.off()

## Decoupling figure
## Fig 2
# pdf("fig2June15.pdf", width = 5.5, height = 4)
ggplot(smry1, aes(y = log10(totBiomassNoVelig), x = log10(totAlBiov), col = season)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.1) +
  facet_wrap(~ year) +
  labs(x = expression(paste("Phytoplankton biovolume, ", µm^3, "/L, log")), 
       y = expression(paste("Zooplankton biomass, ", µg/m^3, " log"))) +
  scale_colour_manual(values = c("green", "blue")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.position = "NONE")
# ggsave("Fig2new.eps", device = "pdf", width = 6, height = 4)
# dev.off()


## Herbivorous zooplankton only figures, to address reviewer comments
smry1$herb <- smry1$totBiomassNoVelig - smry1$Limno - smry1$predClad
smry2 <- smry1[,c(1,12:14,21:23,26:28,29)] # for edible
smry2$lake <- "All lakes"
# smry3 <- rbind(smry1[,c(1,12:14,23:28)], smry2) # All: this is to have all lakes in the same panel
smry3 <- rbind(smry1[,c(1,12:14,21:23,26:28,29)], smry2) # Edible: this is to have all lakes in the same panel

smry3$lake <- recode(smry3$lake, "ER" = "Erie", "HU" = "Huron", "MI" = "Michigan",
                     "ON" = "Ontario", "SU" = "Superior")

smry3$season <- recode(smry3$season, "Spr" = "spring", "Sum" = "summer")

# pdf("ediblePhytoHerbZoo6facetJune15.pdf", width = 5.5, height = 4)
ggplot(smry3, aes(y = log10(herb), x = log10(totAlBiov), col = season)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
  facet_wrap(~ lake) +
  labs(x = expression(paste("Edible phytoplankton biovolume, ", µm^3, "/L, log")), 
       y = expression(paste("Herbivorous zooplankton biomass, ", µg/m^3, " log"))) +
  scale_colour_manual(values = c("green", "blue")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.position = "NONE")
ggsave("FigS5new.eps", device = "pdf", width = 6, height = 4)
# dev.off()
# 
# pdf("ediblePhytoDivHerbZoofacetJune15.pdf", width = 5.5, height = 4)
ggplot(smry3, aes(y = log10(herb), x = algalDiversityH, col = season)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
  facet_wrap(~ lake) +
  labs(x = expression(paste("Edible phytoplankton Shannon diversity")), 
       y = expression(paste("Herbivorous zooplankton biomass, ", µg/m^3, " log"))) +
  scale_colour_manual(values = c("green", "blue")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.position = "NONE")
ggsave("FigS6new.eps", device = "pdf", width = 6, height = 4)
# dev.off()

## pdf("figS1_phytoHZooDiv6facet.pdf", width = 5.5, height = 4) # did not re-do without veligers, see earlier version
# ggplot(smry3, aes(y = zooH, x = algalDiversityH, col = season)) +
#   geom_point(alpha = 0.3, size = 1) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ lake) +
#   xlab("Phytoplankton Shannon diversity") + ylab("Zooplankton Shannon diversity") +
#   scale_colour_manual(values = c("green", "blue")) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      legend.position = c(0.92,0.63),
#                      legend.title=element_blank(),
#                      legend.text = element_text(size = 8),
#                      legend.key.size = unit(0.8, 'cm'),
#                      legend.key.height = unit(0.4, 'cm'),
#                      legend.key.width = unit(0.4, 'cm'),
#                      legend.background = element_blank())
# # dev.off()

# ## Predator-prey ratios
## NB: mirror in previous version, do not use smry3 (which doubles up data with all lakes for faceting)
smry4 <- reshape2::melt(smry1[,c(1,14,23:26)], variable = "metric", id.vars = c("id", "season", "predRatio"))
# x <- smry4[,-4:-5] %>% slice_max(predRatio, n = 50)
# x <- unique(x)
smry4$metric <- recode(smry4$metric, "totAlBiov" = "Algal biovolume, log",
                       "algalDiversityH" = "Algal diversity",
                       "algalRichness" = "Algal richness")
smry4$value <- ifelse(smry4$metric == "Algal biovolume, log", log10(smry4$value), smry4$value)

# pdf("figS2_predPreyJune15.pdf", width = 5.5, height = 4) # Not with Edible!!!!!
ggplot(smry4, aes(y = predRatio, x = value)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "gam", se = TRUE, size = 0.8, alpha = 0.3) +
  facet_wrap(~ metric, scales = "free") +
  # facet_grid(season ~ metric, scales = "free") +
  xlab("Phytoplankton metric") + ylab("Zooplanktivore-grazer ratio") +
  ylim(0, 1.8) +
  # scale_colour_manual(values = c("green", "blue")) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.position = "NONE",
                     legend.title=element_blank(),
                     legend.text = element_text(size = 8),
                     legend.key.size = unit(0.8, 'cm'),
                     legend.key.height = unit(0.4, 'cm'),
                     legend.key.width = unit(0.4, 'cm'),
                     legend.background = element_blank())
ggsave("FigS2new.eps", device = "pdf", width = 6, height = 4) # redoing this because reviewer had trouble opening
# dev.off()

# pdf("figS3_phytoEdZooBiov6facet.pdf", width = 5.5, height = 4)
# ggplot(smry3, aes(y = log10(totBiomass), x = log10(totAlBiov), col = season)) +
#   geom_point(alpha = 0.3, size = 1) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ lake) +
#   xlab("Edible phytoplankton biovolume, µm3/L, log") + ylab("Zooplankton biomass, µg/m3, log") +
#   scale_colour_manual(values = c("green", "blue")) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), 
#                      legend.position = c(0.92,0.63),
#                      legend.title=element_blank(),
#                      legend.text = element_text(size = 8),
#                      legend.key.size = unit(0.8, 'cm'),
#                      legend.key.height = unit(0.4, 'cm'),
#                      legend.key.width = unit(0.4, 'cm'),
#                      legend.background = element_blank())
# dev.off()
# 
# pdf("figS4_phytoHEdZooBiov6facet.pdf", width = 5.5, height = 4)
# ggplot(smry3, aes(y = log10(totBiomass), x = algalDiversityH, col = season)) +
#   geom_point(alpha = 0.3, size = 1) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ lake) +
#   xlab("Edible phytoplankton Shannon diversity") + ylab("Zooplankton biomass, µg/m3, log") +
#   scale_colour_manual(values = c("green", "blue")) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      legend.position = c(0.92,0.63),
#                      legend.title=element_blank(),
#                      legend.text = element_text(size = 8),
#                      legend.key.size = unit(0.8, 'cm'),
#                      legend.key.height = unit(0.4, 'cm'),
#                      legend.key.width = unit(0.4, 'cm'),
#                      legend.background = element_blank())
# dev.off()


## Figure 4
install_load(c("ggExtra", "ggcorrplot"))

# smry2 <- subset(smry, season != "Spr")
smry1$BAP <- smry1$BAC + smry1$BAP
smry2 <- smry1[,c(13,11,2:6,10,23:25,16,17,20)] # all; CRY, CHR???!!! BAC?!!!
# smry2 <- smry1[,c(12:13,11,2:6,10,24:26,17:18,21)] # edible

setnames(smry2, old = c("algalRichness", "algalDiversityH", "totAlBiov" , "CYA", "CHL", "BAP", "totBiomassNoVelig", 
                        "predClad", "Limno"),
              new = c("Algal richness", "Algal diversity", "Algal biovolume","Bluegreen", "Green",
                      "Diatom", "Zooplankton biomass", "Predatory Cladocerans", 
                      "Limnocalanus"))

cols <- c("Algal biovolume","Bluegreen", "Green",
          "Diatom", "Daphnid", "Cyclopoid", "Calanoid", "Bosminid", "Rotifers", "Zooplankton biomass", 
          "Predatory Cladocerans", "Limnocalanus")
smry2[cols] <- log10(smry2[cols]+1)

corr <- cor(smry2, use = "complete.obs", method = "spearman") #method = "spearman" [,c(3,6:9,5,4,11,12,15,16,1)]

pdf("fig4_crossCorr_June16.pdf", width = 6, height = 6) # did not redo in Sept 2023
ggcorrplot(corr, hc.order = F,
           type = "lower",
           lab = TRUE,
           lab_col = "white",
           lab_size = 2.5,
           # method="circle",
           insig = "blank", outline.col = "white",
           colors = c("dark green", "white", "dark blue"), # NB this takes numeric and verbal arguments
           ggtheme = theme_bw) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.2,0.75))
ggsave("Fig4.eps", device = "pdf", width = 6, height = 6)
# dev.off()

## Generate p-values to confirm
# p.mat <- cor_pmat(smry2, use = "complete.obs", method = "spearman")
# write.csv(p.mat, "p.mat.csv")

#########################################
# # RF for zoo predRatio
# ph <- read.csv("allPhytoSummaries.csv")
# zo <- read.csv("allZooSummaries.csv")
# 
# zoph <- join(zo, ph, by = "id", type = "inner") # 1451 obs
# 
# rf <- zoph[,-16:-135] # remove zoo spp data
# rf <- rf[, colSums(rf !=0) > 145] # remove rare spp
# # write.csv(rf, "rf.csv", row.names = F)
# # rf of pred cladocera and Limno omitting predClad, Limno, all zoo spp
# 
# fit <- randomForest(rf$Limno ~ ., dat = rf[,-c(1,15)], ntree = 10000, keep.forest = TRUE,
#                     importance = TRUE, na.action = na.omit)
# 
# fit
# 
# # Model optimization by dropping the worst predictors
# imp <- as.data.frame(importance(fit)) # if dimension error restart
# imp <- cbind(names = rownames(imp), imp)
# setnames(imp, old = "%IncMSE", new = "perMse")
# 
# # Double transpose
# rfT <- as.data.frame(t(rf[,-1])) # do not transpose with colnames in it
# colnames(rfT) <- rf$id
# rfT$names <- factor(row.names(rfT))
# rfT2 <- rfT[-c(5,14),c(1452,1:1451)] # drop pred cladocerans and limnocalanus
# 
# rfBest <- join(imp[,-3], rfT2, by = "names")
# # compare(rfT$names, imp$names)
# 
# rfBest2 <- slice_max(rfBest, perMse, n = 28) # keep only the best predictors
# rfBest2T <- as.data.frame(t(rfBest2[,-1:-2]))
# colnames(rfBest2T) <- rfBest2$names
# rfBest2T$id <- rownames(rfBest2T)
# rfBest2T2 <- join(rfBest2T,rf[c(1,6)], by = "id")
# 
# fit2 <- randomForest(rfBest2T2$Limno ~ ., dat = rfBest2T2[,-c(29)], ntree = 10000, keep.forest = TRUE,
#                     importance = TRUE, na.action = na.omit)
# 
# fit2
# saveRDS(fit, "limnoAllRf.rds")
# saveRDS(fit2, "limnoBestRf.rds")

# pdf("fig4b_limnoBestRf.pdf", width = 4.5, height = 6)
# varImpPlot(fit2, n.var = 15, type = 1, main = "Limnocalanus variation explained 36%")
# dev.off()
# 
# plot(log10(zoph$Limno+1) ~ log10(zoph$BITOLLU+1))
# mod1 <- lm(log10(zoph$Limno+1) ~ log10(zoph$SYNOSTE+1))
# summary(mod1)
# 
# fit3 <- randomForest(rf$predClad ~ ., dat = rf[,-c(1,15)], ntree = 10000, keep.forest = TRUE,
#                     importance = TRUE, na.action = na.omit)
# 
# fit3
# 
# pdf("fig4a_predCladRf.pdf", width = 4.5, height = 6)
# varImpPlot(fit3, n.var = 15, type = 1, main = "Predatory cladoceran variation explained 50%")
# dev.off()
# 
# fit4 <- randomForest(rf$Rotifers ~ ., dat = rf[,-c(1,15)], ntree = 10000, keep.forest = TRUE,
#                      importance = TRUE, na.action = na.omit)
# 
# fit4
# 
# pdf("fig4c_rotiferaRf.pdf", width = 4.5, height = 6)
# varImpPlot(fit4, n.var = 15, type = 1, main = "Rotifer variation explained 26%")
# dev.off()
# 
# plot(log10(zoph$Rotifers+1) ~ log10(zoph$predClad+1))
# mod1 <- lm(log10(zoph$Limno+1) ~ log10(zoph$SYNOSTE+1))
# summary(mod1)

#########################################
# # Filtering to edible algae only, Euan ranked externally; all workflow changed to add it in
# ph <- read.csv("allPhytoSummaries.csv")
# ph <- ph[,-c(2:9,11:12)]
# 
# ph2 <- reshape2::melt(ph, value.name = "biov", id.vars = c('id', 'totAlBiov'))
# ph2$maxRelBiov <- ph2$biov/ph2$totAlBiov*100
# 
# # ph3 <- ph2[,3:5] %>%
# #   group_by(variable)
# #   summarize_each(funs(mean(., na.rm=T), n = sum(!is.na(.)), sd = sd(., na.rm = T)))
# 
# # ph3 <- aggregate(. ~ variable, ph2[,-c(1,2,4)], 
# #                  FUN = function(x) c(max = max(x), freq = sum(x != 0), na.rm = TRUE, na.action = NULL))
#   
# ph3 <- aggregate(. ~ variable, ph2[,-c(1,2,4)], FUN = max, na.rm = TRUE, na.action = NULL)
# ph4 <- aggregate(. ~ variable, ph2[,-c(1,2,4)], function(x) sum(x != 0))
# setnames(ph4, old = "maxRelBiov", new = "frequency")
# 
# ph5 <- join(ph3, ph4)
# setnames(ph5, old = "variable", new = "SPECCODE")
# 
# dt <- read_excel("Phyto compiled 2001-2017.xlsx", sheet = 1) # ignore warnings (poorly named columns)??
# dt <- dt[, c("STATION", "SPECCODE", "DIV", "SAMPLENUM", "SPECIES","CELLML", "BIOVOLUME")]
# dt <- subset(dt, SAMPLENUM != "SAMPLENUM") # full dataset has an aberrant entry SAMPLENUM as a value
# dt$CELLML <- ifelse(dt$SAMPLENUM == "14GH46I72", dt$CELLML/10, dt$CELLML) # an aberrant count for one sample
# 
# dt2 <- read_excel("PHYCNT18.xls", sheet = 1)
# dt2 <- dt2[, c("STATION", "SPECCODE", "DIV", "SAMPLENUM", "SPECIES","CELLML", "BIOVOLUME")]
# dt3 <- merge(dt, dt2, all = TRUE)
# 
# idPhyto1 <- read_excel("sample detailsUpdated.xlsx", sheet = "PHYCNT") # get ID info for phyto sample numbers to match to zoo
# idPhyto18 <- read_excel("PHYSAM18.xls", sheet = "PHYSAM") # Add sample details for 2018
# idPhyto18$YEAR <- 2018
# setnames(idPhyto18, old = c("SEASON", "DEPTH"), new = c("CRUISE", "SAMPLE TYPE"))
# idPhyto <- join(idPhyto1, idPhyto18[,c(5,4,1,19,23)], by = "SAMPLENUM", type = "full")
# 
# idPhyto$STATION <- gsub(" ", "", idPhyto$STATION, fixed = TRUE) # remove spaces in station IDs
# idPhyto$STATION <- sub("M$", "", idPhyto$STATION) # trailing M from master incompletely removed, remove M
# idPhyto$CRUISE <- recode(idPhyto$CRUISE, "SUMMER" = "Summer", "SPRING" = "Spring") # correct spelling for season
# idPhyto$id <- paste(idPhyto$YEAR, idPhyto$CRUISE, idPhyto$STATION, sep =".") # create lake*station*year ID to match to zoo
# 
# joinPhyto <- join(dt3, idPhyto[,4:6], by = "SAMPLENUM")
# uniquePhyto <- unique(joinPhyto[,c(2,3,5)])
# 
# ph6 <- join(ph5, uniquePhyto)
# # write.csv(ph6, "algalSppEdibility.csv", row.names = F)

#########################################
# Which sites lost the most biovolume/time of day analyses Oct 2022
# Adding edible/inedible ratio analyses Jan 2023
ph <- read.csv("allPhytoSummaries.csv")
phEd <- read.csv("allPhytoSummariesEdible.csv")
setnames(phEd, old = "totAlBiov", new = "totAlBiovEd")

al <- join(ph[,1:10], phEd[,1:10], by = "id") # 1451 obs
zo <- read.csv("allZooSummaries.csv")

# # Adding time of day info (added Oct 19 2022)
# day <- read.csv("GLNPO_VisitTable_tidy.txt")
# day$id <- paste(day$Year, day$Season, day$Station, sep = ".")
# day$id <- sub("M$", "", day$id) # how the fuck does this work?
# zo2 <- join(zo, day, by = "id")
# 
# dif <- join(zo2[,c(1,12,145)], al, by = "id", type = "inner") #

dif <- join(zo[,c(1,12:15)], al, by = "id", type = "inner") #
# dif$difference <- dif$totAlBiov - dif$totAlBiovEd

dif$ratio <- dif$totAlBiovEd/dif$totAlBiov*100

dif$season <- as.factor(stri_sub(dif$id,6,8))
dif$lake <- as.factor(stri_sub(dif$id,-4,-3))
dif$year <- as.numeric(stri_sub(dif$id,1,4))

smry1 <- subset(dif, totAlBiov > 100)

# ggplot(smry1, aes(x = log10(difference))) +
#   geom_histogram()

# smry1$inedible <- smry1$totAlBiov - smry1$totAlBiovEd

# smry2 <- subset(smry1, year > 2011)
# smry3 <- subset(smry1, lake == "ER")
# pdf("alBiovEdibleInedibleSeason.pdf", width = 5.5, height = 4)
# ggplot(smry2, aes(y = log10(totAlBiovEd), x = log10(inedible), col = season)) +
#   # geom_density_2d_filled() +
#   geom_point(alpha = 0.3, size = 2) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ lake) +
#   # geom_abline(intercept = 0, slope = 1, size = 0.1) +
#   scale_colour_manual(values = c("green", "blue")) +
#   xlab("Inedible algal biovolume, µm3/L, log") + ylab("Edible algal biovolume, µm3/L, log") +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = "NONE")
# dev.off()

# pdf("biomassEdibleRatioSeasonLake.pdf", width = 5.5, height = 4)
# ggplot(smry1, aes(y = log10(totBiomass), x = ratio, col = season)) +
#   geom_point(alpha = 0.3, size = 2) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ lake) +
#   # geom_abline(intercept = 0, slope = 1, size = 0.1) +
#   scale_colour_manual(values = c("green", "blue")) +
#   xlab("Edible to total biovolume, %") + ylab("Zooplankton biomass, µg/m3, log") +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = "NONE")
# dev.off()
# 
# # complete cases!!!!
# pdf("alBiovEdibleTotalDay.pdf", width = 8, height = 6)
# ggplot(ph, aes(x = algalDiversityH, y = log10(totAlBiov))) +
#   # geom_density_2d_filled() +
#   geom_point(alpha = 0.7, size = 1.5) +
#   geom_smooth(method = "lm", se = F, size = 0.8, alpha = 0.6) +
#   # facet_grid(lake ~ season, scales = "free") +
#   # facet_wrap(~ year) +
#   scale_colour_manual(values = c("green", "blue", "red")) +
#   scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
#   xlab("Edible algal biovolume, µm3/L, log") + ylab("Total zooplankton biomass, µg/m3, log") +
#   labs(col = "") +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = c(0.9,0.13))
# dev.off()

## FigS7. Time of day
# smry1 <- smry1[complete.cases(smry1[,"Sun_category"]),]
dat <- join(smry1[,c(1,13,21,26,27)], unique(zo[,c(8,12)]), by = "id")
dat <- dat[complete.cases(dat[,"Sun_category"]),]

dat$lake <- recode(dat$lake, "ER" = "Erie", "HU" = "Huron", "MI" = "Michigan",
                     "ON" = "Ontario", "SU" = "Superior")

dat$season <- recode(dat$season, "Spr" = "spring", "Sum" = "summer")

# pdf("zooAlBiovEdibleTimeOfDayJune16.pdf", width = 6.5, height = 4)
ggplot(dat, aes(y = log10(totBiomassNoVelig), x = log10(totAlBiov), col = Sun_category)) +
  # geom_density_2d_filled() +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = F, size = 0.8, alpha = 0.6) +
  facet_grid(lake ~ season, scales = "free") +
  # facet_wrap(~ Sun_category) +
  scale_colour_manual(values = c("grey", "orange", "black")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.5)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.5)) +
  labs(x = expression(paste("Edible phytoplankton biovolume, ", µm^3, "/L, log")),
       y = expression(paste("Zooplankton biomass, ", µg/m^3, " log"))) +
  labs(col = "") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.position = "right")
ggsave("FigS7new.eps", device = "pdf", width = 6, height = 4)
# dev.off()

x <- subset(smry1, Sun_category == "Night")
fit <- lm(log10(totAlBiovEd) ~ log10(totAlBiov), x) # night R2 = 0.73; day R2= 0.76; dawn R2 =0.81
summary(fit)

#########################################
## Old figures
## Fig. S3.
# pdf("figS3_phytoZooBiovEdible3.pdf", width = 5.5, height = 4)
# ggplot(smry1, aes(y = log10(totBiomass), x = log10(totAlBiov), col = season)) +
#   # geom_density_2d_filled() +
#   geom_point(alpha = 0.6, size = 1.5) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   # facet_wrap(~ season) +
#   xlab("Edible algal biovolume, µm3/L, log") + ylab("Total zooplankton biomass, µg/m3, log") +
#   scale_colour_manual(values = c("green", "blue")) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = "NONE")
# dev.off()
# 
# smrySum <- subset(smry, season == "Sum")
# fit <- lm(log10(totBiomass) ~ log10(totAlBiov), smrySum) # R2 = 0.2182; 0.08683 for edible
# summary(fit)
# 
# smrySpr <- subset(smry, season == "Spr")
# fit <- lm(log10(totBiomass) ~ log10(totAlBiov), smrySpr) # R2 = 0.06494; 0.01816 for edible
# summary(fit)
# 
# fitAll <- lm(log10(totBiomass) ~ log10(totAlBiov), smry) # R2 = 0.2116
# summary(fitAll)
# 
# # Fig. S2
# smry1 <- subset(smry, totAlBiov > 100)
# # x <- smry %>% slice_min(totAlBiov, n = 5)
# pdf("figS2_phytoZooBiovYrEdible.pdf", width = 7, height = 5)
# ggplot(smry1, aes(y = log10(totBiomass), x = log10(totAlBiov), col = season)) +
#   geom_point(alpha = 0.6, size = 1.5) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ year) +
#   xlab("Edible algal biovolume, µm3/L, log") + ylab("Total zooplankton biomass, µg/m3, log") +
#   scale_colour_manual(values = c("green", "blue")) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = "NONE")
# dev.off()
# 
# Fig. S1
# x <- smry %>% slice_min(totAlBiov, n = 5)
# pdf("figS1_phytoZooBiov_lk_setScales.pdf", width = 7, height = 5)
# ggplot(smry1, aes(y = log10(totBiomass), x = log10(totAlBiov), col = season)) +
#   geom_point(alpha = 0.6, size = 1.5) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ lake) +
#   xlab("Algal biovolume, µm3/L, log") + ylab("Total zooplankton biomass, µg/m3, log") +
#   scale_colour_manual(values = c("green", "blue")) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = c(0.8,0.2))
# dev.off()
# 
# # Fig. 2
# pdf("fig2_algalDivZooBiomLakeEdible.pdf", width = 5.5, height = 4)
# ggplot(smry, aes(y = log10(totBiomass), x = algalDiversityH, col = season)) +
#   geom_point(alpha = 0.6, size = 1.5) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ lake) +
#   scale_colour_manual(values = c("green", "blue")) +
#   xlab("Edible algal Shannon diversity") + ylab("Total zooplankton biomass, µg/m3, log") +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = c(0.8, 0.2))
# dev.off()
# 
# smrySpr <- subset(smry, season == "Spr")
# fit <- lm(log10(totBiomass) ~ algalDiversityH, smrySpr) # R2 = 0.03929
# summary(fit)
# 
# fitAll <- lm(log10(totBiomass) ~ algalDiversityH, smry) # R2 = 0.006317
# summary(fitAll)
# 
# # Plots not included in main text or appendices:
# # Richness vs richness, by lake
# pdf("zoPhRichEdible3.pdf", width = 5.5, height = 4)
# ggplot(smry, aes(y = zooRichness, x = algalRichness)) +
#   # geom_point(alpha = 0.6, size = 1.5) +
#   geom_density_2d_filled() +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ season, scales = "free") +
#   scale_colour_manual(values = c("green", "blue")) +
#   xlab("Edible algal richness") + ylab("Zooplankton richness") +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = "NONE")
# dev.off()
# 
# pdf("zoPhHEdible.pdf", width = 5.5, height = 4)
# ggplot(smry, aes(y = zooH, x = algalDiversityH, col = season)) +
#   geom_point(alpha = 0.6, size = 1.5) +
#   geom_smooth(method = "lm", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ lake, scales = "free") +
#   scale_colour_manual(values = c("green", "blue")) +
#   xlab("Edible algal Shannon diversity") + ylab("Zooplankton Shannon diversity") +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = "NONE")
# dev.off()

# fit <- lm(zooRichness ~ algalRichness, smry) # R2 = 0.02569 
# summary(fit)
# 
# fit <- lm(zooH ~ algalDiversityH, smry) # p-value: 0.6954
# summary(fit)

# pdf("zoPhH.pdf", width = 5.5, height = 4)
# ggplot(smry, aes(y = predRatio*100, x = algalDiversityH, col = season)) +
#   geom_point(alpha = 0.6, size = 1.5) +
#   geom_smooth(method = "gam", se = TRUE, size = 0.8, alpha = 0.1) +
#   facet_wrap(~ lake, scales = "free") +
#   scale_colour_manual(values = c("green", "blue")) +
#   xlab("Algal Shannon diversity") + ylab("Zooplankton Shannon diversity") +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), legend.position = "NONE")
# dev.off()

# fit <- lm(log10(predRatio+1) ~ algalDiversityH, smry) # 0.04368 
# summary(fit)
# 
# fit <- lm(log10(predRatio+1) ~ algalRichness, smry) # 0.03363 
# summary(fit)

## Site info
smry <- zoph # (see above)
smry$st <- stri_sub(smry$id, 1,-4)
smry$season <- as.factor(stri_sub(smry$id,6,8))
smry$lake <- as.factor(stri_sub(smry$id,-4,-3))
smry$year <- as.numeric(stri_sub(smry$id,1,4))

smry2 <- smry[,c(700:703)] %>% 
  group_by(year, lake, season) %>%
  count(st)

write.csv(smry2, "SI Table.csv")

## Figures Fig5 and S8 were redone directly in Adobe to change ug3 to superscript (Sept 28 2023) to address reviewer comments
