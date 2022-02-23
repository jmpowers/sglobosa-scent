library(tidyverse)
library(reshape2)
library(lubridate)
library(hms)
library(vegan)
library(ggvegan)
library(viridis)

#setwd("~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/Schiedea/sglobosa-scent")
source("./read_shimadzu.R")

# Read chromatograms ------------------------------------------------------

#reading in the Shimadzu search output is slow, skip it
#sch.data <- read.shimadzu("schiedea190801.txt")
#sch.data2 <- read.shimadzu("Schiedea_190814.txt")
#sch.data <- rbind(sch.data, sch.data2)
#save(sch.data, file="schiedea_190801_190814.Rdata")
load("./data/schiedea_190801_190814.Rdata")

sch.all <- dcast(sch.data, Filename~Name, sum, value.var="Area")
rownames(sch.all) <- sch.all[,1]
sch.all[,1] <- NULL
sch.cut <- sch.all[,colSums(sch.all)>4e6]#arbitrary cutoff

#which names have which CAS numbers?
sch.data %>% filter(CAS.. != "0 - 00 - 0", Name!="") %>% group_by(Name) %>% summarize(CAS=first(CAS..)) %>% 
  mutate(CAS = str_remove_all(CAS, fixed(" "))) %>%  write_csv("./data/schiedea_Name_CAS.csv")

# k-means  --------------------------------------------------------------------
k <- 40
set.seed(1)
km <- kmeans(decostand(sch.cut, "log"), k, nstart=3)
#save(km, file="km30.Rdata")

sch.km <- tibble(FileName=rownames(sch.all)) %>% 
  mutate(rowSum = rowSums(sch.all),
         Project = str_extract(FileName, "AMBI|Bakeout|Blank|GLOB|KAHO|SCH|STEL") %>% replace_na("sample"),
         Type = fct_collapse(Project, blank=c("Blank","Bakeout"), air="AMBI", other_level = "sample"),
         nameBlank = Type=="blank",
         runYear = str_extract(FileName, "2018|2019|2020")  %>% factor,
         Cluster = km$cluster) %>% # Figure out which k-means clusters are the blanks
  mutate(kBlank = Cluster %in% (count(., nameBlank, Cluster) %>% filter(nameBlank, n>2) %>% pull(Cluster)),
         Mixup = nameBlank != kBlank)

with(sch.km, table(kBlank, nameBlank))

# Schiedea runs -----------------------------------------------------------
#Section moved from  ../Inventory/markes_sequence.R, code above that section outputs this Rdata:
load("./data/markes_sequence.rda")

Powers <- str_detect(sequ.summary$FullName, "Powers RMBL Data")
schiedea.batchids <- sequ.summary %>% filter(Powers) %>% select(id) %>% unique() %>% na.omit()
schgc <- sequ.summary %>% filter(id %in% schiedea.batchids$id | Powers) %>% #get entire batch if it had a sample that matches
  left_join(sch.km %>% select(FileName, nameBlank, Mixup, kBlank, Cluster)) %>% 
  left_join(sch.data %>% rename(FileName=Filename) %>% group_by(FileName) %>% tally(name="n_peaks")) %>% 
  mutate(vial = str_match(FileName, "_V(\\d{1,3})")[,2],
         dupe = !is.na(vial) & vial %in% vial[duplicated(vial)]) %>% 
  left_join(read_tsv("./data/leak_notebook.csv"), by=c("vial"="Vial"), na_matches = "never") %>% 
  mutate(TubeMatch = Tube == SteelTube) %>% 
  write_csv("./data/schiedea_all190815.csv") %>% #for annotation with verdicts
  bind_cols(read_csv("./data/Schiedea RMBL GC-MS Inventory - schiedea_all.csv") %>%  #read_sheet("1_rNiQ3IwKIKQhfxZTBJbKU14ed6L1rGQpe7e_ttJl6Y")
              select(sample, verdict)) %>% #read back the verdicts
  mutate(FileName2 = ifelse(is.na(sample), FileName, sample), #new filename from verdicts
         vial2 = as.integer(str_match(FileName2, "_V(\\d{1,3})")[,2]),
         dupe2 = !is.na(vial2) & vial2 %in% vial2[duplicated(vial2)])
rownames(schgc) <- 1:nrow(schgc)

#match up with metadata
s.nr <- read_csv("./data/Schiedea Volatiles Sampling - Samples.csv")%>% 
  filter(GC=="NR") %>% 
  mutate(Vial = as.integer(Vial),
         DN = factor(ifelse(Start > 16*60*60, "Night", "Day")))
sv.all <- full_join(s.nr, schgc, by=c("Vial"="vial2")) %>% 
  write_csv("./data/s.nr.schgc190815.csv")

#compare metadata (s.nr) and filenames (schgc)
sort(union(s.nr$Vial, schgc$vial2))
sort(setdiff(s.nr$Vial, schgc$vial2)) #not run yet
sort(setdiff(schgc$vial2, s.nr$Vial)) #run but not entered  
broken <- c(52,53,56,86,117,124,125,129,174) #broken traps - not run
empty <- c(7,16,18,28,40,313,314,315,316) # presumed clean - run
notrun <- c(5,8,21,30,33,50,51,55,57,64,66,74,75,80,81,82,83,84,88,89,91,95,100,102,103,104,107,108,109,
            11,12,22,23,24,25,26,58,59,60,61,62,63,58,92,94,96,97,98,9, 90, 79, 87, 72, 99)
setdiff(setdiff(setdiff(s.nr$Vial, schgc$vial2), broken), notrun)

# NMDS of blanks and all samples --------------------------------------------------------------------
nmds.sch <- metaMDS(decostand(sch.cut, "hellinger"), dist="bray", autotransform = FALSE, trymax=1, try=1)
nmds.points <- fortify(nmds.sch) %>% as_tibble() %>% 
  filter(Score=="sites") %>% left_join(sch.km, by=c("Label"="FileName"))

ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=log(rowSum))) + 
  geom_point() + scale_color_viridis_c() + theme_dark()

ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=Cluster, shape=Type)) + geom_point() +
  scale_color_gradientn(colors=turbo(k)) + theme_dark()

ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=nameBlank, shape=Type)) + geom_point()

ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, label=Cluster, color=Type)) + geom_text(size=3)

# CAP of blanks vs all samples ---------------------------------------------------------------------

sch.cap.kblank <- capscale(sch.cut ~ kBlank , distance="bray", metaMDSdist = F, data=sch.km)
sch.cap.kblank.points <-  fortify(sch.cap.kblank) %>% as_tibble() %>% 
  filter(Score=="sites") %>% left_join(sch.km, by=c("Label"="FileName"))
ggplot(sch.cap.kblank.points, aes(x=CAP1, y=MDS2, color=nameBlank)) + geom_point()
enframe(sch.cap.kblank$CCA$v) %>% arrange(value) # volatiles higher in samples compared to blanks


# Get metadata -----------------------------------------------------------------------
sv.full <- tibble(FileName = rownames(sch.cut)) %>%  ##2020-12-03: changed from FileName2 (renamed) to FileName (original)
  left_join(distinct(sv.all, FileName, .keep_all=T)) %>% #TODO adds two rows - a duplicate blank and a vial contaminated with two samples
  mutate(badtube = verdict %in% c("alreadyrun", "empty","leak-blank", "notmine", "notrun", "skip-blank", "skip-notrun","check-chrom","trapblank"),
         type = ifelse(str_detect(FileName2, "Blank|Bakeout") ,"blank", ifelse(Species=="ambient", "ambient", "floral")) %>% replace_na("other"),
         Species = ifelse(type=="blank", "blank", Species),
         sp = toupper(substr(Species, 1, 4)),
         across(c(type, Species, Population, Plant, Cutting, Sex, sp), factor),
         duped = duplicated(FileName))

ggplot(left_join(nmds.points,sv.full), aes(x=NMDS1, y=NMDS2, color=type, shape=badtube)) + geom_point()
ggplot(left_join(nmds.points,sv.full), aes(x=NMDS1, y=NMDS2, label=Cluster, color=paste(type,badtube))) + geom_text(size=3) #shows remaining inconsistencies

sch.cut2 <- sch.cut[!sv.full$badtube & sv.full$type !="other",]
sv <- sv.full %>% filter(!badtube & type !="other") %>% 
  mutate(type = fct_drop(type),
         sp2 = sp %>% replace(Population == "KK", "KAAL") %>% replace(Population == "HH", "HOOK") %>% droplevels(),
         Species2 = Species %>% replace(Population == "KK", "kaalae") %>% replace(Population == "HH", "hookeri"),
         Species2Pop = factor(paste(Species2, Population))) #%>% 
  #write_csv("./data/sv_190815_fixed.csv")

# Inventory
with(sv, table(kBlank, type))
with(sv, table(Species2, DN))
with(sv, table(Species2Pop, DN))
with(sv, table(Species2Pop, paste(Sex,DN)))

# NMDS with metadata -----------------------------------------------------------------------
spcol <- unique(c('#3366cc', '#dc3912', '#ff9900', '#109618', '#990099', '#0099c6', '#dd4477', '#66aa00', '#b82e2e', '#316395', '#994499', '#22aa99', '#aaaa11', '#6633cc', '#e67300', '#8b0707', '#651067', '#329262', '#5574a6', '#3b3eac', '#b77322', '#16d620', '#b91383', '#f4359e', '#9c5935', '#a9c413', '#2a778d', '#668d1c', '#bea413', '#0c5922', '#743411', '#3366cc'))

set.seed(1)
nmds.sch.cut2 <- metaMDS(sqrt(sch.cut2), dist="bray", autotransform = FALSE, trymax=5, try=5)
nmds.sch.cut2.plot <- ordiplot(nmds.sch.cut2, type = "n")#, xlim=c(0.6,1), ylim=c(-.2,0.2)
points(nmds.sch.cut2, display="sites", col=ifelse(is.na(sv$Species),"black",spcol[as.integer(sv$Species)]), pch=as.integer(sv$type)) #pch=as.integer(sv$DN)+15
text(nmds.sch.cut2, display="sites", col=ifelse(is.na(sv$Species),"black",spcol[as.integer(sv$Species)]), labels=sv$FileName)
ordiellipse(nmds.sch.cut2, sv$type, col=spcol, lwd=4)
legend("topleft", levels(sv$Species), fill=spcol, cex=0.7)
identify(nmds.sch.cut2.plot, what="sites", labels=rownames(sch.cut2), cex=0.7)

# CAP - species ---------------------------------------------------------------------
sch.cap <- capscale(sqrt(sch.cut2) ~ sv$Species, distance="bray", metaMDSdist = F)
sch.cap.plot <- plot(sch.cap, type="n")
points(sch.cap, display="sites", col=spcol[as.integer(sv$Species)], pch=as.integer(sv$type)) #pch=as.integer(sv$DN)+15
ordiellipse(sch.cap, na.omit(sv$Species), col=spcol, lwd=4)
legend("topleft", levels(sv$Species), fill=spcol, cex=0.7)
identify(sch.cap.plot, what="sites", labels=rownames(sch.cut2[!is.na(sv$Species),]), cex=0.7)
#View(sch.cap$CCA$v)

# Filtering ---------------------------------------------------------------
library(bouquet)
sch.data.cut2 <- sch.data[sch.data$Filename %in% sv$FileName,]
sch.data.cut2$Filename <- sv$FileName2[match(sch.data.cut2$Filename, sv$FileName)]
sch.data.cut2 <- sch.data.cut2[sch.data.cut2$Name != "", ]
sch.data.cut2$Name <- droplevels(sch.data.cut2$Name)
longdata <- load_longdata(sch.data.cut2, sample="Filename", RT="Ret.Time", name="Name", area="Area", match = "SI", maxmatch=100)
metadata <- load_metadata(sv, date="SampleDate", sample="FileName2", group="Species", type="type", amount="Flrs")
vol.all <- make_sampletable(longdata)
chems <- make_chemtable(longdata, metadata)
chemsf <- filter_RT(chems, 2, 17) %>% 
  filter_match(0.8) %>% 
  filter_freq(0.1, group = FALSE) %>% 
  filter_contaminant(cont.list = c("Caprolactam","2-Undecanone, 6,10-dimethyl-","Octanoic acid, 2-tetrahydrofurylmethyl ester","Nonadecane, 9-methyl-","Decyl pentyl ether","Dodecane, 2,6,11-trimethyl-","2,2-Dimethyleicosane","Sulfurous acid, 2-ethylhexyl tetradecyl ester")) %>% 
  filter_area(min_maximum = 4e6) %>% 
  filter_ambient_ratio(vol.all, metadata, ratio = 4) %>%
  filter_ambient_ttest(vol.all, metadata, alpha = 0.05, adjust = "fdr") %>%
  combine_filters() 
chemsf$filter_final <- with(chemsf, filter_RT == "OK" & filter_ambient_ratio == "OK" & filters_passed > 5 & filter_contaminant == "OK")# &(freq.kaalae > 0.1 | freq.hookeri > 0.1 |freq.kaho > 0.1))
lapply(chemsf[grepl("filter", names(chemsf))], table)
#write.csv(chemsf, "./data/chemsf_190815_kaho_fixed.csv")

vol <-prune_sampletable(vol.all, chemsf, metadata)

#save.image("sch_workspace190815_fixed.Rdata")

# Finish metadata ---------------------------------------------------------
#Only need to run this section to get final vol matrix and svf metadata

# 2021-12-12 discovered that vol and svf from Schiedea/schiedea_gc_20200205.Rdata and Schiedea/svf_vol.Rdata are not lined up 
# (the latter was used for globosa.Rmd analyses until now)
# because V529 is missing from vol, which normally follows V526 (excluded in code because it fell on the floor).
# Potentially vol <- vol[!(svf$Vial %in% vials_exclude),] was run twice, which deletes the V529 row at the position of V526.
# Scheidea/sch_workspace190815_fixed.Rdata", listed in all versions of scheidea_gc.R, has vol lined up correctly but 
# includes V526, needs to be excluded by running below code.
# for diagnosing this:
# tibble(svf.vial=svf$Vial, svf.fn=svf$FileName, svf.sa=svf$sample, 
# vol.rn=rownames(vol), vol.ac=vol$Acetoin, 
# match=as.integer(svf.fn==vol.rn), match.sa=as.integer(svf.sa==vol.rn)) %>% clipr::write_clip()

load("../sch_workspace190815_fixed.Rdata")

svf <- metadata %>% filter(type == "floral") %>% 
  mutate(sp2 = droplevels(sp2), SampleDate = date, Flrs = amount,
         across(c(Page,Sample,Leaves,Inflor,Mphase,Fphase,Closed,Buds,amount,Flrs,ClosedCN,BudsCN,FlrsCN,Temp), as.integer),
         across(where(is.character), factor),
         across(where(is.difftime), as_hms))

#  add sunset times in the greenhouse
add.sunset <- function(s) {
  library(maptools)
  s$SampleDateStart <- paste(s$SampleDate, s$Start)
  is.na(s$SampleDateStart) <- (is.na(s$SampleDate) + is.na(s$Start))>0
  s$SampleDateStart <- as.POSIXct(s$SampleDateStart)
  s$SampleDateStop <- paste(s$SampleDate, s$Stop)
  is.na(s$SampleDateStop) <- (is.na(s$SampleDate) + is.na(s$Stop))>0
  s$SampleDateStop <- as.POSIXct(s$SampleDateStop)
  
  greenhouse <- matrix(c(-117.8475746, 33.6472914), nrow=1)
  sunsets <- do.call("c", lapply(lapply(s$SampleDateStart, sunriset, crds=greenhouse, direction="sunset", POSIXct.out=T), function(x) x$time))
  s$StartSunset <- difftime(s$SampleDateStart, sunsets,units="hours")
  s$StopSunset <- difftime(s$SampleDateStop,sunsets,units="hours")
  solarnoons <- do.call("c", lapply(lapply(s$SampleDateStart, solarnoon, crds=greenhouse, POSIXct.out=T), function(x) x$time))
  s$StartNoon <- difftime(s$SampleDateStart,solarnoons,units="hours")
  return(s)
}
# Wrong year on this sampledate
svf[as.character(svf$SampleDate)=="2019-08-01","SampleDate"] <- svf[as.character(svf$SampleDate)=="2019-08-01","SampleDate"] - 365

svf <- add.sunset(svf)

# final tweaks to metadata

# Inventory says "Genotype 6 Cutting 2 F" and all other cuttings are labelled female - change sex from M to F to match.
svf$Sex[svf$Population=="13886" & svf$Plant=="6F"] <- "F"

#All the 10228 plants have a three-digit number, so merge 12 and 012
svf$Plant[svf$Population=="10228" & svf$Plant=="12"] <- "012"

# Due to a labelling mistake in the greenhouse, 906 24-2 is actually 905 24-2 
svf$Population[svf$Population=="906" & svf$Plant=="24-2"] <- "905"
#Change label of 905 98s 24-2M to merge its genotype with the above
svf$Plant[svf$Population=="905" & svf$Plant=="98s 24-2M"] <- "24-2"

# Therefore the plants labelled "906 16-3 x 24-2" are actually interisland hybrids, "906 16-3 x 905 24-2"
plants_exclude <- str_detect(svf$Plant, "16-3 x 24-2") %>% replace_na(FALSE)

#Vial 526 fell on floor
vials_exclude <- c(526)

samples_exclude <- (svf$Vial %in% vials_exclude) | plants_exclude #total of 7 samples to exclude
#View(svf[samples_exclude,])
vol <- vol[!samples_exclude,]; svf <- svf[!samples_exclude,] #only run this line once! Don't store it and run it again.

save(svf, vol, file="./data/svf_vol.Rdata")


# S. globosa data only ----------------------------------------------------

library(tidyverse)
library(vegan)
load("./data/svf_vol.Rdata") # lined up, unlike Schiedea/svf_vol.Rdata
svglob <-   svf[svf$sp=="GLOB",]
vol.glob <- vol[svf$sp=="GLOB",]

#table of calibration standards and slopes (area/ng) to use for each type/subtype of compound
standards <- left_join(read_csv("./data/glob_chemtypes.csv"), 
                       read_csv("./data/regressions_2019_prexfer_filtered_slopes.csv"))

#get short names for compounds
chems.glob <- data.frame(name=colnames(vol.glob), mean=colMeans(vol.glob), freq.glob=colMeans(decostand(vol.glob, "pa")), row.names=NULL) #%>% write_csv("./output/glob_chemnames.csv")
chems.glob2 <- read.csv("./data/glob_chemnames - shortnames.csv", stringsAsFactors = F) %>% mutate(name3=ifelse(name2=="",ifelse(IUPAC=="",name,IUPAC), name2), name4=make.unique(name3)) %>% 
  left_join(standards)
chems.glob2$name4[131] <- "linalool oxide (furanoid)"
chems.glob2$name4[184] <- "linalool oxide (furanoid).1"
chems.glob2$name4[36]  <- "linalool oxide (furanoid).2"

colnames(vol.glob)  <- chems.glob2$name4
vol.glob <- sweep(vol.glob, 2, chems.glob2$area_per_ng, FUN = '/')  # convert areas to nanograms 
vol.glob <- vol.glob[, colSums(vol.glob)>0 & chems.glob2$IUPAC != ""]#filter zeros and above 3% of samples (chems with names looked up)
vol.glob <- vol.glob / svglob$Total #divide by equilibration + pumping time
fvol.glob <- vol.glob / sqrt(as.integer(svglob$Flrs))
save(chems.glob2, vol.glob, fvol.glob, svglob, file="./data/vol_glob.rda") #just the S. globosa data for analysis

# Table of volatiles ------------------------------------------------------
add_counts_freqs <- function(chemtable, sampletable, groups) {
  for (g in levels(groups)[levels(groups) != ""]) {
    chemtable[, paste0("count.", g)] <- sapply(na.omit(sampletable[groups == 
                                                                     g, ]), function(x) sum(x > 0))
    chemtable[, paste0("freq.", g)] <- sapply(na.omit(sampletable[groups == 
                                                                    g, ]), function(x) sum(x > 0)/length(x))
  }
  return(chemtable)
}
metadata$sp2_DN <- factor(with(metadata, paste(sp2, DN, sep=".")))
chemsdn <- add_counts_freqs(chemsf, vol.all, metadata$sp2_DN) #%>% 
  #write_csv("./data/chemsdn_190815_fixed.csv")


# Ordinations of filtered table -------------------------------------------

###NMDS of filtered table####
set.seed(1)
nmds.vol <- metaMDS(decostand(vol, "hellinger"), dist="bray", autotransform = FALSE, trymax=5, try=5)
ordiplot(nmds.vol, type = "n", xlim=c(-2,2), ylim=c(-1.5,1.5))
points(nmds.vol, display="sites", col=ifelse(is.na(svf$sp2),"black",spcol[as.integer(svf$sp2)]), pch=c(1,19)[as.integer(svf$DN)]) 
ordiellipse(nmds.vol, svf$sp2, col=spcol, lwd=4)
legend("topleft", levels(svf$sp2), fill=spcol, cex=0.9)

##CAP - species of filtered table##
sch.cap <- capscale(decostand(vol, "hellinger") ~ sp2*DN, data=svf, distance="bray", metaMDSdist = F)
sch.cap
#anova(sch.cap, by="term") #slow
plot(sch.cap, type="n", xlim=c(-1.4, 2), ylim=c(-2.5, 1.5))
points(sch.cap, display="sites", col=spcol[as.integer(svf$sp2)], pch=c(1,19)[as.integer(svf$DN)]) 
ordiellipse(sch.cap, svf$sp2, col=spcol, lwd=2)
#text(sch.cap, display="cn", col=spcol, labels=levels(svf$sp2))# prob only works without interaction
legend("topleft", levels(svf$sp2), fill=spcol, cex=0.8)
legend("topright", levels(svf$DN), pch=c(1,19))
View(sch.cap$CCA$v)

library(vegan3d)
ordirgl(sch.cap, col=spcol[as.integer(svf$sp2)], arr.col=spcol)

## CAP - species at night only
schN.cap <- capscale(sqrt(vol[svf$DN=="Night",]) ~ sp2, data=svf[svf$DN=="Night",], distance="bray", metaMDSdist = F)
schN.cap
cairo_pdf("schiedea_CAP_night.pdf", height=15, width=15)
plot(schN.cap, type="n", xlim=c(-1.4, 1.5), ylim=c(-2, 1.5), main="Night")
points(schN.cap, display="sites", col=spcol[as.integer(svf$sp2[svf$DN=="Night"])], pch=1) 
ordihull(schN.cap, svf$sp2[svf$DN=="Night"], col=spcol, lwd=2)
text(schN.cap, display="cn", col=spcol, labels=levels(svf$sp2), font=2)
dev.off()

## CAP - species at day only
schD.cap <- capscale(sqrt(vol[svf$DN=="Day",]) ~ sp2, data=svf[svf$DN=="Day",], distance="bray", metaMDSdist = F)
schD.cap
cairo_pdf("schiedea_CAP_day.pdf", height=15, width=15)
plot(schD.cap, type="n", xlim=c(-1.4, 2), ylim=c(-2, 2.5), main="Day")
points(schD.cap, display="sites", col=spcol[as.integer(svf$sp2[svf$DN=="Day"])], pch=c(19,15,1)[as.integer(svf$Year)]) 
ordihull(schD.cap, svf$sp2[svf$DN=="Day"], col=spcol, lwd=2)
text(schD.cap, display="cn", col=spcol[-9], labels=levels(svf$sp2)[-9], font=2)#No KAHO daytime samples
legend("topright", levels(svf$Year), pch=c(19,15,1))
dev.off()

## CAP - S. verticillata
vert.cap <- capscale(sqrt(vol[svf$sp2=="VERT",]) ~ DN, data=svf[svf$sp2=="VERT",], distance="bray", metaMDSdist = F)
vert.cap
anova(vert.cap, by="term")
barplot(as.matrix(t(vol[svf$sp2=="VERT",])), las=2, col=rainbow(20))
plot(vert.cap, type="n")
points(vert.cap, display="sites", col=spcol[as.integer(svf$sp2[svf$sp2=="VERT"])], pch=c(1,19)[as.integer(svf$DN[svf$sp2=="VERT"])]) 
text(vert.cap, display="species")
View(vert.cap$CCA$v)


