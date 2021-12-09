library(reshape2)
library(lubridate)
library(hms)
library(tidyr)
library(readr)
library(dplyr)
library(vegan)
library(ggplot2)
library(viridis)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

setwd("~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/Schiedea")
source("../read_shimadzu.R")

# Schiedea runs -----------------------------------------------------------
#Section moved from  ../Inventory/markes_sequence.R, needs to be integrated with rest
Powers <- grepl("Powers RMBL Data", sequ.summary$FullName, fixed=T)
schiedea.ids <- sequ.summary %>% filter(Powers) %>% select(id) %>% unique %>% na.omit
schgc <- sequ.summary %>% filter(id %in% schiedea.ids$id | Powers)
schgc$vial <- as.integer(str_match(schgc$FileName, "_V(\\d{1,3})")[,2])
schgc$dupe <- !is.na(schgc$vial) & schgc$vial %in% schgc$vial[duplicated(schgc$vial)]
#Run  kmeans first
schgc <- cbind(schgc, sch.km[,-1])
rownames(schgc) <- 1:nrow(schgc)
setwd("../Schiedea")
notebook <- read.delim("leak_notebook.csv")
schgc <- cbind(schgc, notebook[match(schgc$vial, as.integer(as.character(notebook$Vial)), incomparables=NA),-4])
schgc$TubeMatch <- schgc$Tube == schgc$SteelTube
write.csv(schgc, "schiedea_all190815.csv")
#read back the verdicts
rmbl <- gs_title("Schiedea RMBL GC-MS Inventory")
schgc.verdict <- gs_read(rmbl)
schgc$verdict <- schgc.verdict$verdict
schgc$FileName2 <- ifelse(is.na(schgc.verdict$sample), schgc$FileName, schgc.verdict$sample)
schgc$vial2 <- as.integer(str_match(schgc$FileName2, "_V(\\d{1,3})")[,2])
schgc$dupe2 <- !is.na(schgc$vial2) & schgc$vial2 %in% schgc$vial2[duplicated(schgc$vial2)]
#schgc.use <- schgc[!(schgc$verdict %in% c("alreadyrun", "empty","leak-blank", "notmine", "notrun", "skip-blank", "skip-notrun")),]
#match up with metadata
sampl <-gs_title("Schiedea Volatiles Sampling")
s <- gs_read(sampl, col_types = cols(.default = col_character(),
                                     Flow = col_double(),
                                     Mass = col_double(),
                                     RunDate = col_date(format = ""),
                                     SampleDate = col_date(format = ""),
                                     Start = col_time(format = ""),
                                     Stop = col_time(format = ""),
                                     Equi = col_double(),
                                     Duration = col_double(),
                                     Total = col_double()
))
s.nr <- s[s$GC=="NR",]
s.nr$Vial <- as.integer(s.nr$Vial)
sort(union(s.nr$Vial, schgc$vial2))
sort(setdiff(s.nr$Vial, schgc$vial2)) #not run yet
sort(setdiff(schgc$vial2, s.nr$Vial)) #run but not entered  
broken <- c(52,53,56,86,117,124,125,129,174) #broken traps - not run
empty <- c(7,16,18,28,40,313,314,315,316) # presumed clean - run
notrun <- c(5,8,21,30,33,50,51,55,57,64,66,74,75,80,81,82,83,84,88,89,91,95,100,102,103,104,107,108,109,11,12,22,23,24,25,26,58,59,60,61,62,63,58,92,94,96,97,98,9, 90, 79, 87, 72, 99)
setdiff(setdiff(setdiff(s.nr$Vial, schgc$vial2), broken), notrun)
s.nr$DN <- factor(ifelse(s.nr$Start > 16*60*60, "Night", "Day"))
#s.nr$FileName <- schgc$FileName2[match(s.nr$Vial, schgc$vial2)]
sv.all <- merge(s.nr, schgc, by.x = "Vial", by.y="vial2", all.x = T, all.y = T)
write.csv(sv.all, "s.nr.schgc190815.csv")
##########################

#sch.data <- read.shimadzu("schiedea190801.txt")
#sch.data2 <- read.shimadzu("Schiedea_190814.txt")
#sch.data <- rbind(sch.data, sch.data2)
#save(sch.data, file="schiedea_190801_190814.Rdata")
load("schiedea_190801_190814.Rdata")

#which names have which CAS numbers?
sch.data %>% filter(CAS.. != "0 - 00 - 0", Name!="") %>% group_by(Name) %>% summarize(CAS=first(CAS..)) %>% 
  mutate(CAS = str_remove_all(CAS, fixed(" "))) %>%  write_csv("schiedea_Name_CAS.csv")

sch.all <- dcast(sch.data, Filename~Name, sum, value.var="Area")
rownames(sch.all) <- sch.all[,1]
sch.all[,1] <- NULL
sch.cut <- sch.all[,colSums(sch.all)>4e6]
rs <- rowSums(sch.cut)
isblank <- grepl("Blank|Bakeout", rownames(sch.cut))#if original FileName says blank

####NMDS and kmeans######
nmds.sch.cut <- metaMDS(decostand(sch.cut, "hellinger"), dist="bray", autotransform = FALSE, trymax=5, try=5)
ordiplot(nmds.sch.cut, type = "n")
points(nmds.sch.cut, display="sites", col=viridis(200)[round(200*sqrt(rs/max(rs)))], pch=as.integer(isblank)+1)

k <- 30
set.seed(1)
km <- kmeans(decostand(sch.cut, "log"), k, nstart=3)
#save(km, file="km30.Rdata")
kblank <- km$cluster %in% c(28,4,18)
sch.km <- data.frame(FileName=row.names(sch.cut), nameBlank= isblank, Mixup= isblank!=kblank, kBlank=kblank, Cluster=km$cluster)[match(schgc$FileName, row.names(sch.cut)),]
with(sch.km, table(kBlank, nameBlank))

par(bg="grey40")
ordiplot(nmds.sch.cut, type = "n",xlim=c(-1.8,2.1),ylim=c(-1.3,1.1))
text(nmds.sch.cut, display="sites", labels=km$cluster, col=kblank+1,cex=0.5 )
points(nmds.sch.cut, display="sites", col=ifelse(kblank, "black", rainbow(k)[km$cluster]), pch=as.integer(isblank)+1) #triangle if original FileName says blank

####CAP - kblanks###
sch.cap.kblank <- capscale(sch.cut ~ as.factor(kblank), distance="bray", metaMDSdist = F)
plot(sch.cap.kblank)
View(sch.cap.kblank$CCA$v)

##### Get metadata ###
sv <- sv.all[match(rownames(sch.cut), sv.all$FileName), ] ##12/3/2020: changed from Fileame2 (renamed) to FileName (original)
sv$badtube <- as.factor(sv$verdict) %in% c("alreadyrun", "empty","leak-blank", "notmine", "notrun", "skip-blank", "skip-notrun","check-chrom","trapblank")
sv$type <- ifelse(grepl("Blank|Bakeout", sv$FileName2),"blank", ifelse(sv$Species=="ambient", "ambient", "floral"))
sv$type[is.na(sv$type)] <- "other"
sv$type <- as.factor(sv$type)
sv$Species[sv$type=="blank"] <- "blank"
sv$Species <- as.factor(sv$Species)
sv$Population <- as.factor(sv$Population)
sv$Plant <- as.factor(sv$Plant)
sv$Cutting <- as.factor(sv$Cutting)
sv$Sex <- as.factor(sv$Sex)
sv$sp <- factor(toupper(substr(as.character(sv$Species), 1, 4)))

ordiplot(nmds.sch.cut, type = "n",xlim=c(-1.8,2.1),ylim=c(-1.3,1.1))
points(nmds.sch.cut, display="sites",col=as.integer(sv$type), pch=1+sv$badtube)
text(nmds.sch.cut, display="sites", labels=as.character(km$cluster), col=as.integer(sv$type)+2*sv$badtube) #shows remaining inconsistencies for nmds above

sch.cut2 <- sch.cut[!sv$badtube & sv$type !="other",]
sv <- sv[!sv$badtube & sv$type !="other",]
sv$type <- droplevels(sv$type)
sv$sp2 <- replace(sv$sp, sv$Population == "KK", "KAAL"); sv$sp2 <- droplevels(replace(sv$sp2, sv$Population == "HH", "HOOK"))
sv$Species2 <- replace(sv$Species, sv$Population == "KK", "kaalae"); sv$Species2 <- replace(sv$Species2, sv$Population == "HH", "hookeri")
sv$Species2Pop <- as.factor(paste(sv$Species2, sv$Population))

with(sv, table(kBlank, type))
with(sv, table(Species2, DN))
with(sv, table(Species2Pop, DN))
with(sv, table(Species2Pop, paste(Sex,DN)))
write.csv(sv, "sv_190815_fixed.csv")
#sv$Mixup2 <- with(sv,(kBlank & type != "blank") | (!kBlank & type=="blank"))
#write.csv(sv[!is.na(sv$Mixup2) & sv$Mixup2==TRUE,], "mixup2.csv")

##NMDS with metadata###
library(ggthemes)
#spcol <- c(gdocs_pal()(10), alpha(gdocs_pal()(10), 0.8), alpha(gdocs_pal()(10), 0.5))
spcol <- unique(c('#3366cc', '#dc3912', '#ff9900', '#109618', '#990099', '#0099c6', '#dd4477', '#66aa00', '#b82e2e', '#316395', '#994499', '#22aa99', '#aaaa11', '#6633cc', '#e67300', '#8b0707', '#651067', '#329262', '#5574a6', '#3b3eac', '#b77322', '#16d620', '#b91383', '#f4359e', '#9c5935', '#a9c413', '#2a778d', '#668d1c', '#bea413', '#0c5922', '#743411', '#3366cc'))

set.seed(1)
nmds.sch.cut2 <- metaMDS(sqrt(sch.cut2), dist="bray", autotransform = FALSE, trymax=5, try=5)
nmds.sch.cut2.plot <- ordiplot(nmds.sch.cut2, type = "n")#, xlim=c(0.6,1), ylim=c(-.2,0.2)
points(nmds.sch.cut2, display="sites", col=ifelse(is.na(sv$Species),"black",spcol[as.integer(sv$Species)]), pch=as.integer(sv$type)) #pch=as.integer(sv$DN)+15
text(nmds.sch.cut2, display="sites", col=ifelse(is.na(sv$Species),"black",spcol[as.integer(sv$Species)]), labels=sv$FileName)
ordiellipse(nmds.sch.cut2, sv$type, col=spcol, lwd=4)
legend("topleft", levels(sv$Species), fill=spcol, cex=0.7)
identify(nmds.sch.cut2.plot, what="sites", labels=rownames(sch.cut2), cex=0.7)

##CAP - species##
sch.cap <- capscale(sqrt(sch.cut2) ~ sv$Species, distance="bray", metaMDSdist = F)
sch.cap.plot <- plot(sch.cap, type="n")
points(sch.cap, display="sites", col=spcol[as.integer(sv$Species)], pch=as.integer(sv$type)) #pch=as.integer(sv$DN)+15
ordiellipse(sch.cap, na.omit(sv$Species), col=spcol, lwd=4)
legend("topleft", levels(sv$Species), fill=spcol, cex=0.7)
identify(sch.cap.plot, what="sites", labels=rownames(sch.cut2[!is.na(sv$Species),]), cex=0.7)
View(sch.cap$CCA$v)

####Filtering######
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
write.csv(chemsf, "chemsf_190815_kaho_fixed.csv")

vol <-prune_sampletable(vol.all, chemsf, metadata)

#save.image("sch_workspace190815_fixed.Rdata")
load("sch_workspace190815_fixed.Rdata")

svf <- metadata[metadata$type == "floral",]
svf$sp2 <- droplevels(svf$sp2)
svf$SampleDate <- svf$date
svf$Flrs <- as.integer(svf$amount)
ints <- c("Page","Sample","Leaves", "Inflor", "Mphase", "Fphase", "Closed", "Buds", "amount", "Flrs", "ClosedCN","BudsCN","FlrsCN","Temp")
svf[,ints] <- lapply(svf[,ints], as.integer)
svf<- mutate_if(svf, is.character, as.factor)
str(svf)

#Vial 526 fell on floor
vials_exclude <- c(526)
vol <- vol[!(svf$Vial %in% vials_exclude),]
svf <- svf[!(svf$Vial %in% vials_exclude),]

###add sunset times to metadata
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
library(hms)
svf <- mutate_if(svf, is.difftime, as.hms)
svf[as.character(svf$SampleDate)=="2019-08-01","SampleDate"]<-svf[as.character(svf$SampleDate)=="2019-08-01","SampleDate"]-365#wrong year
svf <- add.sunset(svf)

save(svf, vol, file="svf_vol.Rdata")

#####Table of volatiles###
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
chemsdn <- add_counts_freqs(chemsf, vol.all, metadata$sp2_DN)
write.csv(chemsdn, "chemsdn_190815_fixed.csv")


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


