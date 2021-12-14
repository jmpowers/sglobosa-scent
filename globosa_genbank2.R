#S globosa sequences and metadata from genbank
library(tidyverse)

#Get an sequence metadata from GenBank flat files
#install.packages("BiocManager")
#BiocManager::install(update=F)
#BiocManager::install(c('BiocGenerics', 'Biostrings', 'GenomeInfoDb', 'GenomicRanges', 'IRanges', 'S4Vectors', 'XVector'), update=F)
#remotes::install_github("gschofl/biofiles")
library(BiocGenerics)
library(biofiles)
library(doParallel)

genes <- c("ncpGS", "pepC", "psbM-trnD")
doParallel::registerDoParallel(cores=4)
dss2df <- function(dss) data.frame(width=BiocGenerics::width(dss), seq=as.character(dss), names=names(dss))
gb <-  gb_seq <- list()
gbr <- lapply(genes, function(g) biofiles::gbRecord(paste0("./data/genbank/",g,"/sequence.gb"))) #slow
names(gbr) <- genes
for (g in genes) {
  gb[[g]] <- as.data.frame(t(as.data.frame(lapply(gbr[[g]], function(x) x@features[[1]]@qualifiers))))
  gb_seq[[g]] <- gbr[[g]] %>% lapply(getSequence) %>% lapply(dss2df) %>% bind_rows()
}

gb.meta <- gb %>% 
  lapply(rownames_to_column,var="accession") %>% 
  bind_rows( .id = "gene") %>% 
  mutate(pop.hap = factor(paste(organism, pop_variant, haplotype))) %>%
  select(-mol_type) %>%
  pivot_wider(names_from=gene, values_from=accession) %>% 
  mutate(pop_variant = recode(pop_variant, "Wood-Huelo"="102HU"),# 102HU == Huelo 
         #Could also lump Kanehoalani: 964 is Kawelo is Wood 13885
         individual = factor(paste(pop_variant, ifelse(
           substr(haplotype, nchar(haplotype),nchar(haplotype)) %in% c("A", "B"), 
           substr(haplotype, 1, nchar(haplotype)-1), 
           haplotype)))) %>% as.data.frame 

#the two nuclear genes (pepC and ncpGS) have up to two alleles that are denoted by A and B after the haplotype
#there are 59 individuals of S globosa, 69 individuals total. Number of seqs: ncpGS:116, pepCL106, cp:69
#to get a list of individuals, just look at the chloroplast gene (psbM-trnD)
gb.meta.cp <- gb.meta %>% 
  filter(!is.na(`psbM-trnD`) & organism=="Schiedea globosa") %>%
  arrange(`psbM-trnD`)
write.csv(gb.meta.cp, "./data/genbank/Wallace2009_Sglobosa_individuals.csv")

######Alignments#####
#read in the downloaded alignment in fasta format (.aln)
#Use POFAD to get genetic distances for diploid organisms 
#Joly, S., Bryant, D. and Lockhart, P.J. 2015. Flexible methods for estimating genetic distances from nucleotide data. Methods in Ecology and Evolution, 6, 938–948.  
#remotes::install_github("simjoly/pofadinr")

library(adegenet)
library(pofadinr)
library(ape)
library(vegan)
aln.paths <- setNames(paste0("./data/genbank/",
                      c("ncpGS/FJ496426.1 and 115 other sequences.aln",
                        "pepC/FJ496542.1 and 105 other sequences.aln", 
                        "psbM-trnD/FJ496357.1 and 68 other sequences.aln")),genes)
ids <- aln.meta <- aln.meta.unique <- con <-aln.con <- genpofaddist <- meand <-tn93dist<- list()
aln <- lapply(aln.paths, fasta2DNAbin) #slow
for (g in genes) {
  ids[[g]] <- substr(rownames(aln[[g]]), 4, 11) #Genbank sequence IDs
  aln.meta[[g]] <- gb.meta[match(ids[[g]],gb.meta[,g]),] 
    # %>% filter(!is.na(organism)) %>% mutate(pop_variant = droplevels(pop_variant)) #NAs are other Schiedea species outgroups
  ##for the nuclear loci, need to represent heterozygosity with IUPAC ambiguity codes
  aln.meta.unique[[g]] <- distinct(aln.meta[[g]], individual, .keep_all=T)
  con[[g]] <- lapply(aln.meta.unique[[g]]$individual, 
                     function(ind) pofadinr::consensus.dna(aln[[g]][aln.meta[[g]]$individual==ind,]))
  names(con[[g]]) <- aln.meta.unique[[g]]$individual
  aln.con[[g]] <- as.alignment(con[[g]])
  genpofaddist[[g]] <- pofadinr::dist.snp(as.DNAbin(aln.con[[g]]), model="GENPOFAD")
  tn93dist[[g]] <-          ape::dist.dna(as.DNAbin(aln.con[[g]]), model="TN93")
  meand[[g]] <- meandist(genpofaddist[[g]], aln.meta.unique[[g]]$pop_variant)
}

        
#Recommended reweighting by Joly and Bruneau 2006 to make each gene have equal weight
tn93dist.max <- lapply(tn93dist, function (x) x/max(x))
tn93dist.avg <- (tn93dist.max[[1]]+tn93dist.max[[2]]+tn93dist.max[[3]])/3
pofaddist.max <- lapply(genpofaddist, function (x) x/max(x))
pofaddist.avg <- (pofaddist.max[[1]]+pofaddist.max[[2]]+pofaddist.max[[3]])/3
#take average of three genes, then take means by population
avgd <- meandist(pofaddist.avg, aln.meta.unique[[3]]$pop_variant)#the levels in gene 3 are not W to E to match meand
#plot(as.vector(avgd), as.vector(meand[[1]]+meand[[2]]+meand[[3]]))

#alternate method for eventual bootstrapping- concatenate consensus sequences, then take distance
con.dnabin <- lapply(aln.con, as.DNAbin)
rownames(con.dnabin[[3]])[31] <- "951 103_F1" #otherwise drops this individual for name mismatch
con.concat <- do.call(cbind.DNAbin,con.dnabin)
con.concat <- con.concat[match(labels(pofaddist.avg), rownames(con.concat)),]#get in same order for comparison
pofad.concat <- pofadinr::dist.snp(con.concat, model="GENPOFAD")

save(meand,genpofaddist,pofaddist.avg,avgd,aln.meta.unique,gb.meta,genes, file="./data/genbank/mean_gen_dist.Rdata")

#Palette for plotting
library(RColorBrewer)
library(scales)
library(viridis)
gen.pops.WtoE <- c("964","844","906","10228","102HU","951","905","850","11221","15191","851","866","793","794","862","896","901")
aln.meta.unique[[1]]$pop_variant <- factor(aln.meta.unique[[1]]$pop_variant, levels=gen.pops.WtoE)
gen.pops.pal <- setNames(c(brewer.pal(5,"Blues")[c(3,4,5)], brewer.pal(4,"Greens")[3:4], brewer.pal(5,"Oranges")[3:5], brewer.pal(4,"Reds")[3:4], brewer.pal(9,"Greys")[3:9]),gen.pops.WtoE)
show_col(gen.pops.pal)

#Compare distances
plot(as.vector(tn93dist[[1]]), as.vector(genpofaddist[[1]]))
plot(as.vector(pofaddist.avg), as.vector(pofad.concat))
plot(as.vector(pofaddist.avg), as.vector(tn93dist.avg))

#Compare two nuclear genes
plot(genpofaddist[[1]],genpofaddist[[3]])
plot(as.vector(meand[[1]]),as.vector(meand[[2]]))

#Heatmaps
heatmap(as.matrix(pofaddist.avg), col=viridis(512))
heatmap(as.matrix(pofad.concat), col=viridis(512))
heatmap(as.matrix((meand[[1]]+meand[[2]]+meand[[3]])/3), col=viridis(512))
heatmap(avgd, col=viridis(512))
plot(avgd)

#Make phylogenies
glob.root <- "896 108_1_7"
phy.tn93 <-  root(nj(tn93dist.avg),  glob.root)
phy.pofad <- root(nj(pofaddist.avg), glob.root)
phy.concat <- root(nj(pofad.concat), glob.root)
#phy.tn93$root.edge <- phy.pofad$root.edge <- phy.concat$root.edge <- 0

#Plot phylo
cairo_pdf("./output/trees/nj_pofad.pdf", width=10, height=10)
par(mar=c(0,0,0,0))
plot(phy.pofad, tip.color=gen.pops.pal[as.integer(aln.meta.unique[[1]]$pop_variant)])
dev.off()

#Bootstrapping concatenated tree
#TODO to bootstrap with the avg method, 
  phy.concat.boot <- boot.phylo(phy.concat, con.concat, function(e) root(nj(pofadinr::dist.snp(e, model="GENPOFAD")),"866 2"), B=100)

cairo_pdf("./output/trees/nj_pofad_concat.pdf", width=10, height=10)#concatenated tree looks very different! 
par(mar=c(0,0,0,0))
plot(phy.concat, tip.color=gen.pops.pal[as.integer(aln.meta.unique[[1]]$pop_variant)])
nodelabels(phy.concat.boot, cex=0.7, bg=viridis(101)[101-phy.concat.boot])
dev.off()

comparePhylo(phy.concat, phy.pofad, plot=T, force.rooted = T)

#Bootstrapping by gene and averaging scaled genetic distances

boot.multilocus.phylo <- function(phy, xlist, DISTFUN, TREEFUN=ape::nj, B = 100, root, rooted = is.rooted(phy)) {#modified from ape::boot.phylo
  boot.tree <- vector("list", B)
  for (i in 1:B) {
    gendist <- list()
    for(g in 1:length(xlist)) {
      boot.samp <- xlist[[g]][, sample.int(ncol(xlist[[g]]), replace = TRUE)]
      gendist[[g]] <- DISTFUN(boot.samp)
    }
    gendist.max <- lapply(gendist, function (x) x/max(x))
    gendist.avg <- (gendist.max[[1]]+gendist.max[[2]]+gendist.max[[3]])/length(gendist.max)
    boot.tree[[i]] <- ape::root(TREEFUN(gendist.avg), root)
  }
  pp <- prop.part(boot.tree)
  return(list(bs=prop.clades(phy, part = pp, rooted = rooted), trees=boot.tree))
}
phy.pofad.boot <- boot.multilocus.phylo(phy.pofad, con.dnabin, DISTFUN=function(x) pofadinr::dist.snp(x, model="GENPOFAD"), root=glob.root, B=100)
phy.pofad$node.label <- phy.pofad.boot$bs

cairo_pdf("./output/trees/nj_pofad_boot.pdf", width=10, height=10)#concatenated tree looks very different! 
par(mar=c(0,0,0,0))
plot(phy.pofad, tip.color=gen.pops.pal[as.integer(aln.meta.unique[[1]]$pop_variant)])
nodelabels("", pch=21, cex=1, frame="none", bg=grey.colors(5, rev=T, start=0, end=1)[cut(phy.pofad$node.label, seq(0,100,by=100/5), include.lowest = T)])
legend("bottomleft",title="Bootstrap value", legend=c("0-20%","21-40%", "41-60%","61-80%","81-100%"), pt.bg=grey.colors(5, rev=T, start=0, end=1), pch=21, pt.cex=2)
dev.off()

#Densitree
library(phangorn)
class(phy.pofad.boot$trees) <- "multiPhylo"
cairo_pdf("./output/trees/nj_pofad_boot_densitree.pdf", width=10, height=10)
phangorn::densiTree(phy.pofad.boot$trees, alpha=0.1, consensus=phy.pofad, scaleX=T)
dev.off()

#NeighborNet
nn.pofad <- phangorn::neighborNet(pofaddist.avg)
cairo_pdf("./output/trees/nn_pofad.pdf", width=20, height=20)
set.seed(11)
plot(nn.pofad, "2D", edge.width=0.5, cex=0.5)
dev.off()


#####Plot phylogeny on map####
library(sf)
library(phytools)
pops <- st_read("./data/Schiedea populations.kml")# http://www.google.com/maps/d/u/0/kml?forcekml=1&mid=1nrFQN_CJTffw2Smd3-dD9C-sbqU&lid=2WoGNCc9HTw") #From Google My Map
pops <- cbind(st_coordinates(pops), pops)
#pops <- pops[-38,] #TODO remove this duplicate
rownames(pops) <- pops$Pop_Number
pops<- pops[pops$Pop_Number != "Nuuanu",]
glob.pops.WtoE <- c("13886","13885","964","844","906","10228","102HU","951","905","850","11221","15191")
glob.pops <- pops[glob.pops.WtoE,]
glob.pops$Pop_Number <- factor(glob.pops$Pop_Number, levels=glob.pops.WtoE)
#glob.pops$Topography <- droplevels(glob.pops$Topography)
glob.pops$Name <- gsub("W&S","Weller & Sakai",glob.pops$Name)
glob.pops$Name[3] <- "Weller & Sakai 964"
all.glob.pops <- pops[pops$Species=="globosa",]
#all.glob.pops$Topography <- droplevels(all.glob.pops$Topography)
#gen.locality <- pops[match(levels(aln.meta.glob$pop_variant), pops$Pop_Number), c("Island", "Pop_Number", "Locality")]
ages <- read.delim("./data/island_Ages.csv")
glob.pops.pal <- setNames(c(brewer.pal(5,"Purples")[c(2,3,3,4,5)], brewer.pal(4,"Blues")[3:4], brewer.pal(5,"Greens")[3:5], brewer.pal(4,"Oranges")[3:4]),glob.pops.WtoE)

just.glob <- aln.meta.unique[[1]]$pop_variant %in% gen.pops.WtoE[1:10]
aln.meta.glob <- aln.meta.unique[[1]][just.glob,]
aln.meta.glob$pop_variant <- droplevels(aln.meta.glob$pop_variant)
aln.meta.glob <- cbind(aln.meta.glob, as.data.frame(pops[match(aln.meta.glob$pop_variant, pops$Pop_Number), ]))
#aln.meta.glob$Topography <- droplevels(aln.meta.glob$Topography)
aln.meta.glob$individual2 <- gsub(" ","_", aln.meta.glob$individual, fixed=T)
rownames(aln.meta.glob) <- aln.meta.glob$individual2
phyg.pofad <- ape::keep.tip(phy.pofad, as.character(aln.meta.glob$individual))
phyg.coord <- as.matrix(aln.meta.glob[match(phyg.pofad$tip.label,aln.meta.glob$individual),c("Y","X")])#makes rownames have underscores as side effect
phyg.pofad <- phytools::minRotate(phyg.pofad, setNames(aln.meta.glob$X, aln.meta.glob$individual)) #this is used by phylo to map when rotating and puts underscores in the tip labels!

spcol <- unique(c('#3366cc', '#dc3912', '#ff9900', '#109618', '#990099', '#0099c6', '#dd4477', '#66aa00', '#b82e2e', '#316395', '#994499', '#22aa99', '#aaaa11', '#6633cc', '#e67300', '#8b0707', '#651067', '#329262', '#5574a6', '#3b3eac', '#b77322', '#16d620', '#b91383', '#f4359e', '#9c5935', '#a9c413', '#2a778d', '#668d1c', '#bea413', '#0c5922', '#743411', '#3366cc'))
cols <- glob.pops.pal#spcol#sample(rainbow(17))#gen.pops.pal

library(mapdata)
phyg.map <- phytools::phylo.to.map(phyg.pofad, phyg.coord, rotate=T, 
                         database="worldHires", region="hawaii", xlim = c(-158.3, -154.8), ylim = c(18.9, 21.72))

source("../plotphylomap2.R") #line 140 to color the labels

cairo_pdf("./output/AppendixS2.pdf", width=9.5, height=12)
plot.phylo.to.map2(phyg.map, 
                   colors=setNames(cols[as.character(aln.meta.glob$pop_variant)], aln.meta.glob$individual2), 
                   lty=1, cex.points=c(0,1), 
                   pch=NA)
#points(mg[-(26:28),])
with(all.glob.pops, points(X,Y, pch=c(21,22,23)[factor(Topography)]))
with(glob.pops, points(X,Y, bg=cols[as.character(Pop_Number)], 
                       pch=c(21,22,23)[Topography]))
nodelabels("", pch=21, cex=1, frame="none", 
           bg=grey.colors(5, rev=T, start=0, end=1)[cut(phy.pofad$node.label, seq(0,100,by=100/5), include.lowest = T)])
legend("topright",title="Bootstrap", 
       legend=c("0-20%","21-40%", "41-60%","61-80%","81-100%"), pt.bg=grey.colors(5, rev=T, start=0, end=1), 
       pch=21, pt.cex=2, bty="n",  title.adj=0)
map.scale(x=-155.65, y=18.9, ratio=F)
legend(x=-158.5, y=20, title="Population",
       legend=with(glob.pops, paste0(Island, "\t", Locality,": ", Name)), 
       pt.bg=cols[as.character(glob.pops$Pop_Number)], 
       pch=c(21,22,23)[factor(glob.pops$Topography)], pt.cex=2, bty="n", title.adj=0)
legend(x=-158.5, y=20.5, title="Topography", 
       legend=sort(unique(glob.pops$Topography)), 
       pch=c(21,22,23), pt.cex=2, pt.bg="white", bty="n", title.adj=0)
with(ages, text(x=X, y=Y, labels=paste(Island, Age, "mya")))
#draw a bracket for Maui Nui
#plotrix::draw.arc(-156.6, 21.1, 1, deg1=310, deg2=190, col="black") #text(x=-156.6, y=20, labels="Maui Nui")
br_lat  <- 20.3
br_long1 <- -157.5
br_long2 <- -155.95
br_height <- 0.08
segments(br_long1,br_lat,br_long2,br_lat)
segments(br_long1,br_lat,br_long1,br_lat+br_height)
segments(br_long2,br_lat,br_long2,br_lat+br_height)
text(x=(br_long1+br_long2)/2, y=br_lat-br_height, labels="Maui Nui")
dev.off()
