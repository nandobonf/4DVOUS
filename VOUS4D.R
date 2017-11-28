#!/usr/bin/env Rscript
cat("\n=======================================================================\n")
cat("||                                                                   ||\n")
cat("||  Script developed by: Ferdinando Bonfiglio (nandobonf@gmail.com)  ||\n")
cat("||                                                                   ||\n")
cat("=======================================================================\n")

options(warn=-1)
list.of.packages.cran <- c("optparse")
new.packages <- list.of.packages.cran[!(list.of.packages.cran %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cloud.r-project.org/")
suppressMessages(lapply(list.of.packages.cran, require, character.only = TRUE))

# create a list of options to provide from the command line
option_list = list(
  make_option(c("-c", "--coordinates"), type="character", default=NA,
              help="path to the text file with input coordinates in the format 'chr:start-stop'"),
  make_option(c("-m", "--method"), type="character", default="loose",
              help="Choose overlapping method, options are: loose (match any region with overlaps >=1 BP), stringent (>50% overlap required), default is loose"
              , metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default="VOUS4D.txt",
              help="Output file name '.txt' and '.xls' format are supported, default is VOUS4D.out.txt"
              , metavar="character"),
  make_option(c("-d", "--database"), type="character", default="VOUS4D.db",
              help="Path to the VOUS4D.db file, default assumes it is in the same folder"
              , metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


cat("\n### Loading required packages...\n")
list.of.packages.cran <- c("data.table", "optparse", "plyr", "magrittr", "knitr", "tidyr", "pander", "pkgmaker", "WriteXLS")
new.packages <- list.of.packages.cran[!(list.of.packages.cran %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cloud.r-project.org/")
suppressMessages(lapply(list.of.packages.cran, require, character.only = TRUE))

# set working directory to the script folder ONLY FOR TESTING
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # ONLY FOR TESTING
#opt$coordinates <-c("chr16:100000-110000, chr22:100000-210000, chr1:15000000-15500000")  # ONLY FOR TESTING
#opt$coordinates <-c("chr16:1000000-1100000") # ONLY FOR TESTING
#opt$coordinates<-"chr11:116730000-116750000"
#opt$method <- "loose"
#opt$method <- "stringent"
#opt$coordinates<-"chr1:31584329-31656511"

# MATCH BY OVERLAP 
a <- Sys.time()
options(scipen = 99)
if(is.na(opt$coordinates)) {stop("### Please provide at least one input coordinate set, for help: ./VOUSflow --help")}

coords <- as.data.table(matrix(unlist(strsplit(opt$coordinates, ":|-| |,|, ")),ncol = 3, byrow = T)) %>% 
  set_names(c("input.chr","input.start","input.end")) %>% .[, lapply(.SD, as.numeric), by=input.chr] %>% 
  setkey(input.chr, input.start, input.end)
cat("### input coordinates:")
kable(coords)

cat("\n### Reading 4DGenome and ASD data ...")
load(paste(opt$database))
cat("DONE\n")
cells <- c("H1ESC", "H1 derived mesenchymal stem cell", "H1 derived mesendoderm cell", "H1 derived neural progenitors", "H1 derived trophoblast")
dat <- dat %>% subset(`Cell/Tissue` %in% cells)

if(opt$method == "loose"){
  res <- rbindlist(list(
    foverlaps(dat, coords, by.x=c('InteractorAChr', 'InteractorAStart', 'InteractorAEnd'), type="any", nomatch = 0L),
    foverlaps(dat, coords, by.x=c('InteractorBChr', 'InteractorBStart', 'InteractorBEnd'), type="any", nomatch = 0L))
    , use.names = T) 
} else if(opt$method == "stringent") {stop("\n### Stringent method not implemented yet!\n")
  coords[,mean:=(input.end+input.start)/2]
  coords <- as.data.frame(coords)
  res = data.table()
  dat <- as.data.frame(dat)
  dat$meanA <- (dat$InteractorAStart+dat$InteractorAEnd)/2
  dat$meanB <- (dat$InteractorBStart+dat$InteractorBEnd)/2
  for(i in 1:nrow(coords)) {
    resA <- dat[dat$InteractorAChr==coords[i,1] & coords[i,4] >= dat$InteractorAStart & coords[i,4] <= dat$InteractorAEnd,]
    resB <- dat[dat$InteractorBChr==coords[i,1] & coords[i,4] >= dat$InteractorBStart & coords[i,4] <= dat$InteractorBEnd,]
    resC <- dat[dat$InteractorAChr==coords[i,1] & dat$meanA >= coords[i,2] & dat$meanA <= coords[i,3],]
    resD <- dat[dat$InteractorBChr==coords[i,1] & dat$meanB >= coords[i,2] & dat$meanB <= coords[i,3],]
    res <- rbind(res, rbind(resA,resB, resC, resD))
  }
} else { 
  stop("Please chooose a overlap method between: loose, stringent")
}


# remove rows with overlaps in both interactors
res.clean <- res[!(duplicated(res) | duplicated(res, fromLast = TRUE)),]
if(nrow(res.clean) == 0) {
  stop("### The input region(s) does not contain any interactions in 4DGenome data")
}

# make a vector with mapped genes, remove NAs from the vector and genes without a HUGO name (ie, ENSG0xxx)
mapped.genes <- data.table(gene=c(res.clean$Agene, res.clean$Bgene)
                           , method = c(res.clean$Detection_Method, res.clean$Detection_Method)
                           , cell = c(res.clean$`Cell/Tissue`, res.clean$`Cell/Tissue`)
) %>% unique() %>% .[!is.na(gene)] %>%
  mutate(genes = strsplit(as.character(gene), ",")) %>% 
  unnest(genes) 
mapped.genes <- mapped.genes[grep("ENSG0", mapped.genes$genes, invert = T),]
mapped.genes <- mapped.genes[grep("ENST0", mapped.genes$genes, invert = T),]
mapped.genes$gene = NULL
mapped.genes <-ddply(mapped.genes, .(genes), summarize, cells=paste(cell,collapse=", "), methods=paste(method,collapse=", "))


cat("\n###", length(unique(mapped.genes$genes)),"unique genes mapped using 4Genome data\n") 
if(opt$method == "loose") {
  cat("### Method: LOOSE overlap in H1ESC cells and derived\n")
} else {
  cat("### Method: STRINGENT overlap in H1ESC cells and derived\n")
}
cat("### Method: removed matches when both interactors overlap input coordinates from", nrow(coords),"input region(s)\n")
cat("### Method: removed matches without gene symbol annotation\n")
#kable(mapped.genes)

# match mapped genes with ASD
asd <-merge(mapped.genes, asd, by.x = "genes", by.y = "symbol")
cat("\n###", length(unique(asd$genes)), "input genes found in Genome-wide predictions of autism-associated genes (ASD)","\n")
cat("\n### Top genes (ranked by q-value):")
asd <- asd[order(asd$`q-value`),]
pander(head(asd[,c(1:3,8,9,11)],n = 3), split.cell=15, row.names = F)

# write out file
if(file_extension(paste(opt$outfile) %in% c("xls", "xlsx"))){
  WriteXLS(asd, ExcelFileName = paste(opt$outfile), BoldHeaderRow = T, na = "NA")
} else {fwrite(asd, file = paste(opt$outfile), sep = "\t")}
cat("### Output written to:", opt$outfile, "\n")
cat("### Analysis completed in", Sys.time() - a, "seconds\n\n")
