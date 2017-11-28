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
              help="Choose overlapping method, options are: loose (match any region with overlaps >=1 BP), stringent (>50% overlap required). [default: %default]"
              , metavar="character"),
  make_option(c("-d", "--database"), type="character", default="VOUS4D.db",
              help="Path to the VOUS4D.db file, default assumes it is in the same folder"
              , metavar="character"),
  make_option(c("-a", "--asd"), type="logical", default=FALSE,
              help="Perform match against ASD data? If TRUE an additional file ([outfile].ASD.txt) will be produced. [default: %default]"
              , metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default="4DG.out.txt",
              help="Output file name '.txt' and '.xls' format are supported. [default: %default]"
              , metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


cat("\n### Loading required packages...\n")
list.of.packages.cran <- c("data.table", "optparse", "plyr", "magrittr", "knitr", "tidyr", "pander", "pkgmaker", "WriteXLS", "tools")
new.packages <- list.of.packages.cran[!(list.of.packages.cran %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cloud.r-project.org/")
suppressMessages(lapply(list.of.packages.cran, require, character.only = TRUE))

# set working directory to the script folder ONLY FOR TESTING
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # ONLY FOR TESTING
#opt$method <- "loose"
#opt$method <- "stringent"
#opt$coordinates<-"chr1:31584329-31656511"
#opt$coordinates <- "dummy.coordinates.txt"

# MATCH BY OVERLAP 
a <- Sys.time()
options(scipen = 99)
if(is.na(opt$coordinates)) {stop("### Please provide at least one input coordinate set, for help: ./VOUSflow --help")}

coords <- fread(paste(opt$coordinates), sep = ":", h=F) %>% 
  .[, c("input.start","input.end") := tstrsplit(V2, "-", fixed = T)] %>%
  .[, V2 := NULL] %>% setnames("V1", "input.chr") %>% 
  .[, lapply(.SD, as.numeric), by=input.chr] %>% 
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
} else if(opt$method == "stringent") {
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
    mg <- rbind(resA,resB, resC, resD); if(nrow(mg) == 0) {next}
    res.temp <- data.frame(coords[i,1:3], mg, row.names = NULL)
    res <- rbind(res, res.temp)
  }
} else {stop("Please, chooose a overlap method between: loose, stringent")}

# remove rows with overlaps in both interactors
res.clean <- res[!(duplicated(res) | duplicated(res, fromLast = TRUE)),]
if(nrow(res.clean) == 0) {stop("### The input region(s) does not contain any interactions in 4DGenome data")}

cat("\n###", nrow(unique(res)),"interactions found in 4Genome\n") 
if(opt$method == "loose") {cat("### Method: LOOSE overlap in H1ESC cells and 4 derived cells\n")
} else {cat("### Method: STRINGENT overlap in H1ESC cells and derived\n")}

cat("### Method: removed", nrow(unique(res))-nrow(res.clean),"matches with both interactors overlapping input coordinates\n")

if(file_extension(paste(opt$outfile) %in% c("xls", "xlsx"))){
  WriteXLS(res.clean, ExcelFileName = paste(opt$outfile), BoldHeaderRow = T, na = "NA")
} else {fwrite(res.clean, file = paste(opt$outfile), sep = "\t")}


if(opt$asd == TRUE) {
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


cat("\n###", length(unique(mapped.genes$genes)),"unique genes mapped\n") 
cat("### Method: removed matches without gene symbol annotation\n")
#kable(mapped.genes)

# match mapped genes with ASD
asd <-merge(mapped.genes, asd, by.x = "genes", by.y = "symbol")
cat("\n###", length(unique(asd$genes)), "input genes found in Genome-wide predictions of autism-associated genes (ASD)","\n")
cat("\n### Top genes (ranked by q-value):")
asd <- asd[order(asd$`q-value`),]
pander(head(asd[,c(1:3,8,9,11)],n = 3), split.cell=15, row.names = F)

# write out file
if(paste(file_extension(paste(opt$outfile)) %in% c("xls", "xlsx"))){
  WriteXLS(asd, ExcelFileName = paste0(file_path_sans_ext(paste(opt$outfile)),".ASD.xls"), BoldHeaderRow = T, na = "NA")
  cat("### ASD output file written to:", paste0(file_path_sans_ext(paste(opt$outfile)),".ASD.xls"), "\n")
  } else {
  fwrite(asd, file = paste0(paste0(file_path_sans_ext(paste(opt$outfile)),".ASD.txt")), sep = "\t")
  cat("### ASD output file written to:", paste0(file_path_sans_ext(paste(opt$outfile)),".ASD.txt"), "\n")
  }

}
cat("### 4D Genome output file written to:", opt$outfile, "\n")
cat("### Analysis completed in", Sys.time() - a, "seconds\n\n")