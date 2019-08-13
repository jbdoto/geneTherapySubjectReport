options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
library(gt23)    # https://github.com/everettJK/gt23
library(RMySQL)  # loads DBI
library(gtools) 
library(wordcloud) 
library(ggplot2) 
library(reldist) 
library(vegan) 
library(argparse)
library(GenomicRanges)
library(dplyr)

longitudinalCellTypesNum <- 3
maxWordCloudWords <- 100
ppNum <- function(n){ if(is.na(n)) return(''); format(n,big.mark=",", scientific=FALSE, trim=TRUE) }

parser <- ArgumentParser()
parser$add_argument("--CPUs", type="integer", default=20, help="Number of CPUs to use while compiling data for the report")
parser$add_argument("--trial", type = "character", help = "Trial identifier")
parser$add_argument("--patient", type = "character", help = "Patient identifier")
parser$add_argument("--outputDir", type = "character", default = "output", help = "Output directory")
parser$add_argument("--numClones", type="integer", default=10, help="Number of most abundant clones to show in relative abundance plots")
parser$add_argument("--intSiteDB",  type="character", default="intsites_miseq", help="IntSite database group name")
parser$add_argument("--reportFile", type="character", default="report.Rmd", help="Rmd file used to build reports")
parser$add_argument("--specimenDB", type="character", default="specimen_management", help="Specimen database group name")
parser$add_argument("--legacyData", type="character", default='legacyData.rds', help="Legacy fragment level data with read counts, ie. 454 sequencing, which can be subseted by patient and trial.")
parser$add_argument("--richPopCells", type='character', default='Whole blood,T cells,B cells,NK cells,Neutrophils,Monocytes,PBMC', help='Cell types to report rich populations')
parser$add_argument("--allowedSamples", type="character", default='allowedSamples.list', help="List of GTSP ids which are allowed to be considered in the report.")
parser$add_argument("--cellTypeNameConversions", type = "character", default = "cellTypeNameConversions.tsv", help = "File containing cell type name conversions")
parser$add_argument("--use454ReadLevelRelAbunds", action='store_true',  help="Use the number of reads per intSite rather than unique fragments when calculating relative abundances.")
args <- parser$parse_args()


# IDE overrides for when working within an IDE.
if(! 'args' %in% ls()) args <- list()
args$specimenDB               <- 'specimen_management'
args$intSiteDB                <- 'intsites_miseq'
args$reportFile               <- 'report.Rmd'
args$patient                  <- 'pP9'
args$trial                    <- 'SCID1_Paris_Cavazzana'
args$outputDir                <- 'output'
args$richPopCells             <- 'Whole blood,T cells,B cells,NK cells,Neutrophils,Monocytes,PBMC'
args$cellTypeNameConversions  <- './cellTypeNameConversions.tsv'
### args$allowedSamples           <- 'sampleLists/TET2.PMID29849141.samples'
args$legacyData               <- '/home/everett/projects/project.management.dashboard/data/legacyData.rds'
args$use454ReadLevelRelAbunds <- TRUE


# Convert comma delimited string to vector.
args$richPopCells <- unlist(strsplit(args$richPopCells, ','))


# Read in legacy data.
legacyData <- GRanges()
if(file.exists(args$legacyData)){
  legacyData <- makeGRangesFromDataFrame(subset(readRDS(args$legacyData), 
                                                patient == args$patient & 
                                                  trial == args$trial), 
                                        keep.extra.columns = TRUE)
  legacyData$trial <- NULL
}


# Connect to sample database are retrieve subject data.
dbConn  <- dbConnect(MySQL(), group = args$specimenDB)
q <- paste0('select * from gtsp where patient="', args$patient, '" and trial="', args$trial, '"')
sampleData <- unique(dbGetQuery(dbConn, q))
sampleIDs  <- sampleData$SpecimenAccNum
dbDisconnect(dbConn)


# Merge legacy data if available.
if(nrow(sampleData) == 0){
   intSites <- legacyData
} else {
  intSites <- gt23::getDBgenomicFragments(sampleIDs, 'specimen_management', 'intsites_miseq')
  intSites$dataSource <- 'Illumina'
  
  if(length(intSites) > 0 & length(legacyData) > 0){
     message(paste0('Merging legacy (', length(legacyData), ' sites) and production (', length(intSites), ') data.'))
    
     # Order the metadata columns of intSites and legacy data so that they can be combined.
     d <- data.frame(mcols(intSites)) 
     mcols(intSites) <- d[, order(colnames(d))]
   
     # Reprocess the time points with the same function used in the gt23 pipeline for consistancy.
     d <- data.frame(mcols(legacyData)) 
     d <- dplyr::select(d, -timePointDays, -timePointMonths)
     d <- dplyr::bind_cols(d, expandTimePoints(d$timePoint))
     mcols(legacyData) <- d[, order(colnames(d))]
   
     intSites <- unlist(GRangesList(intSites, legacyData))
   }
}


if(length(intSites) == 0) stop('No intSite data available.')


# If the allowed samples flag is provided, limit the retreived intsites to those samples. 
if(file.exists(args$allowedSamples)){
  allowedSamples <- gsub('\\s+', '', readLines(args$allowedSamples))
  inSites <- subset(intSites, GTSP %in% allowedSamples)
  if(length(intSites) == 0) stop('Error -- all rows removed after applying allowedSamples filter.')
}


# Check time point conversions.
if(! all(is.numeric(intSites$timePointDays))) 
  stop(paste0('All time points were not convered to numeric values :: ',
              paste0(intSites$timePoint, collapse = ','), ' :: ', 
              paste0(intSites$timePointDays, collapse = ','))) 


# Standardize genomic fragments, call intSites, annd annotate sites.
intSites <- gt23::stdIntSiteFragments(intSites) %>%
            gt23::collapseReplicatesCalcAbunds() %>%
            gt23::annotateIntSites()


# Convert GRange object to data frame and correct cellType names
d <- data.frame(intSites)

# Convert cell types to uppercase and remove leading and trailing white spaces.
d$cellType <- gsub('^\\s+|\\s+$', '', toupper(d$cellType))

# Read in cell type conversion table.
cellTypeNameConversions <- read.table(args$cellTypeNameConversions, sep='\t', header = TRUE, strip.white = TRUE)

# Convert cell type conversion table to uppercase.
cellTypeNameConversions <- data.frame(apply(cellTypeNameConversions, 2, toupper))

# Join the conversion table to the intSite table, identify which cell types have a conversion and resassign the corresponding values.
d <- dplyr::left_join(d, cellTypeNameConversions, by = c('cellType' = 'From'))
d[which(! is.na(d$To)),]$cellType <- d[which(! is.na(d$To)),]$To



# Build relative abundace plot data and plots
#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~

d <- dplyr::group_by(d, timePoint, cellType) %>%
     dplyr::mutate(readsPerSample = sum(reads)) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(timePoint, cellType, posid) %>%
     dplyr::mutate(readsRelAbund = (sum(reads) / readsPerSample[1])*100) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(cellType, timePoint, GTSP) %>%
     dplyr::mutate(totalSampleFrags = sum(estAbund)) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(cellType, timePoint) %>%
     dplyr::mutate(include = ifelse(totalSampleFrags == max(totalSampleFrags), 'yes', 'no')) %>%
     dplyr::ungroup()

# Swap in readsRelAbund values for relAbund vlues for 454 data.
if(args$use454ReadLevelRelAbunds) d <- dplyr::mutate(d, relAbund = ifelse(dataSource == '454', readsRelAbund, relAbund))

# Add nearest feature flags.
d <- d %>%
  mutate(labeledNearestFeature = paste0(nearestFeature, ' ')) %>% 
  mutate(labeledNearestFeature = ifelse(inFeature, paste0(labeledNearestFeature, '*'), labeledNearestFeature)) 

if('nearestOncoFeatureDist' %in% names(d))
  d <- mutate(d, labeledNearestFeature = ifelse(abs(nearestOncoFeatureDist) <= 50000, paste0(labeledNearestFeature, '~'), labeledNearestFeature))

if('nearestlymphomaFeatureDist' %in% names(d))
  d <- mutate(d, labeledNearestFeature = ifelse(abs(nearestlymphomaFeatureDist) <= 50000, paste0(labeledNearestFeature, '!'), labeledNearestFeature)) 


# Create label for unique intSites.
d$posidLabel <- paste0(d$labeledNearestFeature, '\n', d$posid)


# Select cell types of interest which will be used to focus the selection of which 
# clones to show in clones expanding over time overview.
args$cellTypes <-
  dplyr::group_by(d, GTSP) %>%
  dplyr::mutate(maxAbund = sum(estAbund)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(maxAbund >= 50) %>%
  dplyr::group_by(cellType) %>%
  dplyr::summarise(timePoints = n_distinct(timePoint),
                   maxRelAbund = floor(max(relAbund))) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(maxRelAbund)) %>%
  dplyr::pull(cellType)

args$cellTypes <- args$cellTypes[1:ifelse(length(args$cellTypes) < longitudinalCellTypesNum, 
                                          length(args$cellTypes), 
                                          longitudinalCellTypesNum)]


# Build summary table.
#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~

calculateUC50 <- function(abund){
  stopifnot(is.vector(abund) & is.numeric(abund))
  abund <- abund[order(abund)]
  accum <- sapply(1:length(abund), function(i){sum(abund[1:i])})
  length(accum[accum >= sum(abund)/2])
}

summaryTable <- group_by(d, GTSP) %>%
                summarise(dataSource    = dataSource[1],
                          timePointDays = timePointDays[1],
                          # Replicates    = n_distinct(sampleName),
                          Patient       = patient[1],
                          Timepoint     = timePoint[1],
                          CellType      = cellType[1],
                          TotalReads    = ppNum(sum(reads)),
                          InferredCells = ppNum(sum(estAbund)),
                          UniqueSites   = ppNum(n_distinct(posid)),
                          Gini          = sprintf("%.3f", gini(estAbund)),
                          Chao1         = ppNum(round(estimateR(estAbund, index='chao')[2], 0)),
                          Shannon       = sprintf("%.2f", diversity(estAbund)),
                          Pielou        = sprintf("%.3f", diversity(estAbund)/log(specnumber(estAbund))),
                          UC50          = ppNum(calculateUC50(estAbund)),
                          Included      = include[1]) %>%
               ungroup() %>%
               arrange(timePointDays) %>%
               select(-timePointDays)


# Now that we have report all the samples, remove those which will not be used for plotting.
d <- subset(d, include == 'yes')


# Create data frame needed to generate relative abundance plots.
abundantClones <- bind_rows(lapply(split(d, d$cellType), function(x){
  
  # Adjust the number of clones to return based on the number of sites per cell type.
  if(nrow(x) < args$numClones) args$numClones <- nrow(x)
  
  # Sort nearest genes by abundance.
  x <- x[order(x$estAbund, decreasing = TRUE),]
  
  # Select clones to report.
  topClones <-  unique(x$posidLabel)[1:args$numClones]
  
  # For each time point, create a data frame for relative abundance plots
  bind_rows(lapply(split(x, x$timePoint), function(x2){
    
    lowAbundData <- dplyr::mutate(x2, posidLabel = 'LowAbund',
                                      totalCells = sum(estAbund),
                                      relAbund   = 100) %>%
                    dplyr::slice(1) %>% 
                    dplyr::select(cellType, timePoint, dataSource, posidLabel, totalCells, timePointDays, relAbund)

    x3 <- subset(x2, posidLabel %in% topClones)
    if(nrow(x3) == 0) return(lowAbundData)
    x3$totalCells <- sum(x2$estAbund)
    
    lowAbundData$relAbund <- 100 - sum(x3$relAbund)
    bind_rows(lowAbundData,  dplyr::select(x3, cellType, timePoint, dataSource, posidLabel, totalCells, timePointDays, relAbund))
  }))
}))



# Create word cloud figures.
#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~

unlink(file.path(args$outputDir, args$patient, 'wordClouds'), recursive = TRUE)
dir.create(file.path(args$outputDir, args$patient, 'wordClouds'), recursive = TRUE)

invisible(lapply(split(d, paste(d$cellType, d$timePoint)), function(x){
  #browser()
  x <- x[order(x$estAbund, decreasing = TRUE),]

  # Handle older data which does not have estAbund values.  
  if(all(is.na(x$estAbund))) x <- x[order(x$reads, decreasing = TRUE),]

  if(nrow(x) < maxWordCloudWords) maxWordCloudWords <- nrow(x)
  x <- x[1:maxWordCloudWords,]
  
  # Handle older data which does not have estAbund values.
  if(all(is.na(x$estAbund))){
    w <- setNames(x$reads, x$labeledNearestFeature)
  } else {
    w <- setNames(x$estAbund, x$labeledNearestFeature)
  }
  
  w <- w[! is.na(w)]
  
  cellType <-  gsub('\\_|\\/', ' ', x$cellType[1])
  
  png(file = file.path(args$outputDir, args$patient, 'wordClouds', paste0(x$patient[1], '_', cellType, '_', x$timePoint[1], '_', 
                    x[maxWordCloudWords,]$estAbund, '_', x[1,]$estAbund, '_wordCloud.png')), res = 150)
  
  wordcloud(names(w), w, random.order=FALSE, colors=colorRampPalette(brewer.pal(12, "Paired"))(maxWordCloudWords), rot.per=0, max.words = 50, scale = c(2, 0.2))
  dev.off()
}))


# Additional data objects that others may find useful

d.wide <- reshape2::dcast(dplyr::select(d, posid, cellType, timePoint, estAbund), posid ~ cellType+timePoint, value.var='estAbund', fill=0)  
names(d.wide) <- stringi::stri_replace_last(str = names(d.wide), regex = "_", replacement = "/")

# Save data objects and render the report.
#~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~o~-~-~-~-~-~-~-~-~-~-~-~-~

d <- d[order(d$timePointDays),]
save(args, sampleData, d, d.wide, summaryTable, abundantClones, file=file.path(args$outputDir, paste0(args$patient, '.RData')))

rmarkdown::render(args$reportFile, 
                  output_file = paste0(args$outputDir, '/', paste0(args$trial, '.', args$patient), '.pdf'),
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = paste0('Analysis of integration site distributions and relative clonal abundance for subject ', args$patient)))
