options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
library(gt23)    # https://github.com/everettJK/gt23
library(RMySQL)  # loads DBI
library(lubridate)
library(gtools) 
library(wordcloud) 
library(ggplot2) 
library(reldist) 
library(vegan) 
library(argparse)
library(GenomicRanges)
library(dplyr)
library(stringr)

longitudinalCellTypesNum <- 3
minRangeWidth <- 10
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
# if(! 'args' %in% ls()) args <- list()
# args$specimenDB               <- 'specimen_management'
# args$intSiteDB                <- 'intsites_miseq'
# args$reportFile               <- 'report.Rmd'
# args$patient                  <- 'pPatient1'
# args$trial                    <- 'betaThalassemia_sloanKettering_Sadelain'
# args$outputDir                <- 'output'
# args$richPopCells             <- 'PBMC,T CELLS,B CELLS'
# args$cellTypeNameConversions  <- './cellTypeNameConversions.tsv'
# args$numClones                <- 10
# args$use454ReadLevelRelAbunds <- FALSE
# args$legacyData               <- '/home/everett/projects/gtVISA_dashboard/data/legacyData.rds'
# ## args$allowedSamples           <- 'sampleLists/TET2.PMID29849141.samples'


# Convert comma delimited string to vector.
args$richPopCells <- unlist(strsplit(args$richPopCells, ','))


# Retrieve subject data.
dbConn  <- dbConnect(MySQL(), group = args$specimenDB)
q <- paste0('select * from gtsp where patient="', args$patient, '" and trial="', args$trial, '"')
sampleData <- unique(dbGetQuery(dbConn, q))
sampleIDs  <- sampleData$SpecimenAc
dbDisconnect(dbConn)


# Retreive all samples with current intSite data.
dbConn        <- dbConnect(MySQL(), group = args$intSiteDB)
o <- dbGetQuery(dbConn, 'select sampleName, miseqid from samples')

sampleToRunDate <- tibble(GTSP = gsub('\\-\\d+$', '', o$sampleName), 
                           date = lubridate::ymd(stringr::str_extract(o$miseqid, '\\d+'))) %>%
                    group_by(GTSP) %>% 
                    summarise(date = max(date)) %>% 
                    ungroup()


intSiteSamples <- unique(gsub('\\-\\d+$', '', o$sampleName))
dbDisconnect(dbConn)

# Limit samples with intSite data to only those samples associated with this trial / patient.
intSiteSamples <- intSiteSamples[intSiteSamples %in% sampleIDs]

# Read in legacy data.
intSites   <- GRanges()
legacyData <- GRanges()

if(file.exists(args$legacyData)){
  legacyData <- makeGRangesFromDataFrame(subset(readRDS(args$legacyData), 
                                                patient == args$patient & 
                                                trial == args$trial), 
                                         keep.extra.columns = TRUE)
  legacyData$trial <- NULL
}




if(any(sampleIDs %in% intSiteSamples)){
  intSites <- gt23::getDBgenomicFragments(sampleIDs, 'specimen_management', 'intsites_miseq')
  intSites$dataSource <- 'Illumina'
} else if (length(legacyData) > 0){
  intSites <- legacyData
  legacyData <- GRanges()
} else {
  stop('There is no Illumina or legacy data to work with.')
}


# Both Illumina and legacy data is available, merge both into intSites.
if(length(intSites) > 0 & length(legacyData) > 0){
  message(paste0('Merging legacy (', length(legacyData), ' sites) and production (', length(intSites), ' sites) data.'))
    
  # Order the metadata columns of intSites and legacy data so that they can be combined.
  d <- data.frame(mcols(intSites)) 
  mcols(intSites) <- d[, order(colnames(d))]
   
  # Reprocess the time points with the same function used in the gt23 pipeline for consistancy.
  d <- data.frame(mcols(legacyData)) 
  d <- dplyr::select(d, -timePointDays, -timePointMonths, -posid)
  d <- dplyr::bind_cols(d, expandTimePoints(d$timePoint))
  mcols(legacyData) <- d[, order(colnames(d))]
   
  intSites <- unlist(GRangesList(intSites, legacyData))
}


# Remove very short ranges because they are likely not real and may break downstream fragment standardization.
intSites <- intSites[width(intSites) >= minRangeWidth]


if(length(intSites) == 0) stop('No intSite data available.')


# If the allowed samples flag is provided, limit the retreived intsites to those samples. 
if(file.exists(args$allowedSamples)){
  allowedSamples <- gsub('\\s+', '', readLines(args$allowedSamples))
  inSites <- subset(intSites, GTSP %in% allowedSamples)
  if(length(intSites) == 0) stop('Error -- all rows removed after applying allowedSamples filter.')
}


# Hot patch
if(any(grepl('baseline', intSites$timePoint, ignore.case = TRUE))){
  intSites[grepl('baseline', intSites$timePoint, ignore.case = TRUE)]$timePointDays   <- 0
  intSites[grepl('baseline', intSites$timePoint, ignore.case = TRUE)]$timePointMonths <- 0
}

# Check time point conversions.
if(! all(is.numeric(intSites$timePointDays))) 
  stop(paste0('All time points were not convered to numeric values :: ',
              paste0(intSites$timePoint, collapse = ','), ' :: ', 
              paste0(intSites$timePointDays, collapse = ','))) 


# Here we break the typical intSite calling and annoation chain in order to capture
# replicate level so that it can be saved to a data file for later anaylses. 

# Standardize genomic fragments, call intSites, and annotate sites.
intSites.std.reps <- gt23::stdIntSiteFragments(intSites) 


# Multi-hits
#--------------------------------------------------------------------------------------------------
sql <- paste0("select samples.sampleName, samples.refGenome, multihitpositions.multihitID, ",
              "multihitlengths.length from multihitlengths left join multihitpositions on ",
              "multihitpositions.multihitID = multihitlengths.multihitID left join samples on ",
              "samples.sampleID = multihitpositions.sampleID where sampleName like ", 
              paste0(sQuote(paste0(sampleIDs, '-%')), collapse = ' or sampleName like '))

dbConn  <- dbConnect(MySQL(), group = args$intSiteDB)
sites.multi <- unique(dbGetQuery(dbConn, sql))
dbDisconnect(dbConn)

if( nrow(sites.multi) > 0 ){
  
  replicateLevelTotalAbunds <- 
    group_by(data.frame(intSites.std.reps), sampleName, posid) %>%
    mutate(estAbund = n_distinct(width)) %>%
    ungroup() %>%
    group_by(sampleName) %>%
    summarise(totalReplicateAbund = sum(estAbund)) %>%
    ungroup()
  
  sites.multi$GTSP <- str_extract(sites.multi$sampleName, 'GTSP\\d+') 
  
  sites.multi <- left_join(sites.multi, replicateLevelTotalAbunds, by = 'sampleName') %>%
                 group_by(sampleName, multihitID) %>%
                 mutate(estAbund = length(unique(length)),
                        relAbund = ((estAbund / sum(totalReplicateAbund[1], estAbund))*100)) %>%
                 summarise(totalReplicateCells = sum(totalReplicateAbund[1], estAbund), relAbund = relAbund[1], GTSP = GTSP[1]) %>%
                 ungroup() %>%
                 left_join(select(sampleData, CellType, Timepoint, SpecimenAccNum), c("GTSP" = "SpecimenAccNum")) %>%
                 filter(relAbund >= 20) %>%
                 select(-GTSP) %>%
                 arrange(desc(relAbund)) %>%
                 mutate(relAbund = sprintf("%.1f%%", relAbund))

  names(sites.multi) <- c('Replicate', 'Multihit id', 'Total cells in replicate', 'Multihit relative Abundance', 'Celltype', 'Timepoint')
}


intSites <- gt23::collapseReplicatesCalcAbunds(intSites.std.reps) %>%
            gt23::annotateIntSites() %>%
            data.frame()

# Identify which samples were processed but did not return any intSites.
failedIntSiteSamples <- intSiteSamples[! intSiteSamples %in% intSites$GTSP]

# Add annotations to replicate level data.
intSites.std.reps <- data.frame(intSites.std.reps)

possibleFields <- c('posid', 'inFeature', 'nearestFeature', 'nearestFeatureStrand', 'inFeatureExon', 
                    'inFeatureSameOrt', 'nearestFeatureDist', 'nearestOncoFeature', 'nearestOncoFeatureDist', 
                    'nearestOncoFeatureStrand', 'nearestlymphomaFeature', 'nearestlymphomaFeatureDist', 
                    'nearestlymphomaFeatureStrand')

d.reps <- left_join(data.frame(intSites.std.reps), dplyr::select(intSites, possibleFields[possibleFields %in% names(intSites)]), by = 'posid')


# Convert GRange object to data frame and correct cellType names
d <- intSites

# Convert cell types to uppercase and remove leading and trailing white spaces.
d$cellType <- gsub('^\\s+|\\s+$', '', toupper(d$cellType))

# Read in cell type conversion table.
# cellTypeNameConversions <- read.table(args$cellTypeNameConversions, sep='\t', header = TRUE, strip.white = TRUE)
#
# Convert cell type conversion table to uppercase.
#cellTypeNameConversions <- data.frame(apply(cellTypeNameConversions, 2, toupper))

# Duplication of From values will result in duplication of intSite records.
#if(any(duplicated(cellTypeNameConversions$From))) stop('There are duplicate "From" column values in the cellType name conversion file.')

# Join the conversion table to the intSite table, identify which cell types have a conversion and resassign the corresponding values.
#d <- dplyr::left_join(d, cellTypeNameConversions, by = c('cellType' = 'From'))
#d[which(! is.na(d$To)),]$cellType <- d[which(! is.na(d$To)),]$To


# Sync the rep data.frame cell types. 
d.reps$cellType <- NULL
o <- unique(dplyr::select(d, cellType, GTSP))
d.reps <- left_join(d.reps, o, by = 'GTSP')



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
                          # Patient       = patient[1],
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

# Add failed samples.  JKE
if(length(failedIntSiteSamples) > 0){
  f <- tibble(GTSP = failedIntSiteSamples, dataSource = NA, Timepoint = NA, CellType = NA,
              TotalReads = NA, InferredCells = NA, UniqueSites = NA, Gini = NA, Chao1 = NA, Shannon = NA, Pielou = NA,
              UC50 = NA, Included = NA)
  
  f$Timepoint <- sampleData[match(f$GTSP, sampleData$SpecimenAccNum),]$Timepoint
  f$CellType <- sampleData[match(f$GTSP, sampleData$SpecimenAccNum),]$CellType
  #f$timePointDays <- gt23::expandTimePoints(f$Timepoint)$timePointDays
  summaryTable <- dplyr::bind_rows(summaryTable, f)
}

summaryTable$runDate <- sampleToRunDate[match(summaryTable$GTSP, sampleToRunDate$GTSP),]$date

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
d.reps <- d.reps[order(d.reps$timePointDays),]

archivePath <- file.path(args$outputDir, args$patient)
dir.create(archivePath)

# Add patient identifiers to select tables so that they can be merged with other patients if needed.
d.wide$patient         <- args$patient
abundantClones$patient <- args$patient
summaryTable$patient   <- args$patient

write.table(sampleData, sep = ',', col.names = TRUE, row.names = FALSE, quote = TRUE, file = file.path(archivePath, 'sampleData.csv'))
write.table(d, sep = ',', col.names = TRUE, row.names = FALSE, quote = TRUE, file = file.path(archivePath, 'intSites.csv'))
write.table(d.wide, sep = ',', col.names = TRUE, row.names = FALSE, quote = TRUE, file = file.path(archivePath, 'intSites_wideView.csv'))
write.table(d.reps, sep = ',', col.names = TRUE, row.names = FALSE, quote = TRUE, file = file.path(archivePath, 'intSites_replicates.csv'))
write.table(summaryTable, sep = ',', col.names = TRUE, row.names = FALSE, quote = TRUE, file = file.path(archivePath, 'summary.csv'))
write.table(abundantClones, sep = ',', col.names = TRUE, row.names = FALSE, quote = TRUE, file = file.path(archivePath, 'abundantClones.csv'))
write(date(), file = file.path(archivePath, 'timeStamp.txt'))

rmarkdown::render(args$reportFile,
                  output_file = paste0(args$outputDir, '/', args$patient, '/', args$patient, '.pdf'),
                  params = list('date'  = format(Sys.Date(), format="%B %d, %Y"),
                                'title' = paste0('Analysis of integration site distributions and relative clonal abundance for subject ', args$patient)))


unlink(file.path(args$outputDir, args$patient, 'wordClouds'), recursive = TRUE)