# alignment_region_dotplot.R 
# takes alignment coordinates as formatted by convert_gnuplot_to_tsv.sh
# produces dotplot of alignment regions 

library(tidyverse)
library(RColorBrewer)

setwd("~/Desktop/prolongata genome paper")

# both infiles generated from mummerplot outputs with convert_gnuplot_to_tsv.sh 
# infile 1 a tsv with R1 Q1 R2 Q2 headers (reg and query start and stop coords)
INFILE_pro = "prolGeneScafs-RAll_QAll-c2500.plot.tsv"
INFILE_rho = "ref_rhop-c1000-Rgenomic_QprolDeDup.plot.tsv"
INFILE_mel = "ref_mel-c500-RMajorScaff_QprolDeDup.plot.tsv"
INFILE_mel_wDups = "melMainScaff_c100_PQ.plot.tsv"
INFILE_DupPairs = "DeDup-Rpairs_Qremoved-c1000.plot.tsv"
# infile 2 a tsv with scaffold breaks in the alignment sequence
INFILE_BREAKS_pro = "prolGeneScafs-RAll_QAll-c2500.breaks.tsv"
INFILE_BREAKS_rho = "ref_rhop-c1000-Rgenomic_QprolDeDup.breaks.tsv"
INFILE_BREAKS_mel = "ref_mel-c500-RMajorScaff_QprolDeDup.breaks.tsv"
INFILE_BREAKS_mel_wDups = "melMainScaff_c100_PQ.breaks.tsv"
INFILE_BREAKS_DupPairs = "DeDup-Rpairs_Qremoved-c1000.breaks.tsv"

# Remove list for prol assembly deduplication
RemovedDups <- read.csv('prolongata_assembly-BUSCO_full_table_RemoveList.tsv', sep = '\t')

# function to get match coordinate and scaffold ref/qry breaks dataframes, in list format
make_l.coords_breaks <- function(INFILE, INFILE_BREAKS){
  # Make df for coordinates of blast hits
  df.coords <- read.csv(INFILE, sep = "\t")
  df.coords <- df.coords[order(df.coords$R1),]
  # get lengths of hits in ref and qry space, and match orientation
  df.coords$ref_length <- df.coords$R2 - df.coords$R1
  df.coords$qry_length <- df.coords$Q2 - df.coords$Q1
  df.coords$orientation <- apply(df.coords, 1, FUN = function(x){
    if( sign(x[5]) == sign(x[6]) ){return("forward")}
    if( sign(x[5]) != sign(x[6]) ){return("reverse")}
  })
  
  # Make df for scaffold positions in reference and query alignment sequences
  df.breaks <- read.csv(INFILE_BREAKS, sep = "\t")
  # using code line from .gp file to find break between Reference and Query scaffolds (gap between lines)
  # find the row with a gap in the code line, will be last line of Reference rows
  breaks.gpLineJumpRow <- which(diff(df.breaks$gpLine)>1)
  df.breaksRef <- df.breaks[1:breaks.gpLineJumpRow,]
  df.breaksQry <- df.breaks[(breaks.gpLineJumpRow+1):nrow(df.breaks),]
  # assign rownames for later value access
  rownames(df.breaksRef) <- df.breaksRef$SeqName
  rownames(df.breaksQry) <- df.breaksQry$SeqName
  # add length for each sequence
  df.breaksRef$SeqLength <- c( diff(df.breaksRef$Position), NaN )
  df.breaksQry$SeqLength <- c( diff(df.breaksQry$Position), NaN )
  
  # Assign Ref and Qry scaffold IDs for each alignment line in df.coords based on position in total alignment
  df.coords$ref_SeqName <- sapply(df.coords$R1, function(x){ 
    df.breaksRef$SeqName[ which( df.breaksRef$Position > x )[1] - 1]
  })
  df.coords$qry_SeqName <- sapply(df.coords$Q1, function(x){ 
    df.breaksQry$SeqName[ which( df.breaksQry$Position > x )[1] - 1]
  })
  
  # drop last row, no Scaffold name (need to keep for previous step assigning ref/qry scaf IDs)
  df.breaksRef <- df.breaksRef[-nrow(df.breaksRef),]
  df.breaksQry <- df.breaksQry[-nrow(df.breaksQry),]
  
  return( list(coords = df.coords, breaksRef = df.breaksRef, breaksQry = df.breaksQry) )
}

# Function to build a df with subset of coordinates from given sequence names in order given, by default returns df.coords
GetScaffoldCoordsDF <- function(L.COORDS, REFERENCE=NULL, QUERY=NULL, 
                                DROP_EMPTY_SCAFFOLDS=T,
                                MIN_REF_LENGTH=0, MIN_QRY_LENGTH=0,
                                MIN_MATCH_LENGTH=0){
  
  COORDS <- L.COORDS$coords
  BREAKS_REF <- L.COORDS$breaksRef
  BREAKS_QRY <- L.COORDS$breaksQry
  
  # Set default REFERENCE and QUERY based on L.COORDS breaks df 
  if( is.null(REFERENCE) ){
    REFERENCE <- BREAKS_REF$SeqName
  }
  if( is.null(QUERY) ){
    QUERY <- BREAKS_QRY$SeqName
  }
  
  # filter qry and ref by length
  REFERENCE <- REFERENCE[BREAKS_REF[REFERENCE,]$SeqLength > MIN_REF_LENGTH]
  QUERY <- QUERY[BREAKS_QRY[QUERY,]$SeqLength > MIN_QRY_LENGTH]

  # get matches for qry and reference
  DF.COORDS_SUBSET <- COORDS[ (COORDS$ref_SeqName %in% REFERENCE) & (COORDS$qry_SeqName %in% QUERY) ,]
  
  # filter based on minimum match length
  DF.COORDS_SUBSET <- DF.COORDS_SUBSET[ (abs(DF.COORDS_SUBSET$ref_length) > MIN_MATCH_LENGTH) &
                                          (abs(DF.COORDS_SUBSET$qry_length) > MIN_MATCH_LENGTH),]
  
  if( DROP_EMPTY_SCAFFOLDS ){
    REFERENCE <- REFERENCE[ REFERENCE %in% DF.COORDS_SUBSET$ref_SeqName ]
    QUERY <- QUERY[ QUERY %in% DF.COORDS_SUBSET$qry_SeqName ]
  }

  # lengths of remaining qry/scaffs
  LENGTHS_REF <- sapply(REFERENCE, function(x) { as.numeric(BREAKS_REF[x,]$SeqLength) }) 
  LENGTHS_QRY <- sapply(QUERY, function(x) { as.numeric(BREAKS_QRY[x,]$SeqLength) })
  
  STARTPOS_REF <- c(1, 1+cumsum(LENGTHS_REF)[-length(LENGTHS_REF)] )
  names(STARTPOS_REF) <- REFERENCE
  STARTPOS_QRY <- c(1, 1+cumsum(LENGTHS_QRY)[-length(LENGTHS_QRY)] )
  names(STARTPOS_QRY) <- QUERY
  

  DF.COORDS_SUBSET$ref_StartPos <- sapply(DF.COORDS_SUBSET$ref_SeqName, function(x) STARTPOS_REF[x])
  DF.COORDS_SUBSET$qry_StartPos <- sapply(DF.COORDS_SUBSET$qry_SeqName, function(x) STARTPOS_QRY[x])
  
  DF.COORDS_SUBSET$ref_SeqLength <- sapply(DF.COORDS_SUBSET$ref_SeqName, function(x) LENGTHS_REF[x])
  DF.COORDS_SUBSET$qry_SeqLength <- sapply(DF.COORDS_SUBSET$qry_SeqName, function(x) LENGTHS_QRY[x])
  
  
  DF.COORDS_SUBSET[,c('R1','R2')] <- DF.COORDS_SUBSET[,c('R1','R2')] - 
    BREAKS_REF[DF.COORDS_SUBSET$ref_SeqName,]$Position + 
    DF.COORDS_SUBSET$ref_StartPos 
  
  DF.COORDS_SUBSET[,c('Q1','Q2')] <- DF.COORDS_SUBSET[,c('Q1','Q2')] - 
    BREAKS_QRY[DF.COORDS_SUBSET$qry_SeqName,]$Position + 
    DF.COORDS_SUBSET$qry_StartPos
  
  return(DF.COORDS_SUBSET)
  
}


PlotDFCoords <- function(L.COORDS, REFERENCE=NULL, QUERY=NULL, 
                         MIN_REF_LENGTH = 0, MIN_QRY_LENGTH = 0, MIN_MATCH_LENGTH = 0, 
                         TIC_LABELS = FALSE, TIC_COUNT = 40, 
                         DROP_EMPTY_SCAFFOLDS = T,
                         ALPHA = 0.25, POINTSIZE = 0.8, COORD_OFFSET = 0.02,
                         LAB_REF = "", LAB_QRY = ""){
  
  COORDS <- L.COORDS$coords
  BREAKS_REF <- L.COORDS$breaksRef
  BREAKS_QRY <- L.COORDS$breaksQry
  
  # Set default REFERENCE and QUERY based on L.COORDS breaks df 
  if( is.null(REFERENCE) ){
    REFERENCE <- BREAKS_REF$SeqName
  }
  if( is.null(QUERY) ){
    QUERY <- BREAKS_QRY$SeqName
  }
  
  # filter qry and ref by length
  REFERENCE <- REFERENCE[BREAKS_REF[REFERENCE,]$SeqLength > MIN_REF_LENGTH]
  QUERY <- QUERY[BREAKS_QRY[QUERY,]$SeqLength > MIN_QRY_LENGTH]
  
  DF_COORDS <- GetScaffoldCoordsDF(L.COORDS = L.COORDS, REFERENCE = REFERENCE, QUERY = QUERY, 
                                   DROP_EMPTY_SCAFFOLDS = DROP_EMPTY_SCAFFOLDS,
                                   MIN_REF_LENGTH = MIN_REF_LENGTH, MIN_QRY_LENGTH = MIN_QRY_LENGTH,
                                   MIN_MATCH_LENGTH = MIN_MATCH_LENGTH)

  
  # By default only plot scaffolds with matches
  if( DROP_EMPTY_SCAFFOLDS ){  
    
    REF_STARTPOS <- unique(DF_COORDS$ref_StartPos)
    REF_SEQNAMES <- unique(DF_COORDS$ref_SeqName)
    QRY_STARTPOS <- unique(DF_COORDS$qry_StartPos)
    QRY_SEQNAMES <- unique(DF_COORDS$qry_SeqName)
    REF_MAX = max(DF_COORDS$ref_StartPos + DF_COORDS$ref_SeqLength)
    QRY_MAX = max(DF_COORDS$qry_StartPos + DF_COORDS$qry_SeqLength)

  } else {
    # otherwise need to replicate some of GetScaffoldCoordsDF() but without filtering for matches
    DF_BREAKSREF <- L.COORDS$breaksRef[REFERENCE,]
    DF_BREAKSQRY <- L.COORDS$breaksQry[QUERY,]
    
    REF_STARTPOS <- c(1, 1+cumsum(DF_BREAKSREF$SeqLength)[-length(REFERENCE)])
    REF_SEQNAMES <- REFERENCE
    QRY_STARTPOS <- c(1, 1+cumsum(DF_BREAKSQRY$SeqLength)[-length(QUERY)])
    QRY_SEQNAMES <- QUERY
    REF_MAX <- DF_BREAKSREF$Position[nrow(DF_BREAKSREF)] + DF_BREAKSREF$SeqLength[nrow(DF_BREAKSREF)]
    QRY_MAX <- DF_BREAKSQRY$Position[nrow(DF_BREAKSQRY)] + DF_BREAKSQRY$SeqLength[nrow(DF_BREAKSQRY)]
    
  }
  
  # find subset of scaffold positions for axis marking, if requested
  if( TIC_COUNT > 0 ){
    REF_TIC_SPACING = as.integer( REF_MAX / TIC_COUNT )
    QRY_TIC_SPACING = as.integer( QRY_MAX / TIC_COUNT )
    
    GetAxisTics <- function(POSITIONS, SPACING) {
      TICGUIDE <- c( seq( min(POSITIONS), max(POSITIONS), by = SPACING ), max(POSITIONS) )
      TICINDEX <- sapply(TICGUIDE, function(guidepoint) which.min( abs(POSITIONS - guidepoint) ) ) %>% unique()
      TICPOSITIONS <- POSITIONS[TICINDEX][unique( c(1, which( diff(POSITIONS[TICINDEX[-length(TICINDEX)]]) > SPACING ), length(TICINDEX)) )] 
      return(TICPOSITIONS[ diff(TICPOSITIONS) > SPACING ])
    }
    
    l.REF_TIC_POSITIONS <- lapply(REF_SEQNAMES, function(seq){ 
      GetAxisTics( sort(c( DF_COORDS[DF_COORDS$ref_SeqName==seq,]$R1, DF_COORDS[DF_COORDS$ref_SeqName==seq,]$R2 )) , REF_TIC_SPACING)
    })
    names(l.REF_TIC_POSITIONS) <- REF_SEQNAMES
    l.REF_TIC_LABELS <- lapply(REF_SEQNAMES, function(seq){ 
      l.REF_TIC_POSITIONS[[seq]] - DF_COORDS[DF_COORDS$ref_SeqName==seq,]$ref_StartPos[1] + 1
    })
    names(l.REF_TIC_LABELS) <- REF_SEQNAMES
    df.REF_TICS <- data.frame(position = unlist(l.REF_TIC_POSITIONS), label = as.character(unlist(l.REF_TIC_LABELS)))
    
    l.QRY_TIC_POSITIONS <- lapply(QRY_SEQNAMES, function(seq){ 
      GetAxisTics(sort(c(DF_COORDS[DF_COORDS$qry_SeqName==seq,]$Q1, DF_COORDS[DF_COORDS$qry_SeqName==seq,]$Q2)), QRY_TIC_SPACING)
    })
    names(l.QRY_TIC_POSITIONS) <- QRY_SEQNAMES
    l.QRY_TIC_LABELS <- lapply(QRY_SEQNAMES, function(seq){ 
      l.QRY_TIC_POSITIONS[[seq]] - DF_COORDS[DF_COORDS$qry_SeqName==seq,]$qry_StartPos[1] + 1
    })
    names(l.QRY_TIC_LABELS) <- QRY_SEQNAMES
    df.QRY_TICS <- data.frame(position = unlist(l.QRY_TIC_POSITIONS), label = as.character(unlist(l.QRY_TIC_LABELS)))
  } else {
    df.REF_TICS <- data.frame(position = c(REF_MAX))
    df.QRY_TICS <- data.frame(position = c(QRY_MAX))
  }
  
  p <- ggplot(DF_COORDS, aes(x=R1, y=Q1))
  p2 <- p + 
    geom_segment(aes(xend=R2, yend=Q2, color=orientation)) + 
    geom_point(aes(color=orientation), alpha=ALPHA, size=POINTSIZE) +
    geom_point(aes(x=R2, y=Q2, color=orientation), alpha=ALPHA, size=POINTSIZE) +
    geom_segment(data = data.frame(x = c(REF_STARTPOS, REF_MAX), xend=c(REF_STARTPOS, REF_MAX), y=0, yend=QRY_MAX), 
                 aes(x=x, xend=xend, y=y, yend=yend), color = "slateblue", linetype='solid', size=0.25) + 
    geom_segment(data = data.frame(x=0, xend=REF_MAX, y = c(QRY_STARTPOS, QRY_MAX), yend = c(QRY_STARTPOS, QRY_MAX)), 
               aes(x=x, xend=xend, y=y, yend=yend), color = "slateblue", linetype='solid', size=0.25) +
    labs(x=LAB_REF, y=LAB_QRY) +
    scale_x_continuous(breaks = c(REF_STARTPOS, df.REF_TICS$position),
                       labels = c(REF_SEQNAMES, rep("", nrow(df.REF_TICS))),
                       minor_breaks = NULL) + #, expand = c(0.05,0), limits = c(-REF_MAX*0.1,REF_MAX*1.1)) +
    scale_y_continuous(breaks = c(QRY_STARTPOS, df.QRY_TICS$position),
                       labels = c(QRY_SEQNAMES, rep("", nrow(df.QRY_TICS))), position = "right" ,
                       minor_breaks = NULL ) + #, expand = c(0.05,0), limits = c(-QRY_MAX*1.1, QRY_MAX)) +
    scale_color_brewer(palette = "Dark2") + 
    theme_minimal() +
    theme(plot.margin = margin(10,0,5,10), axis.title.y.right = element_text(angle = 90),
          axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0, color = "slateblue", face = "bold"), 
          axis.text.y = element_text(angle = 45, hjust = 0, vjust = 0, color = "slateblue", face = "bold") )
  
  if( TIC_LABELS ){
    p2 <- p2 +
      geom_text(data = df.REF_TICS, aes(x = position, y = -1, label = label), angle = -45, hjust = 0, vjust = 1, size = 2) +
      geom_text(data = df.QRY_TICS, aes(x = REF_MAX + 1, y = position, label = label), angle = 45, hjust = 0, vjust = 1, size = 2) +
      coord_cartesian(xlim = c(COORD_OFFSET*REF_MAX, (1+COORD_OFFSET)*REF_MAX), ylim = c(-COORD_OFFSET*QRY_MAX, (1+COORD_OFFSET)*QRY_MAX)) 
  } else {
    p2 <- p2 +
      coord_cartesian(xlim = c(COORD_OFFSET*REF_MAX, REF_MAX*(1-COORD_OFFSET)), ylim = c(COORD_OFFSET*QRY_MAX, (1-COORD_OFFSET)*QRY_MAX)) 
  }
  
  return(p2)
}

# for manuscript fig 3 prol v rhop 
l.coords_rho <- make_l.coords_breaks(INFILE = INFILE_rho, INFILE_BREAKS = INFILE_BREAKS_rho)
# PlotDFCoords(L.COORDS = l.coords_rho, MIN_REF_LENGTH = 2500000, MIN_QRY_LENGTH = 2500000, TIC_COUNT = 0,
#              LAB_REF = "D. rhopaloa scaffolds >2.5Mb", LAB_QRY = "D. prolongata scaffolds >2.5Mb")
PlotDFCoords(L.COORDS = l.coords_rho, MIN_REF_LENGTH = 2500000, MIN_QRY_LENGTH = 1000000, TIC_COUNT = 0,
             LAB_REF = "D. rhopaloa scaffolds >2.5Mb", LAB_QRY = "D. prolongata scaffolds >1Mb")

# for manuscript fig 3 zoom on duplicated regions against rhop
PlotDFCoords(L.COORDS = l.coords_rho, REFERENCE = 'NW_025334990.1', 
             MIN_QRY_LENGTH = 1000000, MIN_MATCH_LENGTH = 7500,
             TIC_COUNT = 40, TIC_LABELS = F,
             LAB_REF = "D. rhopaloa scaffolds", LAB_QRY = "D. prolongata scaffolds")

# for manuscript fig 3 prol v mel 
l.coords_mel <- make_l.coords_breaks(INFILE = INFILE_mel, INFILE_BREAKS = INFILE_BREAKS_mel)
PlotDFCoords(L.COORDS = l.coords_mel, MIN_QRY_LENGTH = 1000000, TIC_COUNT = 0, 
             LAB_REF = "D. melanogaster chromosome arms", LAB_QRY = "D. prolongata scaffolds >1Mb")

# zoom on duplicated regions against mel
PlotDFCoords(L.COORDS = l.coords_mel, REFERENCE = 'NT_033778.4', 
             MIN_QRY_LENGTH = 1000000,
             TIC_COUNT = 40, TIC_LABELS = F, 
             LAB_REF = "D. mel 2R", LAB_QRY = "D. prol scaffolds")

# for manuscript dedup sup fig prol dups and pairs v mel
l.coords_mel_wDups <- make_l.coords_breaks(INFILE = INFILE_mel_wDups, INFILE_BREAKS = INFILE_BREAKS_mel_wDups)
# removed duplicate scaffolds on mel reference
PlotDFCoords(L.COORDS = l.coords_mel_wDups, QUERY = RemovedDups$Remove,
             DROP_EMPTY_SCAFFOLDS = T, TIC_COUNT = 0, COORD_OFFSET = -0.03,
             LAB_REF = "D. mel major scaffolds", LAB_QRY = "Removed Duplicate Scaffolds")

# Duplicate 325 falls on inversion breakpoint w/ mel 2L (but not with rhop)
PlotDFCoords(L.COORDS = l.coords_mel_wDups, QUERY = c('Scaffold_325', 'Scaffold_413'), REFERENCE = 'NT_033779.5',
             DROP_EMPTY_SCAFFOLDS = T, TIC_COUNT = 50, TIC_LABELS = T, COORD_OFFSET = 0.03,
             LAB_REF = "D. mel 2L", LAB_QRY = "")
PlotDFCoords(L.COORDS = l.coords_rho, QUERY = 'Scaffold_413', REFERENCE = 'NW_025335083.1',
             MIN_QRY_LENGTH = 1000000, TIC_COUNT = 50, TIC_LABELS = T, COORD_OFFSET = 0.03,
             LAB_REF = "D. rhop scaffold", LAB_QRY = "D. prol scaffold")


# kept duplicate scaffolds on mel reference (not included in manuscript)
KeptDupPairs <- sub(",.*", "", unique(RemovedDups$Keep)) # have one double pair, but second member is already represented elsewhere
KeptDupPairs.ordered <- KeptDupPairs[order(as.integer(sub("Scaffold_","",KeptDupPairs)))]
PlotDFCoords(L.COORDS = l.coords_mel_wDups, QUERY = KeptDupPairs)

# for manuscript dedup sup fig, removed pairs onto kept pairs
l.coords_DupPairs <- make_l.coords_breaks(INFILE = INFILE_DupPairs, INFILE_BREAKS = INFILE_BREAKS_DupPairs)
KeptDupPairs.ordered <- KeptDupPairs[order(as.integer(sub("Scaffold_","",KeptDupPairs)))]
PlotDFCoords(L.COORDS = l.coords_DupPairs, REFERENCE = KeptDupPairs.ordered, QUERY = RemovedDups$Remove,
             DROP_EMPTY_SCAFFOLDS = T, TIC_COUNT = 0, TIC_LABELS = F,
             LAB_REF = "Retained Scaffolds", LAB_QRY = "Removed Duplicate Scaffolds")


# other example uses
# which prol scaffolds align to mel mitochondria, and where?
PlotDFCoords(L.COORDS = l.coords_mel, REFERENCE = 'NC_024511.2', TIC_LABELS = T,
             LAB_REF = "D. mel mitochonria", LAB_QRY = "D. prol scaffolds", MIN_MATCH_LENGTH = 5000)
# which prol scaffolds align to mel Y chromosome?
PlotDFCoords(L.COORDS = l.coords_mel_wDups, REFERENCE = 'NC_024512.1', TIC_LABELS = F, #MIN_QRY_LENGTH = 1000000,
             LAB_REF = "D. mel Y chromosome", LAB_QRY = "D. prol scaffolds")



# examples for prol v prol alignment
PlotDFCoords(REFERENCE = 'Scaffold_43', QUERY = 'Scaffold_181')
PlotDFCoords(REFERENCE = c('Scaffold_344', 'Scaffold_351'), 
             QUERY = c('Scaffold_116', 'Scaffold_181', 'Scaffold_406', 'Scaffold_412', 'Scaffold_413'),
             MIN_MATCH_LENGTH = 5000)








# Sandbox for defining subsets

LongRefScaffolds = df.breaksRef$SeqName[df.breaksRef$SeqLength > mean(df.breaksRef$SeqLength)]
ScaffoldsOfInterest = c('Scaffold_28', 'Scaffold_69', 'Scaffold_187')
PlotDFCoords(QUERY = ScaffoldsOfInterest, REFERENCE = LongRefScaffolds)
PlotDFCoords(QUERY = LongRefScaffolds, REFERENCE = LongRefScaffolds, MIN_MATCH_LENGTH = 5000)

# This works well for full alignments and subsets of alignment without rearrangements or skipped regions
REFSTART = 12892032.0 # 96312533.0 + 21786421
REFSTOP = 16984303.0 - 1  # 96312533.0 + 21822754
QRYSTART = 61189503.0 # 158561741.0 + 5050478
QRYSTOP = 69406413.0 - 1  # 158561741.0 + 5178153
MIN_MATCH_LENGTH = 5000
ALPHA = 0.25
POINTSIZE = 0.25
LAB_REF = "reference" #"Scaffold_43" # "reference"
LAB_QRY = "query" #"Scaffold_181" # "query"

p <- ggplot(df.coords[df.coords$ref_length > MIN_MATCH_LENGTH,], aes(x=R1, y=Q1))
p + 
  geom_segment(aes(xend=R2, yend=Q2, color=orientation)) + 
  geom_point(aes(color=orientation), alpha=ALPHA, size=POINTSIZE) +
  geom_point(aes(x=R2, y=Q2, color=orientation), alpha=ALPHA, size=POINTSIZE) +
  labs(x=LAB_REF, y=LAB_QRY) + 
#   geom_vline(xintercept = c(REFSTART,REFSTOP), color = "blue") + 
#   geom_hline(yintercept = c(QRYSTART,QRYSTOP), color = "blue") +
  scale_x_continuous(breaks = df.breaksRef$Position, labels = df.breaksRef$SeqName,
                     expand = c(0,0)) +
  scale_y_continuous(breaks = df.breaksQry$Position, labels = df.breaksQry$SeqName,
                     expand = c(0,0), position = "right") +
  scale_color_brewer(palette = "Dark2") + 
  #coord_cartesian(xlim = c(REFSTART, REFSTOP), ylim = c(QRYSTART, QRYSTOP)) + # c(0, max(df.coords$R1))) + # 
  theme_light() +
  theme(plot.margin = margin(10,0,0,40), axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(vjust = 0.5)) 

