extract_plate <- function(path){
  # This function is for importing and cleaning plate output from BioTek
  # Synergy HT plate reader, specifically from the long-format txt export from
  # my <DATE>_hach_sulfide_665.prt runs for the Hach sulfide assay.

  # Dependencies: tidyr
  # For multi-plate files, we need to allow read from string, which takes
  # "text =" argument, so we need to test if path is a vector or scalar
  if(is.atomic(path) && length(path) == 1L){
    plate <- read.table (path, sep = '\t', skip = 51, header = T, fill = T,
                         stringsAsFactors = F, na.strings = '')
  }
  else(plate <- read.table (text = path, sep = '\t', skip = 51, header = T,
                            fill = T, stringsAsFactors = F, na.strings = ''))


  names(plate) <- c('WellID', 'Well', 'StandConc', 'Abs665', 'AbsMinusBlnk',
                    'Concentration', 'NumReps', 'Mean', 'StdDev', 'CV')

  plate <- plate[c('WellID', 'Well', 'StandConc', 'Abs665')]  # Keep 'Concentration' to test our calcs

  # Remove bottom 3 rows, which are standard curve info not matching cols
  plate <- plate[1:(nrow(plate) -3), ]

  # WellID only listed once per set of reps, need to fill in
  # TODO: find a base R cuntion to do this
  plate <- plate %>% tidyr::fill(WellID, .direction = 'down')

  # TO DO (maybe): could make this more generalizable in case we want to
  # include other cols
  plate <- cbind(plate[c('WellID', 'Well')],
                 lapply(plate[c('StandConc', 'Abs665')], numerize_conc,
                        replace = '\\*|<0.000|\\?'))

  plate['SampleType'] <- sub('[0-9]+', '', plate[['WellID']], 1, 3)
  # Another solution to line above in base R. Simple but prob not as robust
  #plate['SampleType'] <- substring(plate[['WellID']], 1, 3)

  plate
}

#extract_plate(paste0(plate_dir, plate_file))
