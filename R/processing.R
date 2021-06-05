######
correct_stds <- function(plate, na2s_mass, vol_ml = 200){
  # Function assumes serial dilutions from this stock
  na2s_mw <- 240.18
  s_mw <- 32.065

  # Ratio of actual to nominal concentration
  t_conc <- na2s_mass*(s_mw / na2s_mw) / (vol_ml/1000) * 10

  plate[['StandConc']] <- plate[['StandConc']] * t_conc

  plate
}

#correct_stds(plate, 0.1528, 200)
######
calc_curve <- function(plate, exclude = c('BLK')){
  stds <- plate[!plate[['WellID']] %in% exclude,]
  mod = lm(Abs665 ~ StandConc, data = plate)
  coeff = c(mod[['coefficients']], summary(mod)[['r.squared']])
  coeff
}

#calc_curve(plate, c('BLK', 'SPL1', 'SPL2'))
######
calc_conc <- function(plate, coeff){
  plate['ConcugL'] <- (plate[['Abs665']] - coeff[1]) / coeff[2]

  plate
}

# Run this on plate dataframe with Concentration column still present
# They should match (and right now (6/4/2020) they do)
#calc_conc(plate, calc_curve(plate, c('BLK', 'SPL1', 'SPL2')))
######
filter_mdl <- function(plate, blanks = c('BLK'), mdl = NULL){
  # Calculate the detection limit based on average of blanks + 3x stDev
  # For now, I'm using this to remove standards (and samples) before
  # calculating the the standard curve, since we have included several
  # standards which are clearly below detection limit, and I would like to
  # avoid skewing the standard curve. Once we have a more formal MDL, we might
  # just enter the official MDL for filtering. We will also probably have
  # standard curve more suitable to the method's range.
  if(is.null(mdl)){
    blnk <- plate[plate[['WellID']] %in% blanks,]

    b_mean <- mean(blnk[['Abs665']])
    b_sd <- sd(blnk[['Abs665']])

    mdl <- b_mean + 3 * b_sd
  }

  # Should this be <= or just < ??
  plate[which(plate[['Abs665']] <= mdl), c('Abs665')] <- NA  # Consider large number if we do post-hoc
  plate
}

#filter_mdl(plate, blanks = c('BLK'))
######
get_run_data <- function(path){
  # This function extracts the run metadata included in the top lines of
  # exported plate data
  dat <- list(FileDetails = NULL, RunDetails = NULL)


  # For multi-plate files, we need to allow read from string, which takes
  # "text =" argument, so we need to test if path is a vector or scalar

  if(is.atomic(path) && length(path) == 1L){
    plate <- read.table(path, sep = '\t', skip = 0, header = F, fill = T,
                        stringsAsFactors = F, na.strings = '')
  }
  else(plate <- read.table(text = path, sep = '\t', skip = 0, header = F, fill = T,
                           stringsAsFactors = F, na.strings = ''))


  dat[['FileDetails']] <- plate[c(1:9),]
  dat[['RunDetails']] <- plate[c(11:15, 136),]

  dat
}

#get_run_data(paste0(plate_dir, plate_file))
######

get_mdl <- function(plate, blanks = c('BLK')){
  # Calculate the detection limit based on average of blanks + 3x stDev
  # Plan to use this as an input for future version of curve calc function
  # That allows a cutoff value for creating the curve. This would allow you
  # to preserve the Abs data (more transparent)
  blnk <- plate[plate[['WellID']] %in% blanks,]

  b_mean <- mean(blnk[['Abs665']])
  b_sd <- sd(blnk[['Abs665']])

  mdl <- b_mean + 3 * b_sd

  mdl
}

#get_mdl(plate, blanks = c('BLK'))
######

process_plate <- function(path, na2s_mass, na2s_vol, metadata = NULL, mdl = NULL){
  # This function is intended to be a wrapper combining several of the plate file
  # processing functions. This is sort of how I'm thinking I could instantiate a
  # plate object in a future version if I make this a package.

  plate_ls <- list(RawData = NULL, RunData = NULL, MDL = NULL,
                   ProcessedData = NULL, Metadata = NULL)

  if(is.character(path)){
    df <- extract_plate(path)
    plate_ls[['RunData']] <- get_run_data(path)
  }
  if(is.list(path)){
    df <- path[['RawData']]
    plate_ls[['RunData']] <- path[['RunData']]
  }

  if(is.null(mdl)){
    mdl = get_mdl(df, blanks = c('BLK'))
  }
  else(mdl)

  # Filter out samples below detection limit
  plate_corrected <- correct_stds(df, na2s_mass, na2s_vol)

  plate_filtered <- filter_mdl(plate_corrected, blanks = c('BLK'), mdl = mdl)

  plate_ls[['CalCurve']] <- calc_curve(plate_filtered,
                                       c('BLK', 'SPL1', 'SPL2'))

  ### Making plottable lm equation from SO link below
  #https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
  #TODO: This could also go inside the standard curve plotting function
  # Is there an advantage to one approach over the other?
  plate_ls[['CurveLab']] <- as.character(
    as.expression(
      substitute(italic(y) == a + b~italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(plate_ls[['CalCurve']][[1]], digits = 3),
                      b = format(plate_ls[['CalCurve']][[2]], digits = 3),
                      r2 = format(plate_ls[['CalCurve']][[3]], digits = 3)))
    ))
  ###

  # Store components of stock standard used for correcting nominal conc
  plate_ls[['StockSoln']] <- c(na2s_mass, na2s_vol)
  names(plate_ls[['StockSoln']]) <- c('na2s_mass_g', 'na2s_vol_ml')

  # Calculate concentrations
  plate_ls[['ProcessedData']] <- calc_conc(plate_filtered,
                                           plate_ls[['CalCurve']])
  plate_ls[['RawData']] <- df

  plate_ls[['MDL']] <- mdl

  if(!is.null(metadata)){
    plate_ls[['Metadata']] <- metadata[which(metadata[['Date']] ==
                                               plate_ls[['RunData']][['FileDetails']][5, 2]),]
  }

  plate_ls
}

#process_plate(paste0(plate_dir, plate_file), 0.1528, 200)
######


split_plates <- function(plt_txt){
  # Produces vector of groupings to apply to splitting text file into separate
  # plates. Takes a text file as input and outputs a vector of length equal to
  # text file, each item denoting the plate that corresponding line of text
  # belongs to. This is meant to help import multi-plate files where the length
  # of plates varies.

  # This is meant to handle plates of varying length within the same file
  # ID the indices of header text indicating start of new plate
  breaks <- which(grepl('Software', plt_txt))

  # Vector of rows occupied by each plate (excluding final plate)
  plt_lengths <- diff(breaks)

  # Get length last plate in th file
  last_length <- length(plt_txt) - sum(plt_lengths)

  # Append length of plate
  plt_lengths <- c(plt_lengths, last_length)




  # Vector defining groupings, used for splitting text file
  grps <- mapply(function(x, y){rep(x, y)}, seq_along(plt_lengths), plt_lengths,
                 SIMPLIFY = F)
  grps <- unlist(grps) # convert matrix to vector
}

#plt_txt <- readLines('microplate_data/20200702_hach_sulfide_665_quick-test.txt')
#print(split_plates(plt_txt = plt_txt))
######
# Had to rethink how to do this, as it produces different levels of nesting
# TODO: Figure out this function so that we can import directly from list of
# File paths with some files having a single and some multiple plates
process_multi_plate <- function(path, na2s_mass, na2s_vol, metadata = NULL, mdl = NULL){
  # This function allows import and processing of palte reader exports with
  # data from multiple plates. It checks for multiple plates in a file and
  # splits accordingly and runs process functions on it. I wonder if it makes
  # more sense to make this a small standalone function that gets wrapped in
  # the process function rather than wrapping the process function in this.


  #fp2 <- 'microplate_data/20200712_hach_sulfide_665.txt'
  #con2 <- file(description = fp2, open = 'rw')

  #txt_lin2 <- readLines(con2)
  #on.exit(close(con2))

  plt_txt <- readLines(path)

  # Get indices of "Software version" header
  # Each plate begins with two blank lines followed by "Software Version..."
  breaks <- which(grepl('Software', plt_txt))

  # Check for multiple plates (maybe we don't actually need this)
  if(length(breaks) > 1){
    # This splits plates according to plate length, but must all be same length
    # as 1st plate
    #plate_length <- diff(breaks)[1]
    #plate_ls <- split(plt_txt, ceiling(seq_along(plt_txt)/plate_length))

    # This allows variable length plates
    plate_ls <- split(plt_txt, split_plates(plt_txt))
  }

  # Need to convert to list so that lapply iterates over each vec instead of
  # vec elements
  else(plate_ls <- list(plt_txt))

  lapply(plate_ls, function(x){process_plate(x, na2s_mass, na2s_vol, metadata, mdl)})

  # could also stream data, reading one line at a time and checking it, then starting new item if it detects software version
}

#fp <- 'microplate_data/20200702_hach_sulfide_665_quick-test.txt'
#process_multi_plate(fp, na2s_mass = 0.075, na2s_vol = 100, metadata = NULL, mdl = NULL)
######

######

reprocess <- function(plate, min_conc, max_conc){
  # This function allows you to reprocess an imported plate list with new
  # parameters. Right now, the only available parameters are the max and min
  # concentrations of standards to use

  raw_data <- plate[['RawData']]

  if(min_conc < plate[['MDL']]){
    mdl <- min_conc
  }
  else(mdl <- plate[['MDL']])

  plate[['RawData']] <- raw_data[which(raw_data[['StandConc']] >= min_conc &
                                         raw_data[['StandConc']] <= max_conc |
                                         is.na(raw_data[['StandConc']])), ]  # preserve exp samples

  new_plate <- process_plate(plate, na2s_mass = plate[['StockSoln']][1],
                             na2s_vol = plate[['StockSoln']][2],
                             metadata = plate[['Metadata']], mdl = mdl)
  new_plate
}

#identical(reprocess(plate_ls[[1]], 50, 1000)[['ProcessedData']], plate_ls[[1]][['ProcessedData']])

######
# TODO: function to spot outliers. Either have message for samples over certain
# SD threshold or message pointing out wells with points outside of some
# threshold



######
combine_multi_plate <- function(path){
  # This function takes a single file path or vector of file paths,
  # imports the text format of each plate file. Data from each plate becomes
  # a list item of just the text. Any files with multiple plates get split into
  # their component plates. This can then be looped over with the process_plate
  # fun. I would like to eventually incorporate this into an import/processing
  # function.
  #TODO: incorporate into an import/processing function


  #fp2 <- 'microplate_data/20200712_hach_sulfide_665.txt'
  #con2 <- file(description = fp2, open = 'rw')

  #txt_lin2 <- readLines(con2)
  #on.exit(close(con2))

  plt_txt <- unlist(lapply(path, function(x){readLines(x)}))

  # Get indices of "Software version" header
  # Each plate begins with two blank lines followed by "Software Version..."
  breaks <- which(grepl('Software', plt_txt))

  pls <- list(length = length(breaks))

  # Check for multiple plates (maybe we don't actually need this)
  if(length(breaks) > 1){
    # This allows variable length plates
    plate_ls <- split(plt_txt, split_plates(plt_txt))
  } else(plate_ls <- list(plt_txt))
  # Need to convert to list so that lapply iterates over each vec instead of
  # vec elements
  plate_ls

}
