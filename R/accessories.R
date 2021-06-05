#####
numerize_conc <- function(vector, replace = '<0.000'){
  # This function is for converting microplate data columns to numeric
  # It first converts '<0.000' to NA to avoid coercing values to NA
  # Though this might be safer, as there will be a warning if some other type
  # of value was in the data. I only expect '<0.000' non-numeric values
  # There will also be '?????' and numbers bracketed by *

  if(length(vector[grepl(replace, vector)] > 0)){
    message(paste0(length(vector[grepl(replace, vector)]),
                   ' values replaced by NA: ',
                   paste(unique(vector[grepl(replace, vector)]),
                         collapse = ', ')
    ))
  }

  vector[grepl(replace, vector)] <- NA
  as.numeric(vector)
}

# numerize_conc() test cases
#numerize_conc(c('1', '1', NA, '<0.000'))
#numerize_conc(c(1, 1, NA, '<0.000'))
#numerize_conc(c(1, 1, NA, '<0.000', '*2.490*'), '\\*|<0.000')
#numerize_conc(c(1, 1, NA, '<0.000', '*2.490*', '????'), '\\*|<0.000|\\?')
#numerize_conc(c(1, 1, NA, '<0.000', '*2.490*'), '\\*')
#numerize_conc(c(1, 1, NA, 50))
######

######


get_paths <- function(dir, fnames){
  # This function creates a vector of filepaths, combining the provided
  # directory and the provided filenames
  # NOTE: This function uses regular expression matching, so you can enter
  # incomplete filenames or patterns to return multiple files
  dir(dir,
      pattern = paste(fnames,
                      collapse = '|'),
      full.names = T)
}

######

######

get_curve_dat <- function(plates){
  # This function retrieves data on the standard curve (slope, r2) and complies
  # in a dataframe with the time (minutes)
  # Input is a list of plate objects/lists
  # TODO: make sure this works on a single plate list as well
  #############################################################################
  # Retrieve the standard curve slope and r2 as well as time stamp (as matrix)
  mat <- do.call('rbind',
                 lapply(plates,
                        function(x){c(x[['CalCurve']][c(2,3)],  # slope, r2
                                      x[['RunData']][['FileDetails']][5,2],  # date
                                      x[['RunData']][['FileDetails']][6,2])}))  # timestamp

  # Convert to dataframe for further manipulation and rename cols
  df <- as.data.frame(mat, row.names = F, stringsAsFactors = F)
  names(df) <- c('slope', 'r2', 'date', 'time')  # Rename columns

  # convert columns to appropriate format
  df[['slope']] <- as.numeric(df[['slope']])
  df[['r2']] <- as.numeric(df[['r2']])

  # Calculate elapsed time (minutes) from the hours and minutes fields
  df[['time']] <- strptime(df[['time']],
                           format = "%I:%M:%S %p")$'min' +
    strptime(df[['time']], format = "%I:%M:%S %p")$'hour' * 60

  df[['time']] <- unsplit(
    lapply(split(df, df[['date']]), function(x){
      x[['time']] <- x[['time']] - min(x[['time']], na.rm = T)
    }), df[['date']]
  )
  #df[['time']] <- df[['time']] - min(df[['time']], na.rm = T)
  df
}

######

# TODO: make this function take additional args for the info to append (eg. MDL)
get_processed <- function(plate_ls){
  do.call('rbind', lapply(plate_ls,
                          function(x){tmp <- x[['ProcessedData']]
                          tmp['Date'] <- x[['RunData']][['FileDetails']][5, 2]  # Extract date and append to processed data
                          tmp <- merge(tmp, x[['Metadata']],
                                       by = c('Date', 'WellID'), all = T)
                          tmp}))
}

######
