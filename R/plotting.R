######

plot_curve <- function(plate, labels = F, wells = F){
  # This function plots the standard curve for a plate
  # Input is a single plate list
  df <- plate[['ProcessedData']]

  exp_df <- df[which(df[['SampleType']] == 'SPL'), ]

  eq_pos <- max(df[['Abs665']], na.rm = T) / 1.5

  p <- ggplot() +
    # Plot regression line (this should be the same as lm results)
    geom_smooth(data = df, aes(x = Abs665, y = StandConc),
                method = 'lm', alpha = 0.2, color = 'black', se = F) +
    # Plot points for standards, including a point for each tech rep
    geom_point(data = df, aes(x = Abs665, y = StandConc,
                              fill = SampleType,
                              shape = SampleType),
               alpha = 0.75, size = 2) +
    # Plot exp samps separate b/c we're using measured conc, not std conc
    geom_point(data = exp_df,
               aes(x = Abs665, y = ConcugL,
                   fill = SampleType, shape = SampleType),
               alpha = 0.75, size = 3) +
    # position = position_jitter(width = 0.02, height = 10
    # Plot standard curve equation
    geom_text(x = eq_pos, y = 100, aes(label = plate[['CurveLab']]), parse = T) +
    scale_shape_manual(values = c(21, 24, 22)) +
    scale_fill_manual(values = c('dodgerblue1','darkgoldenrod1', 'chartreuse3')) +
    labs(fill = 'Sample Type', shape = 'Sample Type') +
    theme_bw()

  if(labels == T){
    p <- p +
      # Plot label for each exp sample
      geom_text(data = exp_df, # Mean Values for samples, standards show indiv replicatess
                aes(x = Abs665, y = ConcugL, label = WellID), hjust = -0.5,
                vjust = 0, size = 3, alpha = 0.5)
  }

  # Show well # of all reps. Meant to help ID where outliers were
  if(wells == T){
    p <- p +
      # Plot label for each exp sample
      geom_text(data = df, # Mean Values for samples, standards show indiv replicatess
                aes(x = Abs665, y = StandConc, label = Well), hjust = -0.5,
                vjust = 0, size = 3, alpha = 0.5)
  }
  p

}


######

make_bar <- function(df, x_var, y_var, runs = NULL, facet_var = NULL,
                     color_var = 'DistillDate'
                     #error_bar
){
  #TODO: for runs arg, allow to specify what col it is matching
  #TODO: get error bars working for this
  if(!is.null(runs)){
    df <- df[which(df[['DistillDate']] %in% runs),]
  }

  p <- ggplot(data = df) +
    geom_col(aes_string(x = x_var, y = y_var, fill = color_var),
             color = 'black') +
    #geom_errorbar(aes_string(x = x_var, ymin = y_var - error_bar,
    #                  ymax = !!y_var + !!error_bar, group = color_var),
    #              width = 0.3) +
    theme_bw() +
    scale_fill_jco()

  if(!is.null(facet_var)){p + facet_grid(.~eval(as.name(facet_var)), scales = 'free')} #special syntax needed to convert string to variable
  else(p)


}

######

make_bar_stack <- function(df, x_var, y_var, runs = NULL, facet_var = NULL,
                           color_var = 'TimeHr'){
  #TODO: for runs arg, allow to specify what col it is matching
  if(!is.null(runs)){
    df <- df[which(df[['DistillDate']] %in% runs),]
  }

  p <- ggplot(data = df) +
    geom_col(aes_string(x = x_var, y = y_var, fill = color_var),
             color = 'black') +
    theme_bw() +
    scale_fill_jco()

  if(!is.null(facet_var)){p + facet_grid(.~eval(as.name(facet_var)), scales = 'free')} #special syntax needed to convert string to variable
  else(p)


}

######
