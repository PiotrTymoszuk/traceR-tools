# This script works on development of a tool fo detecting possible double events
# characterized by a parallel run of their tracks

  library(tidyverse)
  library(furrr)
  library(stringi)
  
# helper functions -----
  
  pair_stats_ <- function(track1, track2) {
    
    ## The heart-piece function calculating stats for two tracks (x1, y1, z1) and (x2, y2, z2) 
    ## Caution: in case the tracks length differ within a cell pair
    ## they're cut down to the shortest one. 

    ## STAT 1: Spearman's correlation between the (x1, x2), (y1, y2), (z1, z2) vectors

    ## STAT 2: the Euclidean distances at the particular time points
    
    ## STAT 3: the mean and SD of the Euclidean distances
    
    
    
    ## checking the input
    
    if(!is.matrix(track1) | !is.matrix(track2)) {
      
      stop('Tracks must be provided as matrices', 
           call. = F)
      
    }
    
    if(!is.numeric(track1) | !is.numeric(track2)) {
      
      stop('Numeric input required', 
           call. = F)
      
    }
    
    if(!identical(colnames(track1), c('t', 'x', 'y', 'z')) | !identical(colnames(track2), c('t', 'x', 'y', 'z'))) {
      
      stop('Corrupted matrix provided. A matrix with columns t, x, y and z required', 
           call. = F)
      
    }
    
    if(nrow(track1) < 3 | nrow(track2) < 3) {
      
      stop('One of the provided tracks is too short. At least three time points required', 
           call. = F)
      
    }
    
    ## adjusting the track length to the shortest of the tracks
    
    rownames(track1) <- track1[, 't']
    rownames(track2) <- track2[, 't']
    
    cmm_rows <- intersect(rownames(track1), 
                          rownames(track2))
    
    if(length(cmm_rows) < 3) {
      
      warning('No sufficient time overalp (at least 3 time points) between the tracks', 
              call. = F)
      
      return(NULL)
      
    }
    
    track1 <- track1[cmm_rows, -1]
    track2 <- track2[cmm_rows, -1]
    
    ## CStatistic 1: creating the x, y and z dimension matrices (ad hoc)
    ## and calculating the correlation coefficients
    
    corr_vec <- map_dbl(c('x', 'y', 'z'), 
                        function(dimension) cor(track1[, dimension], 
                                                track2[, dimension], 
                                                method = 'spearman')) %>% 
      set_names(c('rho_x', 'rho_y', 'rho_z'))
    
    ## Statistic 2 and 3: Euclidean distance calculation
    
    dist_vec <- map_dbl(cmm_rows, 
                        function(timepoint) dist(rbind(track1[timepoint, ], 
                                                       track2[timepoint, ]), 
                                                 method = 'euclidean'))
    
    dist_stats <- c('mean_dist' = mean(dist_vec, na.rm = T), 
                    'sd_dist' = sd(dist_vec, na.rm = T))
                
    return(list(correlation = corr_vec, 
                distances = dist_vec, 
                distance_stats = dist_stats))
    
  }
  
  set_tracks_ <- function(obj_to_name, 
                          obj_names = names(obj_to_name)) {
    
    ## sets custom names: track IDs
    
    name_vec <- stri_split_fixed(obj_names, 
                                 pattern = ':', 
                                 simplify = T)
    
    if(!is.data.frame(obj_to_name)) {
      
      out_tbl <- do.call('rbind', 
                         obj_to_name)
      
      out_tbl <- mutate(as_tibble(out_tbl),
                        track1 = name_vec[, 1], 
                        track2 = name_vec[, 2], 
                        track_pair = obj_names)
      
    } else {
      
      out_tbl <- mutate(obj_to_name,
                        track1 = name_vec[, 1], 
                        track2 = name_vec[, 2], 
                        track_pair = obj_names)
      
      
    }
    
    return(out_tbl) 

  }
  
  sum_abs_ <- function(x, y, z) {
    
    ## calculates a sum of absolute values
    
    abs_sum <- list(x, y, z) %>% 
      map(function(element) ifelse(is.na(element), 
                                   0, 
                                   element)) %>% 
      map(abs) %>% 
      reduce(`+`)
    
    return(abs_sum)
    
  }
  
# class definition and methods ------
  
  new_pair_stat <- function(x = list()) {
    
    ## creates a new S3 'tracer_dbl_check' object
    
    stopifnot(is.list(x))
    stopifnot(identical(names(x), c('correlation', 'distances', 'distance_stats')))
    
    return(structure(x, 
                     class = 'pair_stat'))
    
  }
  
  plot.pair_stat <- function(x, type = NULL, 
                             fill_color = 'cornsilk', 
                             point_shape = 16, 
                             point_color = 'gray40', 
                             scale_x_transf = 'identity', 
                             scale_y_transf = 'identity', ...) {
    
    ## plots the track pair statistics taking a 'pair_stat' object as a compulsory argument
    ## ... are additional arguments passed to the respective geom_*(). 
    
    ## Type:
    ## rho_x - distribution of the Spearman's correlation coefficients for the X positions
    ## rho_y - distribution of the Spearman's correlation coefficients for the Y positions
    ## rho_z - distribution of the Spearman's correlation coefficients for the Z positions
    ## rho_sum - distribution of the sum of absolute values of the Spearman's correlation coefficients for the X, Y and Z positions
    ## rho_density - distribution of the Spearman's correlation coefficients for the X and Y positions
    ## mean_distance - distribution of the mean Euclidean distances between the tracks
    ## sd_distance - distribution of the standard deviations of Euclidean distances between the tracks
    ## distance_stats - SD vs mean Euclidean distance
    
    ## Note: the tracks with low mean distances and low SDs are likely characteristic for the cell doublets
    
    type <- match.arg(type, choices = c('rho_x', 
                                        'rho_y', 
                                        'rho_z', 
                                        'rho_sum', 
                                        'rho_density', 
                                        'mean_distance', 
                                        'sd_distance', 
                                        'distance_stats'))
    
    pair_plot <- switch(type, 
                        
                        ## Distribution of X position correlations
                        
                        'rho_x' = x$correlation %>% 
                          ggplot(aes(x = rho_x)) + 
                          geom_histogram(fill = fill_color, 
                                         color = 'black', ...) + 
                          labs(title = 'Correlation of X positions', 
                               subtitle = expression('Spearman '*rho), 
                               y = '# track pairs', 
                               x = expression(rho*'(X)')), 
                        
                        ## Distribution of Y position correlations
                        
                        'rho_y' = x$correlation %>% 
                          ggplot(aes(x = rho_y)) + 
                          geom_histogram(fill = fill_color, 
                                         color = 'black', ...) + 
                          labs(title = 'Correlation of Y positions', 
                               subtitle = expression('Spearman '*rho), 
                               y = '# track pairs', 
                               x = expression(rho*'(Y)')), 
                        
                        ## Distribution of Z position correlations
                        
                        'rho_z' = x$correlation %>% 
                          ggplot(aes(x = rho_z)) + 
                          geom_histogram(fill = fill_color, 
                                         color = 'black', ...) +  
                          labs(title = 'Correlation of Z positions', 
                               subtitle = expression('Spearman '*rho), 
                               y = '# track pairs', 
                               x = expression(rho*'(Z)')), 
                        
                        ## Distribution of the sum of absolute values of correlation coefficients
                        
                        'rho_sum' = x$correlation %>% 
                          ggplot(aes(x = rho_sum)) + 
                          geom_histogram(fill = fill_color, 
                                         color = 'black', ...) + 
                          labs(title = 'Sum of absolute correlation coefficients', 
                               subtitle = expression('Spearman '*rho), 
                               y = '# track pairs', 
                               x = expression(Sigma*' '*rho)), 
                        
                        ## Density of the X and Y correlation coefficients
                        
                        'rho_density' = x$correlation %>% 
                          ggplot(aes(x = rho_x, 
                                     y = rho_y)) + 
                          geom_density_2d_filled(...) + 
                          labs(title = 'Correlation of X and Y positions', 
                               subtitle = expression('Spearman '*rho), 
                               y = expression(rho*'(X)'), 
                               x = expression(rho*'(Y)')),
                        
                        ## Mean distance between the tracks
                        
                        'mean_distance' = x$distance_stats %>% 
                          ggplot(aes(x = mean_dist)) + 
                          geom_histogram(fill = fill_color, 
                                         color = 'black', ...) + 
                          labs(title = 'Mean distance between the tracks', 
                               subtitle = expression('Euclidean distance'), 
                               y = '# track pairs', 
                               x = 'Mean distance'), 
                        
                        ## SD of the distances between the tracks
                        
                        'sd_distance' = x$distance_stats %>% 
                          ggplot(aes(x = sd_dist)) + 
                          geom_histogram(fill = fill_color, 
                                         color = 'black', ...) + 
                          labs(title = 'Mean distance between the tracks', 
                               subtitle = expression('Euclidean distance'), 
                               y = '# track pairs', 
                               x = 'SD distance'), 
                        
                        ## SD versus mean distance
                        
                        'distance_stats' = x$distance_stats %>% 
                          ggplot(aes(x = mean_dist, 
                                     y = sd_dist)) + 
                          geom_point(shape = point_shape, 
                                     color = point_color, 
                                     fill = fill_color, 
                                     alpha = 0.3) +
                          geom_density_2d(color = 'firebrick', ...) + 
                          labs(title = 'SD and mean distance between the tracks', 
                               subtitle = expression('Euclidean distance'), 
                               y = 'SD distance', 
                               x = 'Mean distance')
                        
                        )
    
    pair_plot <- pair_plot + 
      theme_classic() + 
      theme(panel.grid.major = element_line(color = 'gray90'), 
            plot.tag.position = 'bottom') + 
      scale_x_continuous(trans = scale_x_transf) + 
      scale_y_continuous(trans = scale_y_transf) + 
      labs(tag = paste('\nTrack pairs: n =', 
                       nrow(x$correlation)))
    
    return(pair_plot)
    
  }
  
  summary.pair_stat <- function(object, ...) {
    
    ## obtains summary statistics of correlation coefficients
    ## and distances between the stats for the 'pair_stat' object
    
    return(list(n_pairs = nrow(object$correlation), 
                correlation = summary(object$correlation[c('rho_x', 'rho_y', 'rho_z', 'rho_sum')], ...), 
                distances = summary(object$distance_stats[c('mean_dist', 'sd_dist')], ...)))
    
  }
  
  names.pair_stat <- function(x) {
    
    ## retrieves track pair names from the pair_stat object
    
    return(x$correlation$track_pair)
    
  }
  
  print.pair_stat <- function(x) {
    
    ## prints a summary of a 'pair_stat' object
    
    print_text <- paste('Pair stat object for', 
                        nrow(x$correlation), 
                        'track pairs')
    
    print(print_text)
    
    invisible(nrow(x$correlation))
    
  }
  
# calculation of the statistics between the track pairs -----
  
  get_pair_stats <- function(track_object, .parallel = F) {
    
    ## the function takes a track object (presumably a S3 build at the top of a list)
    ## and calculates peir statistics
    ## Caution: in case the tracks length differ within a cell pair
    ## they're cut down to the shortest one. 
    ## .parallel enables parallel computing via furrr package
    
    ## a S3 object 'pair_stat' with the following components is returned:
    ### correlation: a tibble with Spearman's rho for the x, y and z coordinate for each cell pair
    ## distances: a table with cell-cell distances at particular time points
    ## distance_stats: a table with mean and sd of the distances for each cell pair
    
    if(class(track_object) != 'tracks') {
      
      stop("The function requires a valid object of the 'tracks' class", 
           call. = F)
      
    }
    
    ## user info
    
    start_time <- Sys.time()
    message(paste('Checking for duplicates for', 
                  length(names(track_object)), 
                  ' tracks'))
    on.exit(message(paste('Elapsed:', 
                          Sys.time() - start_time)))
    
    
    ## getting cell pairs (names of the track)
    ## and listing the tracks
    
    cell_pairs <- combn(names(track_object), 
                        m = 2)
    
    pair_lst <- tibble(track1 = cell_pairs[1, ], 
                       track2 = cell_pairs[2, ])
    
    pair_names <- map2_chr(cell_pairs[1, ], 
                           cell_pairs[2, ], 
                           paste, 
                           sep = ':')
    
    track_lst <- unclass(track_object)
    
    ## getting the correlations, distances and distance stats
    
    if(.parallel) {
      
      plan('multisession')
      
      dbl_check <- future_map2(pair_lst$track1, 
                               pair_lst$track2, 
                               function(x, y) pair_stats_(track1 = track_lst[[x]], 
                                                           track2 = track_lst[[y]]))
      
      plan('sequential')
      
    } else {
      
      dbl_check <- map2(pair_lst$track1, 
                        pair_lst$track2, 
                        function(x, y) pair_stats_(track1 = track_lst[[x]], 
                                                    track2 = track_lst[[y]]))
      
    }
    
    dbl_check <- set_names(dbl_check, 
                           pair_names)
    
    dbl_check <- transpose(compact(dbl_check))
    
    ## generating a nice table output and the final 'tracer_dbl_check' object
    
    dbl_check[c('correlation', 
                'distance_stats')] <- dbl_check[c('correlation', 
                                                  'distance_stats')] %>% 
      map(set_tracks_)
    
    dbl_check$distances <- set_tracks_(tibble(distances = dbl_check$distances), 
                                       obj_names = names(dbl_check$distances))
    
    ## obtaining the sum of correlation absolute coefficients for the X, Y and Z positions
    
    dbl_check$correlation <- mutate(dbl_check$correlation, 
                                    rho_sum = sum_abs_(rho_x, rho_y, rho_z))
    return(new_pair_stat(dbl_check))
    
  }
  
  
# testing ------
  
  ### loading the tracer objects provided by Clemens
  
  load('tracks.RData') ## real experiment tracks
  
  ## calculation of the stats
  
  test_stats <- get_pair_stats(track_object = well_1_tracks_filtered, 
                               .parallel = T) ## I cannot guarantee that the parallel computing option works on all platforms...
  
  ## testing the plots
  
  test_plot <- plot(test_stats, 
                    type = 'distance_stats', 
                    scale_x_transf = 'log10', 
                    scale_y_transf = 'log10')
  
  test_plots <- c('rho_x', 
                  'rho_y', 
                  'rho_z', 
                  'rho_sum', 
                  'rho_density', 
                  'mean_distance', 
                  'sd_distance') %>% 
    map(plot, 
        x = test_stats)
  
  ## testing the summary output
  
  summary(test_stats)
  
  ## testing the name retrieval
  
  names(test_stats)
  
  ## printing
  
  test_stats
  
# END -----