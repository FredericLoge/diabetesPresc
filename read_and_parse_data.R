#' @title Read and parse the UCI diabetes data
#' This function first checks that the data is indeed stored where it is expected.
#' If so it loads it, parse it and does some refactoring.
#' If it does not find it, jsut returns an error.
read_and_parse_diabetes_data <- function(){
  
  # checking out repository
  repos <- './datasets/diabetes_uci_data/'
  repos_exists <- file.exists(repos)
  
  ##
  if( repos_exists == FALSE ){
    
    stop( paste0('Repository "', repos, '" does not exist. Data could not be loaded.') )
    
  }else{
    
    # read data from repository
    df_list <- list()
    for(i in 1:70){
      j <- as.character(i)
      if(i<10) j <- paste0('0',i)
      temp <- tryCatch({
        read.table(paste0(repos, 'data-',j) , stringsAsFactors = FALSE, na.strings = c('0Hi','0Lo', 'NA', "0''")) 
      }, error = function(e) return(NULL))
      if(!is.null(temp)){
        df_list[[length(df_list)+1]] <- cbind('id' = i, temp)
      }
    }
    df <- do.call(rbind.data.frame, df_list)
    colnames(df) <- c('id', 'date', 'time', 'code', 'value')
    
    # parse date and time into proper date format
    df$timestamp <- as.POSIXct(x = paste0(df$date, ' ', df$time), format = '%m-%d-%Y %H:%M')
    
    # mapping code and labels
    code_labels <- c(
      '4' = '1?',
      '36' = '2?',
      '56' = '3?',
      '33' = 'Regular insulin dose ',
      '34' = 'NPH insulin dose ',
      '35' = 'UltraLente insulin dose ',
      '48' = 'Blood glucose measurement ',
      '57' = 'Unspecified blood glucose measurement2',    # 48 = 57 ??
      '58' = 'Pre-breakfast blood glucose measurement ',
      '59' = 'Post-breakfast blood glucose measurement ',
      '60' = 'Pre-lunch blood glucose measurement ',
      '61' = 'Post-lunch blood glucose measurement ',
      '62' = 'Pre-supper blood glucose measurement ',
      '63' = 'Post-supper blood glucose measurement', 
      '64' = 'Pre-snack blood glucose measurement',
      '65' = 'Hypoglycemic symptoms ',
      '66' = 'Typical meal ingestion ',
      '67' = 'More-than-usual meal ingestion' ,
      '68' = 'Less-than-usual meal ingestion' ,
      '69' = 'Typical exercise activity' ,
      '70' = 'More-than-usual exercise activity' ,
      '71' = 'Less-than-usual exercise activity' ,
      '72' = 'Unspecified special event'
    )
    
    # refactoring codes
    df$code_refact <- df$code
    df$code_refact[df$code_refact %in% c(4, 36, 48, 56:64)] <- 48
    index_still_present <- which(names(code_labels) %in% unique(df$code_refact))
    df$code_refact <- factor(df$code_refact, 
                             levels = names(code_labels)[index_still_present], 
                             labels = (as.character(code_labels))[index_still_present])
    
    # pick only specific measures and aggregate by average when double measurements
    tmp <- df[df$code %in% c(33:35, 48, 56:64), c('id', 'timestamp', 'code_refact', 'value')]
    tmp <- aggregate(x = tmp$value, 
                     by = list('id' = tmp$id, 'timestamp' = tmp$timestamp, 'code' = tmp$code_refact),
                     mean, na.rm = TRUE)
  
    # refactor code variable, only 4 instances left 
    tmp$code <- factor(as.character(tmp$code))
    
    # rename column from 'x' to 'value' ('x' is due to using aggregate) 
    colnames(tmp)[4] <- 'value'
    
    # reshape dataset - one line per (id, timestamp) key
    tmp <- reshape(tmp, idvar = c('id', 'timestamp'), timevar = 'code', v.names = 'value', direction = 'wide')
    
    # uniformize column names to avoid problems
    colnames(tmp) <- gsub('\\.| ', '\\_', tolower(colnames(tmp)))
    
    # compute daytime
    tmp$daytime <- as.numeric(format(x = tmp$timestamp, '%H')) + as.numeric(format(x = tmp$timestamp, '%M')) / 60
    
    # sort dataset
    tmp <- tmp[order(tmp$id, tmp$timestamp),]
    
    # add discretized daytime : four periods of the day
    tmp$daytime_disc <- discretize(x = tmp$daytime, method = 'cluster', categories = 4)

    # return dataset
    return(tmp)
    
  }
  
}
