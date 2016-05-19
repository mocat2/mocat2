
# This script will start the shiny files located within the same folder

# define load function
cat ('Loading functions...\n')
load.function <- function(x) { 
  x <- as.character(substitute(x)) 
  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    cat(paste('Loading', x, '...'), sep='')
    eval(parse(text=paste("require(", x, ")", sep=""))) 
    cat (' OK!\n')
  } else {
    cat('Updating packages...')
    update.packages(repos="http://mirrors.softliste.de/cran/", ask=FALSE)    
    cat (' OK!\n')
    cat('Installing packages...')
    eval(parse(text=paste("install.packages('", x, "', repos=\"http://mirrors.softliste.de/cran/\")", sep=""))) 
    cat (' OK!\n')
    cat(paste('Loading', x, '...'), sep='')
    eval(parse(text=paste("require(", x, ")", sep="")))
    cat (' OK!\n')
  } 
}

# Prepare Shiny
load.function('shiny')

# Run shiny
