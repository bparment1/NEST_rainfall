###Loading R library and packages                                                      
#library(gtools)    # loading some useful tools 

library(sp)
library(raster)
library(rgdal)
require(rgeos)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
require(RCurl)
require(stringr)
require(XML)

function(start_date, end_date,var_name,num_cores){
  #
  download_prism <- function(date_selected,var_name){
    url_product <- "http://services.nacse.org/prism/data/public/4km/" #URL is a constant...
    var_name <- "tmin" #tmax,ppt
    url_product_var <- paste(url_product,var_name,sep="")
    #system("wget --content-disposition $base_url/$clim_var/$day");
    cmd_str <- paste("wget --content-disposition ",url_product_var,date_selected,sep="")
    df <- data.frame(var=var_name,date=date_selected,url=url_product)
    return(df)
  }
  
  start_date <- as.Date(start_date,format="%Y/%m/%d") #start date
  end_date <- as.Date(end_date,format="%Y/%m/%d") #end date
  dates_range <- seq.Date(start_date, end_date, by="1 day") #sequence of dates
  
  l_df_download <- mclapply(dates_range,FUN=download_prism,var_name=var_name,
                            mc.preschedule=FALSE,mc.cores = num_cores)
  
  #wget http://services.nacse.org/prism/data/public/4km/tmin/20090405
  #
  #
}

#################### END OF SCRIPT #####################

# use strict;
# use warnings;
# use DateTime;
# day <- 
# my $clim_var = 'ppt';
# my $base_url = 'http://services.nacse.org/prism/data/public/4km';
# my $start = DateTime->new( day => 1, month => 10, year => 1999 );
# my $stop = DateTime->new( day => 30, month => 9, year => 2000 );
# while($start <= $stop) {
#   $day = $start->strftime('%Y%m%d'); #place date in proper form
#   system("wget --content-disposition $base_url/$clim_var/$day");
#   sleep 2; #to be nice to our server
#   $start->add(days => 1);
# }
# 
# 
# dates_queried <- format(ll,"%Y.%m.%d") #formatting queried dates
# 



