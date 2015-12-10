

http://services.nacse.org/prism/data/public/4km/tmin/20090405

wget --content-disposition http://services.nacse.org/prism/data/public/4km/tmin/20090405


st <- as.Date(start_date,format="%Y.%m.%d") #start date
en <- as.Date(end_date,format="%Y.%m.%d") #end date
ll <- seq.Date(st, en, by="1 day") #sequence of dates
dates_queried <- format(ll,"%Y.%m.%d") #formatting queried dates

url_product <-paste("http://e4ftl01.cr.usgs.gov/MOLT/",MODIS_product,"/",sep="") #URL is a constant...
#url_product <- file.path("http://e4ftl01.cr.usgs.gov/MOLT/",MODIS_product)

dates_available <- extractFolders(url_product)  #Get the list of available driectory dates for the product, from 2000 to now

list_folder_dates <- intersect(as.character(dates_queried), as.character(dates_available)) #list of remote folders to access
#list_folder_dates <-setdiff(as.character(dates_available), as.character(dates_queried))

#step 2: list content of specific day folder to obtain specific file...  #parse by tile name!!!

url_folders_str <-paste(url_product,list_folder_dates,"/",sep="") #url for the folders matching dates to download

## loop over
#debug(extractFiles)
#generate outdir for each tile!!!

list_folders_files <- vector("list",length(url_folders_str))
d_files <- vector("list",length(list_folders_files))
file_format<-c("hdf","xml")
#can make this faster using parallelization...
for (i in 1:length(url_folders_str)){
  #debug(extractFiles)
  list_folders_files[[i]] <- try(extractFiles(url_folders_str[i], list_tiles)[file_format]) 
  #d_files[[i]] <- list_folders_files[[i]][[file_format]]                      
}
#list_folders_files <- lapply(url_folders_str,extractFiles,list_tiles_str=list_tiles)
#Now remove error objects...
d_files_tmp <-remove_from_list_fun(l_x=list_folders_files,condition_class ="try-error")
d_files <- as.character(unlist(d_files_tmp$list)) #all the files to download...
n_dir_available <- length(d_files)
n_dir_requested <- length(list_folders_files)

#Step 3: download file to the directory 
#browser()
#prepare files and directories for download
out_dir_tiles <- file.path(out_dir,list_tiles)
list_files_tiles <- vector("list",length(list_tiles))
for(j in 1:length(out_dir_tiles)){
  if (!file.exists(out_dir_tiles[j])){
    dir.create(out_dir_tiles[j])
  }
  list_files_tiles[[j]] <- grep(pattern=list_tiles[j],x=d_files,value=TRUE) 
}

#Now download per tiles: can be parallelized
for (j in 1:length(list_files_tiles)){ #loop around tils
  file_items <- list_files_tiles[[j]]
  for (i in 1:length(file_items)){
    file_item <- file_items[i]
    download.file(file_item,destfile=file.path(out_dir_tiles[j],basename(file_item)))
    #download.file(file_item,destfile="test.hdf")
  }
}

#Prepare return object: list of files downloaded with http and list downloaded of files in tiles directories

list_files_by_tiles <-mapply(1:length(out_dir_tiles),FUN=function(i,x){list.files(path=x[[i]],pattern="*.hdf$",full.names=T)},MoreArgs=(list(x=out_dir_tiles))) #Use mapply to pass multiple arguments
#list_files_by_tiles <-mapply(1:length(out_dir_tiles),FUN=list.files,MoreArgs=list(pattern="*.hdf$",path=out_dir_tiles,full.names=T)) #Use mapply to pass multiple arguments

colnames(list_files_by_tiles) <- list_tiles #note that the output of mapply is a matrix
download_modis_obj <- list(list_files_tiles,list_files_by_tiles)
names(download_modis_obj) <- c("downloaded_files","list_files_by_tiles")
return(download_modis_obj)


