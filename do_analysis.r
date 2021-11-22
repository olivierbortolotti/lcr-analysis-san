library("stringr")
do_analysis<-function(file_name){
  file_name<<-str_sub(file_name,1,-5)
  source("sparks analysis.r")
}
list_datasets<-list.files("~/sparks-analysis/data")

lapply(list_datasets[1:length(list_datasets)-1], do_analysis )
