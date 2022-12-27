#The script below will take the standard mothur taxonomy file as an input and create a
#dataframe with separate variables for each level (kingdom, phylum, class, etc.)
#replace with the name of the taxonomy mothur output file you are using

# Input: filename
library(stringr)
library(dplyr)

mothur_tax_munging_func = function(filename){

taxonomy<-read.table(file = filename, header = T, sep = "\t")

taxonomy$OTU<-as.character(taxonomy$OTU)
taxonomy$Taxonomy<-as.character(taxonomy$Taxonomy)
taxonomy$taxonomy<-taxonomy$Taxonomy

taxonomy$kingdom<-str_sub(taxonomy$Taxonomy,1, str_locate(taxonomy$Taxonomy, pattern = ";")[,1])
taxonomy$Taxonomy<-str_sub(taxonomy$Taxonomy, str_locate(taxonomy$Taxonomy, pattern = ";")[,1]+1, nchar(taxonomy$Taxonomy))
taxonomy$kingdom<-str_sub(taxonomy$kingdom, 1, str_locate(taxonomy$kingdom, pattern = "\\(")[,1]-1)

taxonomy$phylum<-str_sub(taxonomy$Taxonomy,1, str_locate(taxonomy$Taxonomy, pattern = ";")[,1])
taxonomy$Taxonomy<-str_sub(taxonomy$Taxonomy, str_locate(taxonomy$Taxonomy, pattern = ";")[,1]+1, nchar(taxonomy$Taxonomy))
taxonomy$phylum<-str_sub(taxonomy$phylum, 1, str_locate(taxonomy$phylum, pattern = "\\(")[,1]-1)


taxonomy$class<-str_sub(taxonomy$Taxonomy,1, str_locate(taxonomy$Taxonomy, pattern = ";")[,1])
taxonomy$Taxonomy<-str_sub(taxonomy$Taxonomy, str_locate(taxonomy$Taxonomy, pattern = ";")[,1]+1, nchar(taxonomy$Taxonomy))
taxonomy$class<-str_sub(taxonomy$class, 1, str_locate(taxonomy$class, pattern = "\\(")[,1]-1)


taxonomy$order<-str_sub(taxonomy$Taxonomy,1, str_locate(taxonomy$Taxonomy, pattern = ";")[,1])
taxonomy$Taxonomy<-str_sub(taxonomy$Taxonomy, str_locate(taxonomy$Taxonomy, pattern = ";")[,1]+1, nchar(taxonomy$Taxonomy))
taxonomy$order<-str_sub(taxonomy$order, 1, str_locate(taxonomy$order, pattern = "\\(")[,1]-1)

taxonomy$family<-str_sub(taxonomy$Taxonomy,1, str_locate(taxonomy$Taxonomy, pattern = ";")[,1])
taxonomy$Taxonomy<-str_sub(taxonomy$Taxonomy, str_locate(taxonomy$Taxonomy, pattern = ";")[,1]+1, nchar(taxonomy$Taxonomy))
taxonomy$family<-str_sub(taxonomy$family, 1, str_locate(taxonomy$family, pattern = "\\(")[,1]-1)

taxonomy$genus<-str_sub(taxonomy$Taxonomy,1, str_locate(taxonomy$Taxonomy, pattern = ";")[,1])
taxonomy$Taxonomy<-str_sub(taxonomy$Taxonomy, str_locate(taxonomy$Taxonomy, pattern = ";")[,1]+1, nchar(taxonomy$Taxonomy))
taxonomy$genus<-str_sub(taxonomy$genus, 1, str_locate(taxonomy$genus, pattern = "\\(")[,1]-1)

taxonomy<-dplyr::select(taxonomy, -Taxonomy)
return(taxonomy)

}

# pre-filtering (MixOmics processing)
low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent=0.01 # cutoff chosen
){
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}
