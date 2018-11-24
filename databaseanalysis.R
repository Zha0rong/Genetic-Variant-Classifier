#Read in the Variant Call Form Database
arguments=commandArgs(TRUE)
Database=toString(arguments[1])
library('vcfR')
library('tidyr')
vcf=read.vcfR('clinvar.vcf')
fix=data.frame(vcf@fix,stringsAsFactors = FALSE)
fix=fix[,-c(ncol(fix),ncol(fix)-1,ncol(fix)-2)]
info=data.frame(extract_info_tidy(vcf),stringsAsFactors = FALSE)
largedata=data.frame(cbind(fix,info),stringsAsFactors = FALSE)
rm(vcf,fix,info)
#Note on dataset index: 1. Chromosome number
#                       2. Position
#                       3. Reference base
#                       4. Alternative Base
#                       5. Allele ID
#                       6. Disease Associated
#                       7. Clinical Significance
#                       8. secondary Clinical Significance
#                       9. Tetiary Significance
#                       10. Mutation tye
#                       11. Consequence type
#This for-loop strips off the SOI index (journal-related) from the data
workingdataset=largedata[,c(1,2,4,5,10,11,17,18,19,20,25)]
for (i in 1:nrow(workingdataset)) {
  if(!is.na(workingdataset[i,11])) {
  workingdataset[i,11] = gsub('SO:\\d{7}\\|','',workingdataset[i,11])
  }
  print(i/nrow(workingdataset))
}
#This chunck of for loops collects the amount of types of consequences in the database
consequence=data.frame(data.frame(table(unlist(workingdataset[,11])),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
consequenceofvariant=list()
for (i in 1:nrow(consequence)) {
  if (grepl(',',consequence[i,1])) {
    listofconsequence=data.frame(strsplit(toString(consequence[i,1]),','),stringsAsFactors = FALSE)
    colnames(listofconsequence)='consequences'
    for (j in 1:nrow(listofconsequence)) {
      if (listofconsequence[j,1] %in% consequenceofvariant) {
      }
      else {
        consequenceofvariant[toString(listofconsequence[j,1])] = toString(listofconsequence[j,1])
      }
    }
  }
  else {
    if (consequence[i,1] %in% consequenceofvariant) {
    }
    else {
      consequenceofvariant[toString(consequence[i,1])] = toString(consequence[i,1])
    }
    
  }
}
rm(i,j,listofconsequence,consequence)
#This chunck of for loops collects the amount of types of diseases in the database
disease=data.frame(data.frame(table(unlist(workingdataset[,6])),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
diseaseofvariant=list()
for (i in 1:nrow(disease)) {
  if (grepl('\\|',disease[i,1])) {
    listofdisease=data.frame(strsplit(toString(disease[i,1]),'\\|'),stringsAsFactors = FALSE)
    colnames(listofdisease)='diseases'
    for (j in 1:nrow(listofdisease)) {
      if (listofdisease[j,1] %in% diseaseofvariant) {
      }
      else {
        diseaseofvariant[toString(listofdisease[j,1])] = toString(listofdisease[j,1])
      }
    }
    print(i/nrow(disease))
    
  }
  else {
    if (disease[i,1] %in% diseaseofvariant) {
    }
    else {
      diseaseofvariant[toString(disease[i,1])] = toString(disease[i,1])
    }
    print(i/nrow(disease))
    
  }
}
rm(i,j,listofdisease,disease)