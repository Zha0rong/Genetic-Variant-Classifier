#Read in the Variant Call Form Database
arguments=commandArgs(TRUE)
Database=toString(arguments[1])
library('vcfR')
library('tidyr')
library('Biostrings')
library('biomaRt')
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = c('chromosome_name','start','end')
attributes=c('external_gene_name','gene_biotype','name_1006') #ensembl_gene_id, name_1006
vcf=read.vcfR(Database)
fix=data.frame(vcf@fix,stringsAsFactors = FALSE)
fix=fix[,-c(ncol(fix),ncol(fix)-1,ncol(fix)-2)]
info=data.frame(extract_info_tidy(vcf),stringsAsFactors = FALSE)
largedata=data.frame(cbind(fix,info),stringsAsFactors = FALSE)
rm(vcf,fix,info)
#write.table(largedata,'Database.tbl',quote = FALSE,sep = '\t')
print('Finish loading in the vcf format database.')
#Note on dataset index: 1. Chromosome number
#                       2. Position
#                       3. Reference base
#                       4. Alternative Base
#                       5. Allele ID
#                       6. Disease Associated
#                       7. Clinical Significance
#                       8. secondary Clinical Significance(COnflict interpretation)
#                       9. Tetiary Significance
#                       10. Mutation tye
#                       11. Consequence type
#This for-loop strips off the SOI index (journal-related) from the data
getting_information=function(Database) {
  workingdataset=Database[,c(1,2,4,5,10,11,17,18,19,20,25)]
  workingdataset[,11] = gsub('SO:\\d{7}\\|','',workingdataset[,11])
  workingdataset[,6] = gsub('_',' ',workingdataset[,6])
  workingdataset[,7] = gsub('_',' ',workingdataset[,7])
  return(workingdataset)
}
workingdataset=getting_information(largedata)
print('Finish loading the features and prior')
rm(largedata)

#This chunck of codes will handle the genetic variants as 'Conflicting of interpretation: These are the variants that researchers are not sure whether they are pathogenic or not.
#However, the number of reports of different interpretations are provided
Conflict_testing_data_set=subset(workingdataset,grepl('Conflicting interpretations of pathogenicity',workingdataset$CLNSIG,ignore.case = TRUE))
Training_data_set=subset(workingdataset,!grepl('Conflicting interpretations of pathogenicity',workingdataset$CLNSIG,ignore.case = TRUE))

#This chunck of codes collects the amount of types of variants in the database
Typesofvariant=data.frame(table(unlist(Training_data_set$CLNVC)),stringsAsFactors = FALSE)
colnames(Typesofvariant)[1]='Types'

#This chunck of for loops collects the amount of types of consequences in the database
consequence=data.frame(data.frame(table(Training_data_set[,11]),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
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
disease=data.frame(data.frame(table(unlist(Training_data_set[,6])),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
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
  }
  else {
    if (disease[i,1] %in% diseaseofvariant) {
    }
    else {
      diseaseofvariant[toString(disease[i,1])] = toString(disease[i,1])
    }
  }
}
rm(i,j,listofdisease,disease)
print('Finish loading the types of diseases associated with genetic variants.')
#This Chunk of codes formatted the number of significance.
variant=data.frame(data.frame(table(unlist(Training_data_set[,7])),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
significanceofvariant=list()
for (i in 1:nrow(variant)) {
  if (grepl(', ',variant[i,1])) {
    listofvariant=data.frame(strsplit(toString(variant[i,1]),', '),stringsAsFactors = FALSE)
    colnames(listofvariant)='variants'
    for (j in 1:nrow(listofvariant)) {
      if (grepl('/',listofvariant[j,1])) {
        list=data.frame(strsplit(toString(listofvariant[j,1]),'/'),stringsAsFactors = FALSE)
        colnames(list)='variants'
        for (k in 1:nrow(list)) {
          if (list[k,1] %in% significanceofvariant) {
          }
          else {
            significanceofvariant[toString(list[k,1])] = toString(list[k,1])
          }
        }
      }
      else{
        if (listofvariant[j,1] %in% significanceofvariant) {
        }
        else {
          significanceofvariant[toString(listofvariant[j,1])] = toString(listofvariant[j,1])
        }}
    }
  }
  else {
    if (grepl('/',variant[i,1])) {
      a=data.frame(strsplit(toString(variant[i,1]),'/'),stringsAsFactors = FALSE)
      colnames(a)='variants'
      for (l in 1:nrow(a)) {
        if (a[l,1] %in% significanceofvariant) {
        }
        else {
          significanceofvariant[toString(a[l,1])] = toString(a[l,1])
        }
      }
    }
    else {
      if (variant[i,1] %in% significanceofvariant) {
      }
      else {
        significanceofvariant[toString(variant[i,1])] = toString(variant[i,1])
      }}
    
  }
}
rm(i,j,listofvariant,variant,list,k,l,a)
#This chunck of codes calculates the prior probabilities.
Initial_Prior_Probability=data.frame(table(Training_data_set[,7]),stringsAsFactors = FALSE)
Prior_Probability=data.frame(matrix(significanceofvariant),stringsAsFactors = FALSE)
Prior_Probability$Probability=rep(0,nrow(Prior_Probability))
colnames(Prior_Probability)[1]='Class'
for (i in 1:nrow(Prior_Probability)) {
  for (j in 1:nrow(Initial_Prior_Probability)) {
    if (grepl(Prior_Probability[i,1],Initial_Prior_Probability[j,1],ignore.case = TRUE)) {
      Prior_Probability[i,2] = Prior_Probability[i,2] + Initial_Prior_Probability[j,2]
    }
  }
}
Prior_Probability$Probability=Prior_Probability$Probability/nrow(Training_data_set)
rm(i,j,Initial_Prior_Probability)
#Get the first feature likelihood.
Pathogenic=subset(Training_data_set,grepl('pathogenic|likely pathogenic',Training_data_set$CLNSIG,ignore.case = TRUE))
Benign=subset(Training_data_set,grepl('benign|likely benign',Training_data_set$CLNSIG,ignore.case = TRUE))
#The results fit with the prior probability, ready to proceed.
Feature_type_of_variant=data.frame(Typesofvariant[,1],stringsAsFactors = FALSE)
colnames(Feature_type_of_variant)[1]='Variant'
Feature_type_of_variant$Pathogenic=rep(0,nrow(Feature_type_of_variant))
Feature_type_of_variant$Benign=rep(0,nrow(Feature_type_of_variant))
p=data.frame(table(Pathogenic$CLNVC))
b=data.frame(table(Benign$CLNVC))
for (i in 1:nrow(Feature_type_of_variant)) {
  for (j in 1:nrow(p)) {
    if (grepl(p[j,1],Feature_type_of_variant[i,1],ignore.case = TRUE)) {Feature_type_of_variant[i,2]=p[j,2]}
  }
  for (k in 1:nrow(b)) {
    if (grepl(b[k,1],Feature_type_of_variant[i,1],ignore.case = TRUE)) {Feature_type_of_variant[i,3]=b[k,2]}
  }
}
rm(i,j,k,p,b)
Feature_type_of_variant$Pathogenic=(Feature_type_of_variant$Pathogenic+1)/(nrow(Pathogenic)+1*nrow(Feature_type_of_variant))
Feature_type_of_variant$Benign=(Feature_type_of_variant$Benign+1)/(nrow(Benign)+1*nrow(Feature_type_of_variant))
Feature_type_of_consequence=data.frame(matrix(unlist(consequenceofvariant)),stringsAsFactors = FALSE)
colnames(Feature_type_of_consequence)='Consequence'
Feature_type_of_consequence$Pathogenic=rep(0,nrow(Feature_type_of_consequence))
Feature_type_of_consequence$Benign=rep(0,nrow(Feature_type_of_consequence))
for (i in 1:nrow(Feature_type_of_consequence)) {
  c=Feature_type_of_consequence[i,1]
  for (j in 1:nrow(Pathogenic)) {
    if (grepl(c,Pathogenic$MC[j],ignore.case = TRUE)) {
      Feature_type_of_consequence[i,2] = Feature_type_of_consequence[i,2]+1
    }
  }
  for (k in 1:nrow(Benign)) {
    if (grepl(c,Benign$MC[k],ignore.case = TRUE)) {
      Feature_type_of_consequence[i,3] = Feature_type_of_consequence[i,3]+1
    }
  }
}
Feature_type_of_consequence$Pathogenic=(Feature_type_of_consequence$Pathogenic+1)/(nrow(Pathogenic)+1*nrow(Feature_type_of_consequence))
Feature_type_of_consequence$Benign=(Feature_type_of_consequence$Benign+1)/(nrow(Benign)+1*nrow(Feature_type_of_consequence))
rm(c,i,j,k)


####################################################
Pathogen_ensembel_info=data.frame(getBM(attributes,filters = filters,values = list(Pathogenic[1,1],Pathogenic[1,2],Pathogenic[1,2]),mart = ensembl,uniqueRows = TRUE),stringsAsFactors = FALSE)
Pathogen_ensembel_info[1,1]=paste(toString(Pathogen_ensembel_info[,1]))
Pathogen_ensembel_info[1,2]=paste(toString(Pathogen_ensembel_info[,2]))
Pathogen_ensembel_info[1,3]=paste(toString(Pathogen_ensembel_info[,3]))
Pathogen_ensembel_info=Pathogen_ensembel_info[1,]
for ( i in 1:nrow(Pathogenic)) {
  values= list(Pathogenic[i,1],Pathogenic[i,2],Pathogenic[i,2])
  a=data.frame(getBM(attributes,filters = filters,values = values,mart = ensembl,uniqueRows = TRUE),stringsAsFactors = FALSE)
  if (nrow(a)==1) {
    Pathogen_ensembel_info=rbind(Pathogen_ensembel_info,a)
  }
  else {
    a[1,1]=paste(toString(a[,1]))
    a[1,2]=paste(toString(a[,2]))
    a[1,3]=paste(toString(a[,3]))
    Pathogen_ensembel_info=rbind(Pathogen_ensembel_info,a[1,])
  }
  print(i)
}
rm(i,values,a)
Pathogenic=cbind(Pathogenic,Pathogen_ensembel_info)
#####################################################
Benign_ensembel_info=data.frame(getBM(attributes,filters = filters,values = list(Benign[1,1],Benign[1,2],Benign[1,2]),mart = ensembl,uniqueRows = TRUE),stringsAsFactors = FALSE)
Benign_ensembel_info[1,1]=paste(toString(Benign_ensembel_info[,1]))
Benign_ensembel_info[1,2]=paste(toString(Benign_ensembel_info[,2]))
Benign_ensembel_info[1,3]=paste(toString(Benign_ensembel_info[,3]))
Benign_ensembel_info=Benign_ensembel_info[1,]
for ( i in 1:nrow(Benign)) {
  values= list(Benign[i,1],Benign[i,2],Benign[i,2])
  a=data.frame(getBM(attributes,filters = filters,values = values,mart = ensembl,uniqueRows = TRUE),stringsAsFactors = FALSE)
  if (nrow(a)==1) {
    Benign_ensembel_info=rbind(Benign_ensembel_info,a)
  }
  else {
    a[1,1]=paste(toString(a[,1]))
    a[1,2]=paste(toString(a[,2]))
    a[1,3]=paste(toString(a[,3]))
    Benign_ensembel_info=rbind(Benign_ensembel_info,a[1,])
  }
  print(i)
}
Benign=cbind(Benign,Benign_ensembel_info)
####################################################
