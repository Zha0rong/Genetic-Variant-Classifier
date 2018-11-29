#Read in the Variant Call Form Database
######Prep Work######
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

######This chunck of codes will handle the genetic variants as 'Conflicting of interpretation: These are the variants that researchers are not sure whether they are pathogenic or not. However, the number of reports of different interpretations are provided.#######
Conflict_testing_data_set=subset(workingdataset,grepl('Conflicting interpretations of pathogenicity',workingdataset$CLNSIG,ignore.case = TRUE))
Training_data_set=subset(workingdataset,!grepl('Conflicting interpretations of pathogenicity',workingdataset$CLNSIG,ignore.case = TRUE))
Pathogenic=subset(Training_data_set,grepl('pathogenic|likely pathogenic',Training_data_set$CLNSIG,ignore.case = TRUE))
Benign=subset(Training_data_set,grepl('benign|likely benign',Training_data_set$CLNSIG,ignore.case = TRUE))
######This chunck of codes collects the amount of types of variants in the database######
Typesofvariant=data.frame(table(unlist(Training_data_set$CLNVC)),stringsAsFactors = FALSE)
colnames(Typesofvariant)[1]='Types'
######Collects types of variants, consequences and diseases######
Types_collectors = function(Data,column,sep) {
  data=data.frame(data.frame(table(Data[,column]),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
  result=list()
  for (i in 1:nrow(data)) {
    if (grepl(sep,data[i,1])) {
      listofdata=data.frame(strsplit(toString(data[i,1]),split = sep),stringsAsFactors = FALSE)
      colnames(listofdata)='datas'
      for (j in 1:nrow(listofdata)) {
        if (listofdata[j,1] %in% result) {
        }
        else {
          result[toString(listofdata[j,1])] = toString(listofdata[j,1])
        }
      }
    }
    else {
      if (data[i,1] %in% result) {
      }
      else {
        result[toString(data[i,1])] = toString(data[i,1])
      }
    }
  }
  return(result)
}
Types_collectors_multi_sep = function(Data,column) {
  variant=data.frame(data.frame(table(unlist(Data[,column])),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
  significanceofvariant=list()
  for (i in 1:nrow(variant)) {
    if (grepl(',',variant[i,1])) {
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
  return(significanceofvariant)
}
consequenceofvariant=Types_collectors(Training_data_set,11,',')
diseaseofvariant=Types_collectors(Training_data_set,6,'\\|')
significanceofvariant=Types_collectors_multi_sep(Training_data_set,7)
Typesofvariant=data.frame(table(unlist(Training_data_set$CLNVC)),stringsAsFactors = FALSE)
colnames(Typesofvariant)[1]='Types'

######This chunck of codes calculates the prior probabilities.#######
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
######Get the features' likelihood.######
feature_calculator = function (Pathogenic,Benign,column,feature,output) {
  result=data.frame(matrix(unlist(feature)),stringsAsFactors = FALSE)
  colnames(result)=output
  result$Pathogenic=rep(0,nrow(result))
  result$Benign=rep(0,nrow(result))
  for (i in 1:nrow(result)) {
    c=result[i,1]
    for (j in 1:nrow(Pathogenic)) {
      if (grepl(c,toString(Pathogenic[j,column]),ignore.case = TRUE)) {
        result[i,2] = result[i,2]+1
      }
    }
    for (k in 1:nrow(Benign)) {
      if (grepl(c,toString(Benign[k,column]),ignore.case = TRUE)) {
        result[i,3] = result[i,3]+1
      }
    }
  }
  result$Pathogenic=(result$Pathogenic+1)/(nrow(Pathogenic)+1*nrow(result))
  result$Benign=(result$Benign+1)/(nrow(Benign)+1*nrow(result))
  return(result)
}
Feature_type_of_variant=feature_calculator(Pathogenic,Benign,10,Typesofvariant[1],'types')
Feature_type_of_consequence=feature_calculator(Pathogenic,Benign,10,consequenceofvariant,'consequence')
######Ensembl Information retrieve, do not use now.######
ensembl_info_retrieve=function(Data,ensembl,filters,attributes){
  Result=data.frame(getBM(attributes,filters = filters,values = list(Data[1,1],Data[1,2],Data[1,2]),mart = ensembl,uniqueRows = TRUE),stringsAsFactors = FALSE)
  Result[1,1]=paste(toString(Result[,1]))
  Result[1,2]=paste(toString(Result[,2]))
  Result[1,3]=paste(toString(Result[,3]))
  Result=Result[1,]
  for ( i in 1:nrow(Data)) {
    values= list(Data[i,1],Data[i,2],Data[i,2])
    a=data.frame(getBM(attributes,filters = filters,values = values,mart = ensembl,uniqueRows = TRUE),stringsAsFactors = FALSE)
    if (nrow(a)==1) {
      Result=rbind(Result,a)
    }
    else {
      a[1,1]=paste(toString(a[,1]))
      a[1,2]=paste(toString(a[,2]))
      a[1,3]=paste(toString(a[,3]))
      Result=rbind(Result,a[1,])
    }
  }
  rm(i,values,a)
  Data=cbind(Data,Result)
  return(Data)
}
Pathogenic=ensembl_info_retrieve(Pathogenic,ensembl,filters,attributes)
Benign=ensembl_info_retrieve(Benign,ensembl,filters,attributes)
#Optional Training_data_set=ensembl_info_retrieve(Training_data_set,ensembl,filters,attributes)
Slimer=function(Data,column) {
  Data=Data
  for (i in 1:nrow(Data)) {
    a=data.frame(data.frame(table(unlist(strsplit(Data[i,column],', '))),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
    if(nrow(a) > 1) {
      b=''
      for (i in 1:nrow(a)) {
        b=paste(b,a[i,1],',')
      }
      Data[i,column]=b

    }
    else{Data[i,column]=toString(a[1,1])}
  }
  return(Data)
}
######Continue to get features' likelihood######
#If you have run 'Training_data_set=ensembl_info_retrieve(Training_data_set,ensembl,filters,attributes)'
# Gene_type_of_Variant = Types_collectors(Training_data_set,13,', ')
Gene_type_of_Variant_Pathogenic=Types_collectors(Pathogenic,13,', ')
Gene_type_of_Variant_Benign=Types_collectors(Benign,13,', ')
Gene_type_of_Variant=list()
for (i in 1:length(Gene_type_of_Variant_Benign)) {
  a=toString(Gene_type_of_Variant_Benign[i])
  if (a %in% Gene_type_of_Variant) {}
  else {Gene_type_of_Variant[a]=a}
}
for (i in 1:length(Gene_type_of_Variant_Pathogenic)) {
  a=toString(Gene_type_of_Variant_Pathogenic[i])
  if (a %in% Gene_type_of_Variant) {}
  else {Gene_type_of_Variant[a]=a}
}
rm(Gene_type_of_Variant_Benign,Gene_type_of_Variant_Pathogenic)
Feature_gene_type_of_variant=feature_calculator(Pathogenic,Benign,13,Gene_type_of_Variant,'Gene_type')
