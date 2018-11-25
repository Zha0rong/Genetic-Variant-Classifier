#Read in the Variant Call Form Database
arguments=commandArgs(TRUE)
Database=toString(arguments[1])
library('vcfR')
library('tidyr')
library('BSgenome')
library('GenomicRanges')
library('GenomicFeatures')
library('Biostrings')
UCSCgenome=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = c('chromosome_name','start','end')
attributes=c('external_gene_name','gene_biotype')
print('Finish loading GRCh 38 human genome UCSC version.')
vcf=read.vcfR(Database)
fix=data.frame(vcf@fix,stringsAsFactors = FALSE)
fix=fix[,-c(ncol(fix),ncol(fix)-1,ncol(fix)-2)]
info=data.frame(extract_info_tidy(vcf),stringsAsFactors = FALSE)
largedata=data.frame(cbind(fix,info),stringsAsFactors = FALSE)
rm(vcf,fix,info)
write.table(largedata,'Database.tbl',quote = FALSE,sep = '\t')
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
workingdataset=largedata[,c(1,2,4,5,10,11,17,18,19,20,25)]
workingdataset[,11] = gsub('SO:\\d{7}\\|','',workingdataset[,11])
workingdataset[,6] = gsub('_',' ',workingdataset[,6])
workingdataset[,7] = gsub('_',' ',workingdataset[,7])
write.table(workingdataset,'working.tbl',quote = FALSE,sep = '\t')
print('Finish loading the features and prior')
rm(largedata)

#This part is a duplication, later i will determine which one i will keep
ensembl_info=matrix(attributes,c(1,2))
for (i in 1:nrow(workingdataset)) {
values=c(list(workingdataset[i,1]),list(workingdataset[i,2]),list(workingdataset[i,2]))
ensembl_info = rbind(ensembl_info,matrix(getBM(attributes = attributes,filters = filters,values = values,mart = ensembl),c(1,2)))
print(i/nrow(workingdataset))
}
ensembl_info=data.frame(matrix(unlist(ensembl_info),ncol = 2))
colnames(ensembl_info)=ensembl_info[1,]
ensembl_info=ensembl_info[-1,]








Typesofvariant=data.frame(table(unlist(workingdataset$CLNVC)),stringsAsFactors = FALSE)
colnames(Typesofvariant)[1]='Types'
#This chunck of for loops collects the amount of types of consequences in the database
consequence=data.frame(data.frame(table(workingdataset[,11]),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
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
print('Finish loading the consequences of variants.')
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
variant=data.frame(data.frame(table(unlist(workingdataset[,7])),stringsAsFactors = FALSE)[,1],stringsAsFactors = FALSE)
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
print('Finish loading the clinical significance of genetic variants.')
#This chunck of codes will handle the genetic variants as 'Conflicting of interpretation: These are the variants that researchers are not sure whether they are pathogenic or not.
#However, the number of reports of different interpretations are provided
Conflict_testing_data_set=subset(workingdataset,grepl('Conflicting interpretations of pathogenicity',workingdataset$CLNSIG,ignore.case = TRUE))
Training_data_set=subset(workingdataset,!grepl('Conflicting interpretations of pathogenicity',workingdataset$CLNSIG,ignore.case = TRUE))
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
print('Finish calculating likelihood of type of variant feature.')
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
print('Finish calculating likelihood of consequence feature.')
#Patho=Pathogenic[,c(1,2)]
#Benig=Benign[,c(1,2)]
#P=matrix(attributes,c(1,2))
#for (i in 1:nrow(Patho)) {
  #P = rbind(P,matrix(getBM(attributes = attributes,filters = filters,values = c(list(Patho[i,1]),list(Patho[i,2]),list(Patho[i,2])),mart = ensembl),c(1,2)))
  #print(i/nrow(Patho))
#}
#P=data.frame(matrix(unlist(P),ncol = 2))
#colnames(P)=P[1,]
#P=P[-1,]
#B=matrix(attributes,c(1,2))
#for (i in 1:nrow(Benig)) {
  #B = rbind(B,matrix(getBM(attributes = attributes,filters = filters,values = c(list(Benig[i,1]),list(Benig[i,2]),list(Benig[i,2])),mart = ensembl),c(1,2)))
  #print(i/nrow(Benig))
#}
#B=data.frame(matrix(unlist(B),ncol = 2))
#colnames(B)=B[1,]
#B=B[-1,]
#rm(i)
