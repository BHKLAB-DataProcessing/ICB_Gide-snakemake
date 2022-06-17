library(data.table)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")

#############################################################################
#############################################################################

expr = read.table( file.path(input_dir, "mel_gide19_exp_data.csv") , sep="," , header=TRUE , stringsAsFactors = FALSE )
rownames(expr) = expr[ , 1 ]
expr = expr[ , -1 ]

#############################################################################
#############################################################################
## Get Clinical data

clinical = read.table( file.path(input_dir, "CLIN_GIDE.txt") , sep="\t" , header=TRUE , stringsAsFactors = FALSE)
clinical = clinical[ clinical$RNA.Sequencing %in% c( "PRE" , "PRE and EDT" ) , ]

clinical$Overall.Survival..Days. = as.numeric( as.character( clinical$Overall.Survival..Days. ) ) /30
clinical$Last.Followup.Status = ifelse( clinical$Last.Followup.Status %in% "Alive" , 0 , 
								ifelse( clinical$Last.Followup.Status %in% c( "Dead" , "Dead, melanoma" ) , 1 , 0 ) )

clinical = clinical[ order( clinical$Overall.Survival..Days. , clinical$Last.Followup.Status ) , c("Patient.no." , "Age..Years." , "Sex" , "Treatment" , "Best.RECIST.response" , "Progression.Free.Survival..Days." , "Overall.Survival..Days." , "Progressed" , "Last.Followup.Status") ]

#############################################################################
#############################################################################

surv = read.table( file.path(input_dir, "mel_gide19_survival_data.csv") , sep="," , header=TRUE , stringsAsFactors = FALSE , dec=",")
rownames(surv) = surv[ , "X" ]

clin = read.table( file.path(input_dir, "mel_gide19_clin_data.csv") , sep="," , header=TRUE , stringsAsFactors = FALSE , dec=",")
clin = clin[ clin[ , "treatment" ] %in% "pre" , ]

clin = as.data.frame( cbind( clin[ , "X" ] , clin[ , "response" ] , "PD-1/PD-L1" , "Melanoma" , NA , NA , NA , NA , NA , NA , NA , NA , NA ) )
colnames(clin) = c( "patient" , "recist" , "drug_type" , "primary" , "age" , "histo" , "response" , "pfs" ,"os" , "t.pfs" , "t.os" , "stage" , "sex" )
rownames(clin) = clin$patient

clin$t.os = as.numeric( as.character( surv[ rownames(clin) , "OS" ] ) )
clin$os = as.numeric( as.character( surv[ rownames(clin) , "status" ] ) )

clin = clin[ order( clin$t.os , clin$os ) , c("patient" , "t.os" , "os") ]

clinical$Patient.no. = clin$patient
clinical = as.data.frame( cbind( clinical , "Melanoma" , NA , NA , NA , NA , NA , NA ) )
colnames(clinical) = c( "patient" , "age" , "sex" , "drug_type" , "recist" , "t.pfs" , "t.os" , "pfs" ,"os" , "primary" , "histo" , "response" , "stage" , "response.other.info" , "dna" , "rna" )
clinical$drug_type = "PD-1/PD-L1"

clinical$t.pfs = as.numeric( as.character( clinical$t.pfs ) ) /30
clinical$pfs = as.numeric( as.character( ifelse( clinical$pfs %in% "Yes" , 1 , ifelse( clinical$pfs %in% "No" , 0 , NA ) ) ) )

clinical$recist = as.character( clinical$recist )

clinical$response = Get_Response( data = clinical )
clinical$rna = "tpm"

clinical = clinical[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

#############################################################################
#############################################################################

clin = clinical
rownames(clin) = clin$patient

patient = intersect( colnames(expr) , rownames(clin) )
clin = clin[ patient , ]
expr =  expr[ , patient ]

case = cbind( patient , 0 , 0 , 1 )
colnames(case ) = c( "patient" , "snv" , "cna" , "expr" )

write.table( case , file = file.path(output_dir, "cased_sequenced.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( clin , file = file.path(output_dir, "CLIN.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( expr , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )

