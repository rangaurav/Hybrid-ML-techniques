library(Peptides)
library(peptider)
library(seqinr)
library(Biostrings)


sequence<-readDNAStringSet("proteinsequences.fasta")
i=1
for (i in 1:length(sequence)){
  
  cat("\n Completed: ",i,"/", length(sequence))
  
  peptide         <- paste(unlist(sequence[i]), collapse='')
  
  ############################################
  #     Part I: Single Properties
  ############################################
  
  # F1: aliphaticIndex
  F1_aliphaticIndex  <- aIndex(peptide)
  
  # F2: bomanIndex
  F2_bomanIndex      <- boman(peptide)

  # F3: instaIndex
  F3_instaIndex      <- instaIndex(peptide)
  
  # F4: probabilityDetectionPeptide
  F4_probabilityDetectionPeptide <- ppeptide(peptide,libscheme = "NNK", N=10^8)
  
  # F5: numberNeighbors
  #F5_numberNeighbors <- getNofNeighbors(peptide,blosum = 1,method = "peptide", libscheme = "NNK")

  # Merging of Part 1 results
  resultPart1=data.frame(F1_aliphaticIndex,F2_bomanIndex,F3_instaIndex,F4_probabilityDetectionPeptide)#,F5_numberNeighbors)
  
  ############################################
  #     Part II: Double Properties
  ############################################
  
  # F6: homentIndex
  F6_homentIndex1     <- hmoment(seq = peptide, angle = 100, window = 11)
  F6_homentIndex2     <- hmoment(seq = peptide, angle = 160, window = 11)

  # F7: molecularWeight
  F7_molecularWeight1 <-mw(seq = peptide,monoisotopic = FALSE)
  F7_molecularWeight2 <-mw(seq = peptide,monoisotopic = TRUE)

  # Merging of Part 2 results
  resultPart2=data.frame(F6_homentIndex1, F6_homentIndex2,F7_molecularWeight1,F7_molecularWeight2)
  
  
  ############################################
  #     Part III: Multiple Properties
  ############################################
  
  # F8: peptideCharge
  pKscale=c("Bjellqvist", "Dawson", "EMBOSS", "Lehninger", "Murray", "Rodwell", "Sillero", "Solomon", "Stryer")
  F8_peptideCharge=c()
  for (j in 1:length(pKscale)){
    x=charge(seq= peptide,pH= seq(from = 5,to = 9,by = 1), pKscale= pKscale[j])
    F8_peptideCharge = c(F8_peptideCharge,x)
  }
  names(F8_peptideCharge)<-paste("F8_pCharge",c(1:length(F8_peptideCharge)),sep='')
  

  # F9: Hydrophobibity for 44 scales
  scale=c("Aboderin", "AbrahamLeo", "Argos", "BlackMould", "BullBreese", "Casari", "Chothia", "Cid", "Cowan3.4", "Cowan7.5", "Eisenberg", "Engelman", "Fasman", "Fauchere", "Goldsack", "Guy", "HoppWoods", "Janin", "Jones", "Juretic", "Kidera", "Kuhn", "KyteDoolittle", "Levitt", "Manavalan", "Miyazawa", "Parker", "Ponnuswamy", "Prabhakaran", "Rao", "Rose", "Roseman", "Sweet", "Tanford", "Welling", "Wilson", "Wolfenden", "Zimmerman", "interfaceScale_pH8", "interfaceScale_pH2", "octanolScale_pH8", "octanolScale_pH2", "oiScale_pH8","oiScale_pH2")
  F9_hydro=c()
  for (j in 1:length(scale)){
    x=  hydrophobicity(seq = peptide,scale = scale[j])
    F9_hydro=c(F9_hydro,x)
  }
  names(F9_hydro)<-paste("F9_hydro",c(1:length(F9_hydro)),sep='')

  
  # F10: isoElectricPoint at 9 pKscale
  pKscale=c("Bjellqvist", "EMBOSS", "Murray", "Sillero", "Solomon", "Stryer", "Lehninger", "Dawson","Rodwell")
  F10_isoElectricPoint=c()
  for (j in 1:length(pKscale)){
    x=  pI(peptide, pKscale = pKscale[j])
    F10_isoElectricPoint=c(F10_isoElectricPoint,x)
  }
  names(F10_isoElectricPoint)<-paste("F10_isoEP",c(1:length(F10_isoElectricPoint)),sep='')
  

  # F11: kideraFactors
  F11_kFactors       <- as.numeric(unlist(kideraFactors(seq = peptide)))
  names(F11_kFactors)<-paste("F11_kFactors",c(1:length(F11_kFactors)),sep='')


  # F12: aaComp
  F12_aaComp          <- as.numeric(unlist(aaComp(peptide)))
  names(F12_aaComp)   <-paste("F12_aaComp",c(1:length(F12_aaComp)),sep='')

  
  # F13: aaDescriptors: Mean, SD, Var
  aaDescriptors   <- c(as.numeric(unlist(aaDescriptors(peptide))))
  aaDescriptors   <- matrix(aaDescriptors, nrow = nchar(peptide), byrow = TRUE)

  F13.1_aaDescriptorsMean<-colMeans(aaDescriptors, na.rm = FALSE, dims = 1)
  names(F13.1_aaDescriptorsMean)<-paste("F13.1_aaDescriptorsMean",c(1:length(F13.1_aaDescriptorsMean)),sep='')

  F13.2_aaDescriptorsSd <-sapply(as.data.frame(aaDescriptors), sd)
  names(F13.2_aaDescriptorsSd)<-paste("F13.2_aaDescriptorsSd",c(1:length(F13.2_aaDescriptorsSd)),sep='')

  F13.3_aaDescriptorsVar<-sapply(as.data.frame(aaDescriptors), var)
  names(F13.3_aaDescriptorsVar)<-paste("F13.3_aaDescriptorsVar",c(1:length(F13.3_aaDescriptorsVar)),sep='')

  # Merging of Part 2 results
  resultPart3=c(F8_peptideCharge,F9_hydro,F10_isoElectricPoint,F11_kFactors,
                F12_aaComp,F13.1_aaDescriptorsMean,F13.2_aaDescriptorsSd,
                F13.3_aaDescriptorsVar)

  finalResult = c(resultPart1,resultPart2,resultPart3)
  
  if(i==1){
    write.table(finalResult, "dataSet.csv", sep = ",", row.names=F,col.names = T)
  }
  else{
    write.table(finalResult, "dataSet.csv", sep = ",", row.names=F, col.names = F, append = T)
  }

}

cat("\n Done")