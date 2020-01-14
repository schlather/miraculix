# This file has been created automatically by 'rfGenerateConstants'


 ## from  ../../RandomFieldsUtils/RandomFieldsUtils/src/AutoRandomFieldsUtils.h

 MAXUNITS 	<- as.integer(4)
 MAXCHAR 	<- as.integer(18)
 RFOPTIONS 	<- "RFoptions"
 isGLOBAL 	<- as.integer(NA)

 WARN_UNKNOWN_OPTION_ALL 	<- as.integer(3)
 WARN_UNKNOWN_OPTION_SINGLE 	<- as.integer(2)
 WARN_UNKNOWN_OPTION_CAPITAL 	<- as.integer(1)
 WARN_UNKNOWN_OPTION_NONE 	<- as.integer(0)
 WARN_UNKNOWN_OPTION 	<- as.integer(10000)
 WARN_UNKNOWN_OPTION_CONDSINGLE 	<- as.integer((WARN_UNKNOWN_OPTION_SINGLE-WARN_UNKNOWN_OPTION))
 WARN_UNKNOWN_OPTION_DEFAULT 	<- as.integer(WARN_UNKNOWN_OPTION_ALL)



 ## from  src/AutoMiraculix.h

 AutoMiraculix_H 	<- as.integer(1)

 AutoCoding 	<- as.integer(0)
 NoSNPcodingR 	<- as.integer(1)
 NoSNPcodingAVX 	<- as.integer(2)
 NoSNPcoding 	<- as.integer(3)
 ThreeBit 	<- as.integer(4)
 Hamming2 	<- as.integer(5)
 Hamming3 	<- as.integer(6)
 Shuffle 	<- as.integer(7)
 Shuffle256 	<- as.integer(8)
 TwoBit 	<- as.integer(9)
 Packed 	<- as.integer(10)
 Packed256 	<- as.integer(11)
 Multiply 	<- as.integer(12)
 Multiply256 	<- as.integer(13)
 CaseCount 	<- as.integer(14)
 unused14 	<- as.integer(15)
 unused15 	<- as.integer(16)
 unused16 	<- as.integer(17)
 unused17 	<- as.integer(18)
 unused18 	<- as.integer(19)
 unused19 	<- as.integer(20)
 unused20 	<- as.integer(21)
 unused21 	<- as.integer(22)
 unused22 	<- as.integer(23)
 unused23 	<- as.integer(24)
 unused24 	<- as.integer(25)
 unused25 	<- as.integer(26)
 unused26 	<- as.integer(27)
 unused27 	<- as.integer(28)
 unused28 	<- as.integer(29)
 unused29 	<- as.integer(30)
 Haplo 	<- as.integer(31)
 UnknownSNPcoding 	<- as.integer(32)


 nr_snpcoding 	<- as.integer((UnknownSNPcoding+1))
 FirstMoBPSmethod 	<- as.integer(Shuffle)
 LastMoBPSmethod 	<- as.integer(Multiply256)
 FirstGenuineMethod 	<- as.integer(NoSNPcoding)
 LastGenuineMethod 	<- as.integer(LastMoBPSmethod)

 CURRENT_VERSION 	<- as.integer(2)

 VERSION 	<- as.integer(0)
 SNPS 	<- as.integer(1)
 INDIVIDUALS 	<- as.integer(2)
 ADDR0 	<- as.integer(3)
 ADDR1 	<- as.integer(4)
 ALIGNADDR0 	<- as.integer(5)
 ALIGNADDR1 	<- as.integer(6)
 SUMGENO 	<- as.integer(7)
 SUMGENO_E9 	<- as.integer(8)
 METHOD 	<- as.integer(9)
 ALIGNMENT 	<- as.integer(10)
 SNPxIND 	<- as.integer(11)
 BITSPERCODE 	<- as.integer(12)
 BYTESPERBLOCK 	<- as.integer(13)
 CODESPERBLOCK 	<- as.integer(14)
 HEADER 	<- as.integer(15)
 DOUBLEINDIV 	<- as.integer(16)
 LEADINGCOL 	<- as.integer(17)
 MEMinUNITS0 	<- as.integer(18)

 MEMinUNITS1 	<- as.integer(19)
 ALIGNEDUNITS0 	<- as.integer(20)
 ALIGNEDUNITS1 	<- as.integer(21)

 UNITSPERINDIV 	<- as.integer(22)

 INFO_GENUINELY_LAST 	<- as.integer(UNITSPERINDIV)
 INFO_LAST 	<- as.integer(63)

 CURRENT_SNPS 	<- as.integer(HEADER)

 GENOMICMATRIX 	<- "genomicmatrix"
 HAPLOMATRIX 	<- "haplomatrix"
 ORIGINVECTOR 	<- "origindata"



 ## from  src/AutoMiraculix.cc

 SNPCODING_NAMES <-
c( "AutoCoding","NoSNPcodingR","NoSNPCodingAVX","NoSNPcoding","ThreeBit","Hamming2","Hamming3","Shuffle","Shuffle256","TwoBit","Packed","Packed256","Multiply","Multiply256","CaseCount","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","Haplo","unknown coding" )


 INFO_NAMES <-
c( "version","snps","individuals","addr0","addr1","align0","align1","sumgeno","sumgenoE9","method","alignment","isSNPxInd","bitspercode","bytesperblock","codesperblock","header","DoubledIndividuals","leadingcolumns","memInUnits0","meminUnits1","AlignedUnits0","AlignedUnits1","unitsperindiv","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused","unused" )

