#  #Bead level information
setClass("blum",
         representation=representation( #cannot be ExpressionSet because exprs is not a matrix
           ## Contains information about samples           
           phenoData="AnnotatedDataFrame",
           ## contains information about analytes           
           featureData="AnnotatedDataFrame",
           ## list of bead level data
           ## Stored as samples -> analytes
           exprs="data.table")
         )

setClass("slum", 
  contains="ExpressionSet", 
  representation(unit="character", formula="formula", inv="function", fit="data.frame"))

setClassUnion("blumORslum", members=c("blum", "slum"))
