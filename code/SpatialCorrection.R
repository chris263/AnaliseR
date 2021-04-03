prepare_data <- function(va_centered, Traits, Treatments, Location, method){
  library(package = "tidyverse")
  library(package = "INLA")

  Dat = va_centered %>%filter(Location == UQ(Location))
  Rep = unique(Dat$Replication)

  rowCol <- matrix(data=0,nrow = nrow(Dat), ncol = 3)
  colnames(rowCol) <- c("colNumber","rowNumber","plot")

  create_rowCol <- function(col_number, row_number, start_col, max_col, counter){
    for (i in 1:nrow(Dat_Rep)){
      if(col_number > max_col){
        col_number = start_col
        row_number = row_number+1
      }
      rowCol[counter,1] <<- col_number
      rowCol[counter,2] <<- row_number
      rowCol[counter,3] <<- counter+100
      col_number = col_number+1
      counter=counter+1
    }
  }

  for (i in 1:length(Rep)){
    Dat_Rep = Dat %>%filter(Replication == (UQ(Rep[i])))
    start_col = (i*2)-1
    max_col = start_col+1
    row_number = 1

    col_number = start_col
    counter = (length(Treatments)*i) - (length(Treatments)-1)
    create_rowCol(col_number, row_number, start_col, max_col, counter)
  }
  centeredRowCol <- cbind(Dat[,1:2],rowCol,Dat[,3:ncol(Dat)])
  centeredRowCol$intercept = 1

  return(centeredRowCol)

}

spatial_correction <- function(inputSPDE, Treatments, Traits, Location, method){
  # ---- MODEL Single trial - RowCol AR1xAR1 ----
  meshCoor <- inputSPDE %>%filter(Location == UQ(Location))
  meshCoor$intercept = 1

  if(method == "RowCol"){
    # Make formula
    for(i in 1:length(Traits)){
      colnames(meshCoor)[6+i] <- "trait"
      if(i > 1){
        colnames(meshCoor)[6+i-1] <- "Oldtrait"
      }
      model_AR1_cent = trait ~ f(plot, model = "ar1") + f(rowNumber, model = "ar1")  + f(rowNumber, model = "ar1")
      Fit_model_AR1_cent  = inla(formula = as.formula(model_AR1_cent),  data = meshCoor)
      result_RowCol <- data.frame(Fit_model_AR1_cent$summary.random$plot$ID, Fit_model_AR1_cent$summary.random$plot$mean)
      # colnames(result_RowCol) <- c("plot", get(Traits[i]))

      # Structuring the result
      # Final matrix with results
      if(i==1){
        final_matrix <- cbind(meshCoor[,1:6], result_RowCol[,2])
      }else{
        final_matrix <- cbind(final_matrix, result_RowCol[,2])
      }
    }
  }else if(method == "SPDE"){
    # Make the mesh
      # All trials have the same distribution. So, we need to create the mesh just once
      Locations = as.matrix(cbind(2*meshCoor$colNumber, meshCoor$rowNumber))
      nLoc = nrow(Locations) # Number of observations
      Mesh = inla.mesh.2d(loc = Locations, max.edge = c(4,6), offset = c(5,5), cutoff = c(1))
      plot(Mesh)
      points(Locations, col = "blue", pch = 20)

      SpdeStat <- inla.spde2.matern(mesh = Mesh, alpha = 2)

      MeshIndex <- inla.spde.make.index(name = "FieldI", n.spde = SpdeStat$n.spde)
      APlots <- inla.spde.make.A(Mesh, loc = Locations, index = c(1:nLoc))

      for(i in 1:length(Traits)){
          Pheno <- Traits[i]
          # Make formula
          ModelSTSPDE  = Pheno ~ intercept-1 + f(plot, model = "ar1") + f(FieldI, model = SpdeStat)
          stack = inla.stack(data = list( Pheno = meshCoor[,i+6]), A = list(APlots,1),
                        effects = list(c(MeshIndex, list(intercept = 1)),
                        c(list(plot = meshCoor$plot))), tag = "Estimate")
          # Run inla
          FitSTSPDE  = inla(formula = as.formula(ModelSTSPDE),  data = inla.stack.data(stack),
                                control.predictor = list(A = inla.stack.A(stack)))

          result_SPDE <- data.frame(FitSTSPDE$summary.random$plot$ID, FitSTSPDE$summary.random$plot$mean)

          # Structuring the result
          # Final matrix with results
          if(i ==1){
              final_matrix <<- cbind(inputSPDE[,1:6], result_SPDE[,2])
          }else{
              final_matrix <<- cbind(final_matrix, result_SPDE[,2])
              }
      }
  }
  return(final_matrix)
}

#
# }
#
#   # ---- MODEL Single trial - SPDE ----
#



#     # Define the spde model and make index for spde object
#
#
#
#     # Store the data, the A matrix and the effects in a stack object
#
#
#
#   }
# }
#
# colnames(final_matrix) <- c("germplasmName","Replication","colNumber","rowNumber","plot","Location",Traits)
#
#
#


#
#
# head(final_matrix)
# library("lattice")
# library("viridisLite")
# # data with no correction
# pre_pheno = data.frame(as.numeric(centeredRowCol$rowNumber),as.numeric(centeredRowCol$colNumber), as.numeric(centeredRowCol$kgHa))
# dat_rows <- max(pre_pheno[,1])
# dat_cols <- max(pre_pheno[,2])
# pheno_matrix = matrix("NA", nrow=dat_rows,ncol=dat_cols)
# dim(pheno_matrix)
# i=1
# for (i in 1:nrow(pre_pheno)) {
# ref_row = pre_pheno[i,1]
# ref_col = pre_pheno[i,2]
# pheno_matrix[ref_row,ref_col] <- pre_pheno[i,3]
# }
# # rownames(pheno_matrix) <- paste("", 1:dat_rows, sep = "")
# levelplot( t(pheno_matrix), xlab="Columns", ylab="Rows",col.regions=heat.colors(100))
#
