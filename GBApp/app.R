##GBApp
library("rjson")
library(dplyr)
require('rstudioapi') 
require('nloptr')
require('Matrix')
require('MASS')
require('lpSolve')
library(shiny)
library(readODS)
library(shinyMatrix)
library(reactable)
library(shinydisconnect)
library(htmltools)
library(shinyjs)
library(shinybusy)
library(apexcharter)
library(reshape2)
library(tidyr)
library(ggplot2)
library(shinyalert)
library(rintrojs)
library(later)

directory="/app/"
setwd(directory) 
# UI section of the GBApp
ui <- c(
  
  htmlTemplate("www/index.html",table1=reactableOutput("data_table"),table2=reactableOutput("data_kcat"),
               table11=reactableOutput("data_table1"),table21=reactableOutput("data_kcat1"))
  
  
  
  
  
)


# Server section of the GBApp

server <- function(input, output, session) {

  
 #Give me a tour section 
  tour_steps <- reactive({
    list(
      list(element = "#downloadBtn", intro = "First things first! Download the GBA template and create your model in the provided ODS file before uploading it."),
      list(element = "#Modelpreview", intro = "Here is an example model already uploaded for your reference. You can run and check out this model first, to understand how GBA works. Once you're familiar, you can prepare and upload your own model. After uploading, you will see a preview of your model here and you can modify it."),
      list(element = "#fileuploadsection", intro = "Now, you can upload your own model based on the GBA template you downloaded earlier. Click 'Choose File' to select your ODS file."),
      list(element = "#createmodel", intro = "Alternatively, you can create your model directly within GBApp in this section. Here, you can build, analyze, and download your created model without needing an external file."),
      list(element = "#pathwayviz", intro = "Here, you can see the pathway visualization of your model. This visualization updates dynamically whenever you modify and re-upload your model."),
      list(element = "#loading-container", intro = "Now you can check the validity of your model by clicking 'Check Model'. This ensures that all essential components of the model are correctly defined. If everything is fine after validation, click 'Run' to start the numerical optimization."),
      list(element = "#interactiveplotsintro", intro = "After running the optimization, the corresponding plots will appear here. You can visualize the results, select different axes to explore various aspects of your model, and download the results for further analysis."),
      list(element = "#tutorialinfo", intro = "For more guidance, visit the tutorial section of the app. Here, you'll find detailed explanations and additional resources to help you get the most out of GBApp.")
             )
    
  })
  # Start the tour when button is clicked
  observeEvent(input$start_tour, {
    introjs(session, options = list(steps = tour_steps()))
  })


  #Refreshing the app by clicking on the reset button
  observeEvent(input$refresh, {
    session$reload()
  })
  observeEvent(input$refresh1, {
    session$reload()
  })
  
  #Downloading the Template file, it's in the WWW directory
  output$downloadBtn <- downloadHandler(
    filename = function() {
      # Specify the name of the file to be downloaded
      "GBAmodel.ods"
    },
    content = function(file) {
      # Copy the file from the www directory to the temporary directory for download
      file.copy("www/GBAmodel.ods", file)
    }
  )
  
  ###The creating model section, where different matrices for M, Km,KI,KA,kcat and condition is defined and updated
  observeEvent(input$matrix1 ,{
    matrix_1_col=colnames(as.matrix(input$matrix1))
    matrix_1_col[length(matrix_1_col)]="Ribosome"
    matrix_1_row=rownames(as.matrix(input$matrix1))
    matrix_1_row[length(matrix_1_row)]="Protein"
    matrix5_data=as.vector(input$matrix5)
    updateMatrixInput(session, "matrix1", value = matrix(input$matrix1, as.numeric(input$nrow), as.numeric(input$ncol),dimnames = list(matrix_1_row,matrix_1_col) ))
    updateMatrixInput(session, "matrix2", value = matrix(input$matrix2, as.numeric(input$nrow), as.numeric(input$ncol),dimnames = list(matrix_1_row,matrix_1_col) ))
    updateMatrixInput(session, "matrix3", value = matrix(input$matrix3, as.numeric(input$nrow), as.numeric(input$ncol),dimnames = list(matrix_1_row,matrix_1_col) ))
    updateMatrixInput(session, "matrix4", value = matrix(input$matrix4, as.numeric(input$nrow), as.numeric(input$ncol),dimnames = list(matrix_1_row,matrix_1_col) ))
    updateMatrixInput(session, "matrix5", value = matrix(input$matrix5, 2, as.numeric(input$ncol),dimnames = list(c("kcat_f","kcat_b"),matrix_1_col) ))
    
  })
  
  
  # Update matrix input widgets based on numeric input values
  observeEvent(as.numeric(input$nrow), {
    
    updateMatrixInput(session, "matrix1", value = matrix(0, as.numeric(input$nrow), as.numeric(input$ncol)))
    updateMatrixInput(session, "matrix2", value = matrix(0, as.numeric(input$nrow), as.numeric(input$ncol)))
    updateMatrixInput(session, "matrix3", value = matrix(0, as.numeric(input$nrow), as.numeric(input$ncol)))
    updateMatrixInput(session, "matrix4", value = matrix(0, as.numeric(input$nrow), as.numeric(input$ncol)))
    updateMatrixInput(session, "matrix5", value = matrix(0, 2, as.numeric(input$ncol)))
  })
  
  observeEvent(as.numeric(input$ncol), {
    updateMatrixInput(session, "matrix1", value = matrix(0, as.numeric(input$nrow), as.numeric(input$ncol)))
    updateMatrixInput(session, "matrix2", value = matrix(0, as.numeric(input$nrow), as.numeric(input$ncol)))
    updateMatrixInput(session, "matrix3", value = matrix(0, as.numeric(input$nrow), as.numeric(input$ncol)))
    updateMatrixInput(session, "matrix4", value = matrix(0, as.numeric(input$nrow), as.numeric(input$ncol)))
    updateMatrixInput(session, "matrix5", value = matrix(0, 2, as.numeric(input$ncol)))
  })
  
  observe({
    mat1 <- as.matrix(input$matrix1)
    if (is.null(mat1)) return()  
    
    reactant <- rownames(mat1)
    x_reactants <- reactant[grep("^x_", reactant)]  
    
    rownames_6 <- c("Rho", x_reactants)
    
    ncol_new <- as.numeric(input$ncol6)
    if (is.na(ncol_new) || ncol_new < 1) ncol_new <- 1
    
    old_matrix6 <- as.matrix(input$matrix6)
    nrow_new <- length(rownames_6)
    
    if (is.null(old_matrix6)) {
      new_matrix6 <- matrix(
        0,
        nrow = nrow_new,
        ncol = ncol_new,
        dimnames = list(rownames_6, seq_len(ncol_new))
      )
      updateMatrixInput(session, "matrix6", value = new_matrix6)
      return()
    }
    
    nrow_old <- nrow(old_matrix6)
    ncol_old <- ncol(old_matrix6)
    
    new_matrix6 <- matrix(
      0,
      nrow = nrow_new,
      ncol = ncol_new,
      dimnames = list(rownames_6, seq_len(ncol_new))
    )
    
    min_r <- min(nrow_old, nrow_new)
    min_c <- min(ncol_old, ncol_new)
    
    if (min_r > 0 && min_c > 0) {
      new_matrix6[1:min_r, 1:min_c] <-
        old_matrix6[1:min_r, 1:min_c, drop = FALSE]
    }
    
 
    updateMatrixInput(session, "matrix6", value = new_matrix6)
  })
  
  ##Km data table for both sections
  load(file = "/app/www/Kkm.Rdata")
  
  output$data_table=renderReactable(reactable(Kkm,filterable = TRUE,
                                              resizable = TRUE,class = "table-dark",columns = list(
                                                fieldInfo=colDef(minWidth = 200),
                                                EC1=colDef(minWidth = 50),
                                                EC2=colDef(minWidth = 50),
                                                EC3=colDef(minWidth = 50),
                                                EC4=colDef(minWidth = 50),
                                                organism=colDef(minWidth = 120),
                                                Km=colDef(minWidth = 70,align = "center"),
                                                `molecular weight`=colDef(minWidth = 100),
                                                unit=colDef(minWidth = 50)),fullWidth = TRUE))
  output$data_table1=renderReactable(reactable(Kkm,filterable = TRUE,
                                              resizable = TRUE,class = "table-dark",columns = list(
                                                fieldInfo=colDef(minWidth = 200),
                                                EC1=colDef(minWidth = 50),
                                                EC2=colDef(minWidth = 50),
                                                EC3=colDef(minWidth = 50),
                                                EC4=colDef(minWidth = 50),
                                                organism=colDef(minWidth = 120),
                                                Km=colDef(minWidth = 70,align = "center"),
                                                `molecular weight`=colDef(minWidth = 100),
                                                unit=colDef(minWidth = 50)),fullWidth = TRUE))
  
  ##kcat data tables
  load(file = "/app/www/Kkcat.Rdata")
  
  output$data_kcat=renderReactable(reactable(Kkcat,filterable = TRUE,
                                             resizable = TRUE,class = "table-dark",columns = list(
                                               Substrate=colDef(minWidth = 200),
                                               EC1=colDef(minWidth = 50),
                                               EC2=colDef(minWidth = 50),
                                               EC3=colDef(minWidth = 50),
                                               EC4=colDef(minWidth = 50),
                                               organism=colDef(minWidth = 120),
                                               kcat=colDef(minWidth = 70,align = "center"),
                                               `molecular weight`=colDef(minWidth = 100),
                                               unit=colDef(minWidth = 50)),fullWidth = TRUE))
  output$data_kcat1=renderReactable(reactable(Kkcat,filterable = TRUE,
                                             resizable = TRUE,class = "table-dark",columns = list(
                                               Substrate=colDef(minWidth = 200),
                                               EC1=colDef(minWidth = 50),
                                               EC2=colDef(minWidth = 50),
                                               EC3=colDef(minWidth = 50),
                                               EC4=colDef(minWidth = 50),
                                               organism=colDef(minWidth = 120),
                                               kcat=colDef(minWidth = 70,align = "center"),
                                               `molecular weight`=colDef(minWidth = 100),
                                               unit=colDef(minWidth = 50)),fullWidth = TRUE))
  
  
  #Uploading the model section, if nothing is uploaded, use template
  output$example_model_placeholder <- renderUI({
    if (is.null(input$file)) {
      tagList(
        tags$span("(Example Model)")
        
      )
    } else {
      tags$div()  # blank if a file is uploaded
    }
  })
  source("/app/readmodelods.R")
  data <- reactive({
    
     defaultFile <- file.path("www", "GBAmodel.ods")
  #defaultFile <- "www/GBAmodel.ods"
    # Check input type
    if (!is.null(input$file)) {
      # Read ods file
      
      req(input$file)
      
      tryCatch({
      source("/app/readmodelods.R")
      processFile(input$file$datapath)
      #updateMatrixInput(session, "matrixMtotal", value = Mtotal)
      
      }, error = function(e) {
        # If there's an error, show it in a modal
        showModal(shinyalert(
          title = "Error!",
          text = paste("An error has occurred:", e$message,
                       "\nPlease check that all entries are numeric and contain no missing (NA) values.",
                       "\nPlease check that your model has at least one internal metabolites.",
                       "\nPlease check that all matrices have the same size."),
          type = "error",
          showConfirmButton = TRUE,
          confirmButtonText = "Reset",
          callbackR = function(x) {
            # Reload the Shiny session
            if(!is.null(x)) {
              session$reload()
            }
          }
        ))
      })
    } 
    
    else {
      processFile(defaultFile)
      # Return list of matrix input values
      #source("~/Downloads/GBA R/shiny app/manual_input.R")
      #processmanual(input$matrix1,input$matrix2, input$matrix3, input$matrix4,input$matrix5,input$matrix6)
    }
  })
  ##entries in the created model
  datamanual <- reactive({
  
    source("/app/manual_input.R")
    processmanual(input$matrix1,input$matrix2, input$matrix3, input$matrix4,input$matrix5,input$matrix6)
    })
  
  #update entries of the uploaded model in the preview section by users
  observe( {
    req(data())
    updateMatrixInput(session, "matrixMtotal", value = as.matrix(data()$Mtotal))
    updateMatrixInput(session, "matrixK", value = data()$K)
    updateMatrixInput(session, "matrixKI", value = data()$KI)
    updateMatrixInput(session, "matrixKA", value = data()$KA)
    updateMatrixInput(session, "matrixKcatf", value = data()$kcat_matrix)
    updateMatrixInput(session, "matrixcond", value = data()$conditiontab)

  })

  data2 <- reactive({
    source("/app/data2react.R")
    processmanual2(input$matrixMtotal,input$matrixK, input$matrixKI, input$matrixKA,input$matrixKcatf,
                   input$matrixcond)
    
  })
  

  ##updated model download
  output$downloadModel <- downloadHandler(
    filename = function() {
      paste("updated_data", Sys.Date(), ".ods", sep = "")
    },
    content = function(file) {
      updated_data <- list(
        M= as.data.frame(input$matrixMtotal),
        K=as.data.frame(input$matrixK),
        KI=as.data.frame(input$matrixKI),
        KA=as.data.frame(input$matrixKA),
        kcat =as.data.frame(input$matrixKcatf),
        conditions = as.data.frame(input$matrixcond)
      )
      write_ods(updated_data, path = file,row_names = T,col_names = T)
    }
  )
  
  ###created model download
  output$downloadModel1 <- downloadHandler(
    filename = function() {
      paste("updated_data", Sys.Date(), ".ods", sep = "")
    },
    content = function(file) {
      updated_data <- list(
        M= as.data.frame(input$matrix1),
        K=as.data.frame(input$matrix2),
        KI=as.data.frame(input$matrix3),
        KA=as.data.frame(input$matrix4),
        kcat =as.data.frame(input$matrix5),
        conditions = as.data.frame(input$matrix6)
      )
      write_ods(updated_data, path = file,row_names = T,col_names = T)
    }
  )
  

 
  ##Visualize the pathway with 1 second delay, as the entries should be updated first in the matrices
    later(function() {
      observe( {
   req(data2())
   
   
   
   
   Mtotal <- data2()$Mtotal
   K <- data2()$K
   KI <- data2()$KI
   KA <- data2()$KA
   kcatf <- data2()$kcatf
   kcatb <- data2()$kcatb
   condition <- data2()$condition
   rho_cond <- data2()$rho_cond
   x_cond <- data2()$x_cond
   # index for external reactants
   reaction <- colnames(data2()$Mtotal)
   reactant <- rownames(data2()$Mtotal)
   rho <- rho_cond
  source("/app/map_visualization_pre2.R")
  # visualization
  lb <- rep(0, length(kcatb))
  ub <- rep(0, length(kcatf))
  
  # Update ub and lb based on kcatf and kcatb values
  ub[kcatf > 0] <- rho[1]
  lb[kcatb > 0] <- -1000
  #print(lb)
  ###fluxes
  #metabolite_fluxes=opt_state[, reactant][1,]
  #reaction_fluxes=opt_state[grep("pf_", colnames(opt_state))][1,]
  # Convert to COBRA model
  cobra_model <- convert_stoichiometric_matrix_to_cobra2(Mtotal, reactant, reaction, lb, ub)
  
  json_data <- toJSON(cobra_model)
  #print(json_data)
  # Send the JSON data to the UI
  session$sendCustomMessage(type = 'jsondata', message = json_data)})}, delay = 1)
 
 
##Check model section for Uploaded model
    observeEvent(input$alert, {
      if (input$alert) {
        req(data2())
        print(data2())
        tryCatch({
          mat1 <- as.matrix(data2()$Mtotal) 
          if (is.null(mat1)) return()
          # Check if all values are between -1 and 1
          if (any(mat1 < -1 | mat1 > 1, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Values!",
              text = "All values in the Mass Fraction (M) matrix must be between -1 and 1. Please correct them.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()  # Stop execution if values are invalid
          }
          reactant <- rownames(mat1)
          x_reactants <- reactant[grep("^x_", reactant)]  # Extract rows starting with "x_"
          
          if (length(x_reactants) == 0) {
            shinyalert(
              title = "Missing External Reactants!",
              text = "No external reactants found in the Mass Fraction (M) matrix. Please add at least one external reactant (row starting with 'x_') or fix the model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()  # Stop execution if no external reactants are found
          }
          # 1. Check that matrices M, K, KI, and KA have the same row and column names
          if (!identical(rownames(data2()$Mtotal), rownames(data2()$K)) ||
              !identical(colnames(data2()$Mtotal), colnames(data2()$K)) ||
              !identical(rownames(data2()$Mtotal), rownames(data2()$KI)) ||
              !identical(colnames(data2()$Mtotal), colnames(data2()$KI)) ||
              !identical(rownames(data2()$Mtotal), rownames(data2()$KA)) ||
              !identical(colnames(data2()$K), colnames(data2()$KI)) ||
              !identical(rownames(data2()$K), rownames(data2()$KI)) ||
              !identical(colnames(data2()$K), colnames(data2()$KA)) ||
              !identical(rownames(data2()$K), rownames(data2()$KA)) ||
              !identical(colnames(data2()$KI), colnames(data2()$KA)) ||
              !identical(rownames(data2()$KI), rownames(data2()$KA)) ||
              !identical(colnames(data2()$Mtotal), colnames(data2()$KA))) {
            shinyalert(
              title = "Mismatched Matrix Names!",
              text = "All matrices (M, K, KI, KA) must have identical row and column names. Please fix your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          
          last_col <- mat1[, ncol(mat1)]  # Extract last column of M
          if (all(last_col == 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Ribosome Reaction!",
              text = "The last column in the Mass Fraction (M) matrix, representing the Ribosome reaction, cannot be entirely zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          # Check if there is at least one internal metabolite
          internal_metabolites <- setdiff(reactant, c(x_reactants, "Protein"))  # All rows except "x_" and "Protein"
          
          if (length(internal_metabolites) == 0) {
            shinyalert(
              title = "Missing Internal Metabolites!",
              text = "Your Mass Fraction (M) matrix must contain at least one internal metabolite. Ensure there are rows that are NOT external reactants (x_ prefix) or 'Protein'.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          # Check if all values in matrixK (Michaelis Constant Km) are greater than or equal to zero
          if (any(data2()$K < 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Values in Km Matrix!",
              text = "All values in the Michaelis Constant (Km) matrix must be greater than or equal to zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
      
          if (any(data2()$kcatf < 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Values in kcat Matrix!",
              text = "All values in the Turnover number (kcat) matrix must be greater than or equal to zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          if (any(data2()$kcatb < 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Values in kcat Matrix!",
              text = "All values in the Turnover number (kcat) matrix must be greater than or equal to zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          if (is.null(rownames(mat1)) || any(rownames(mat1) == "")) {
            shinyalert(
              title = "Missing Reactant Names!",
              text = "All rows in the Mass Fraction (M) matrix must have a name. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          
          if (is.null(colnames(mat1)) || any(colnames(mat1) == "")) {
            shinyalert(
              title = "Missing Reaction Names!",
              text = "All columns in the Mass Fraction (M) matrix must have a name. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          
          # 3️⃣ Ensure M matrix is not entirely zeros
          if (all(mat1 == 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid M Matrix!",
              text = "The Mass Fraction (M) matrix cannot be entirely zero. Please enter valid values.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          # 4️⃣ Check if each column in M matrix sums to zero
          col_sums <- colSums(mat1, na.rm = TRUE)
          if (!all(abs(col_sums) < 1e-10)) {  # Tolerance for numerical precision
            shinyalert(
              title = "The Model is not Mass Balanced!",
              text = "Each column in the Mass Fraction (M) matrix must sum to zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          #Ensure kcat matrix is not entirely zeros
          if (all(data2()$kcatf == 0, na.rm = TRUE) && all(data2()$kcatb == 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Turnover number (kcat) Matrix!",
              text = "The Turnover number (kcat) matrix cannot be entirely zero. Please enter at least one value.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          if (all(data2()$x_cond == 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid External concentrations!",
              text = "The External concentrations (Condition matrix) matrix cannot be entirely zero. Please enter at least one value.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
      
          validate(
            need(all(data2()$x_cond >= 0), "All entries in External conditions must be greater equal than zero."),
            need(all(data2()$rho_cond > 0), "All entries in Cell density (Rho) must be greater than zero."),
            #need(all(data2()$K[data2()$Mtotal < 0]>0),"There are missing Michaelis Constant (Km) values in the model")
            # need(data2()$x_cond>0, "There must be at least one reaction.")
            #need(all(sapply(data2(), is.numeric)), "All data columns must be numeric."),
            #need(all(!is.na(sapply(data2(), function(x) sum(is.na(x))))), "Data columns must not contain NA values.")
          )
          # If all validations pass, show success alert
          shinyalert(
            "Validation Successful",
            "Your GBA model is valid.",
            type = "success",
            closeOnClickOutside = TRUE,
            closeOnEsc = TRUE
          
       
            
          )
        }, error = function(e) {
          # If there's an error, show it in a modal
          showModal(shinyalert(
            title = "Error!",
            text = paste("An error has occurred:", e$message
            ),
            type = "error"
            
          ))
        })
        tryCatch({
          validate(
            #need(all(data2()$x_cond > 0), "All entries in External conditions must be greater than zero."),
            #need(all(data2()$rho_cond > 0), "All entries in Cell density (Rho) must be greater than zero."),
            need(all(data2()$K[data2()$Mtotal < 0]>0),paste("There is/are missing Michaelis Constant (Km) values in the model.",
                                                            "\nIf you do not provide the missing Km, a low value of 0.1 will be considered."))
            # need(data2()$x_cond>0, "There must be at least one reaction.")
            #need(all(sapply(data2(), is.numeric)), "All data columns must be numeric."),
            #need(all(!is.na(sapply(data2(), function(x) sum(is.na(x))))), "Data columns must not contain NA values.")
          )
        }, error = function(e) {
          # If there's an error, show it in a modal
          showModal(shinyalert(
            title = "Warning!",
            text =  e$message,
            type = "warning"
            
          ))
        })
        # Show a modal when the button is pressed
        #shinyalert("Validation succesful", "Your GBA model is valid", type = "success")
      }})
    
    ##check model for created model section
    observeEvent(input$alert2, {
      if (input$alert2) {
        req(datamanual())
        
        tryCatch({
          mat1 <- as.matrix(datamanual()$Mtotal) 
          if (is.null(mat1)) return()
          # Check if all values are between -1 and 1
          if (any(mat1 < -1 | mat1 > 1, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Values!",
              text = "All values in the Mass Fraction (M) matrix must be between -1 and 1. Please correct them.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()  # Stop execution if values are invalid
          }
          reactant <- rownames(mat1)
          x_reactants <- reactant[grep("^x_", reactant)]  # Extract rows starting with "x_"
          
          if (length(x_reactants) == 0) {
            shinyalert(
              title = "Missing External Reactants!",
              text = "No external reactants found in the Mass Fraction (M) matrix. Please add at least one external reactant (row starting with 'x_') or fix the model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()  # Stop execution if no external reactants are found
          }
          # 1. Check that matrices M, K, KI, and KA have the same row and column names
          if (!identical(rownames(datamanual()$Mtotal), rownames(datamanual()$K)) ||
              !identical(colnames(datamanual()$Mtotal), colnames(datamanual()$K)) ||
              !identical(rownames(datamanual()$Mtotal), rownames(datamanual()$KI)) ||
              !identical(colnames(datamanual()$Mtotal), colnames(datamanual()$KI)) ||
              !identical(rownames(datamanual()$Mtotal), rownames(datamanual()$KA)) ||
              !identical(colnames(datamanual()$K), colnames(datamanual()$KI)) ||
              !identical(rownames(datamanual()$K), rownames(datamanual()$KI)) ||
              !identical(colnames(datamanual()$K), colnames(datamanual()$KA)) ||
              !identical(rownames(datamanual()$K), rownames(datamanual()$KA)) ||
              !identical(colnames(datamanual()$KI), colnames(datamanual()$KA)) ||
              !identical(rownames(datamanual()$KI), rownames(datamanual()$KA)) ||
              !identical(colnames(datamanual()$Mtotal), colnames(datamanual()$KA))) {
            shinyalert(
              title = "Mismatched Matrix Names!",
              text = "All matrices (M, K, KI, KA) must have identical row and column names. Please fix your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          
          last_col <- mat1[, ncol(mat1)]  # Extract last column of M
          if (all(last_col == 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Ribosome Reaction!",
              text = "The last column in the Mass Fraction (M) matrix, representing the Ribosome reaction, cannot be entirely zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          # Check if there is at least one internal metabolite
          internal_metabolites <- setdiff(reactant, c(x_reactants, "Protein"))  # All rows except "x_" and "Protein"
          
          if (length(internal_metabolites) == 0) {
            shinyalert(
              title = "Missing Internal Metabolites!",
              text = "Your Mass Fraction (M) matrix must contain at least one internal metabolite. Ensure there are rows that are NOT external reactants (x_ prefix) or 'Protein'.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          # Check if all values in matrixK (Michaelis Constant Km) are greater than or equal to zero
          if (any(datamanual()$K < 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Values in Km Matrix!",
              text = "All values in the Michaelis Constant (Km) matrix must be greater than or equal to zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          
          if (any(datamanual()$kcatf < 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Values in kcat Matrix!",
              text = "All values in the Turnover number (kcat) matrix must be greater than or equal to zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          if (any(datamanual()$kcatb < 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Values in kcat Matrix!",
              text = "All values in the Turnover number (kcat) matrix must be greater than or equal to zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          if (is.null(rownames(mat1)) || any(rownames(mat1) == "")) {
            shinyalert(
              title = "Missing Reactant Names!",
              text = "All rows in the Mass Fraction (M) matrix must have a name. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          
          if (is.null(colnames(mat1)) || any(colnames(mat1) == "")) {
            shinyalert(
              title = "Missing Reaction Names!",
              text = "All columns in the Mass Fraction (M) matrix must have a name. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          # 3️⃣ Ensure M matrix is not entirely zeros
          if (all(mat1 == 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid M Matrix!",
              text = "The Mass Fraction (M) matrix cannot be entirely zero. Please enter valid values.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          # 4️⃣ Check if each column in M matrix sums to zero
          col_sums <- colSums(mat1, na.rm = TRUE)
          if (!all(abs(col_sums) < 1e-10)) {  # Tolerance for numerical precision
            shinyalert(
              title = "The Model is not Mass Balanced!",
              text = "Each column in the Mass Fraction (M) matrix must sum to zero. Please check your model.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
 
          #Ensure kcat matrix is not entirely zeros
          if (all(datamanual()$kcatf == 0, na.rm = TRUE) && all(datamanual()$kcatb == 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid Turnover number (kcat) Matrix!",
              text = "The Turnover number (kcat) matrix cannot be entirely zero. Please enter at least one value.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          if (all(datamanual()$x_cond == 0, na.rm = TRUE)) {
            shinyalert(
              title = "Invalid External concentrations!",
              text = "The External concentrations (Condition matrix) matrix cannot be entirely zero. Please enter at least one value.",
              type = "error",
              showConfirmButton = TRUE
            )
            return()
          }
          validate(
            need(all(datamanual()$x_cond >= 0), "All entries in External conditions must be greater equal than zero."),
            need(all(datamanual()$rho_cond > 0), "All entries in Cell density (Rho) must be greater than zero."),
            #need(all(data2()$K[data2()$Mtotal < 0]>0),"There are missing Michaelis Constant (Km) values in the model")
            # need(data2()$x_cond>0, "There must be at least one reaction.")
            #need(all(sapply(data2(), is.numeric)), "All data columns must be numeric."),
            #need(all(!is.na(sapply(data2(), function(x) sum(is.na(x))))), "Data columns must not contain NA values.")
          )
          # If all validations pass, show success alert
          shinyalert(
            "Validation Successful",
            "Your GBA model is valid.",
            type = "success",
            closeOnClickOutside = TRUE,
            closeOnEsc = TRUE
            
            
            
          )
        }, error = function(e) {
          # If there's an error, show it in a modal
          showModal(shinyalert(
            title = "Error!",
            text = paste("An error has occurred:", e$message
            ),
            type = "error"
            
          ))
        })
        tryCatch({
          validate(
            #need(all(data2()$x_cond > 0), "All entries in External conditions must be greater than zero."),
            #need(all(data2()$rho_cond > 0), "All entries in Cell density (Rho) must be greater than zero."),
            need(all(datamanual()$K[mat1 < 0]>0),paste("There is/are missing Michaelis Constant (Km) values in the model.",
                                                            "\nIf you do not provide the missing Km, a low value of 0.1 will be considered."))
            # need(data2()$x_cond>0, "There must be at least one reaction.")
            #need(all(sapply(data2(), is.numeric)), "All data columns must be numeric."),
            #need(all(!is.na(sapply(data2(), function(x) sum(is.na(x))))), "Data columns must not contain NA values.")
          )
        }, error = function(e) {
          # If there's an error, show it in a modal
          showModal(shinyalert(
            title = "Warning!",
            text =  e$message,
            type = "warning"
            
          ))
        })
        # Show a modal when the button is pressed
        #shinyalert("Validation succesful", "Your GBA model is valid", type = "success")
      }})
    
    # Observe changes in the Advanced option for solvers
    observeEvent(input$input1, {
      # Define the choices based on the selection of the first input
      choices <- switch(input$input1,
                        "Non-linear Solver" = c("SLSQP","LBFGS", "MMA"))
                       # "Optimization Mode" = c("Steady-state mode"))
      
      # Update the second selectInput with the new choices
      updateSelectInput(session, "input2", choices = choices)
    })
    observeEvent(input$input11, {
      # Define the choices based on the selection of the first input
      choices <- switch(input$input11,
                        "Non-linear Solver" = c("SLSQP","LBFGS", "MMA"))
                        #"Optimization Mode" = c("Steady-state mode"))
      
      # Update the second selectInput with the new choices
      updateSelectInput(session, "input21", choices = choices)
    })
    
    
    ##Running the uploaded model optimization
  observeEvent(input$runModelBtn | input$runModelBtn3, {
    if (input$runModelBtn | input$runModelBtn3 ) {
      
      req(data2())


     

      Mtotal <- data2()$Mtotal
      K <- data2()$K
      KI <- data2()$KI
      KA <- data2()$KA
      kcatf <- data2()$kcatf
      kcatb <- data2()$kcatb
      condition <- data2()$condition
      rho_cond <- data2()$rho_cond
      x_cond <- data2()$x_cond
      # index for external reactants
      reaction <- colnames(data2()$Mtotal)
      reactant <- rownames(data2()$Mtotal)
      n <- 1:dim(data2()$x_cond)[1]
      
      # internal matrix M
      M <- data2()$Mtotal[-n,]
      
      # number of external reactants
      nx <- dim(data2()$x_cond)[1]
      
      # number of growth conditions
      n_conditions <- dim(data2()$x_cond)[2]
      
      # names of internal reactants
      i_reactant <- data2()$reactant[-n]
      
      # number of internal reactants
      p <- dim(M)[1]
      
      # number of reactions
      r <- dim(M)[2]
     # print(r) 

      # the sum of each M column 
      sM <- colSums(M)
      
      # delete numerical artifacts
      sM[abs(sM) < 1e-10] <- 0
      
      # indexes for reactions: s (transport), e (enzymatic), and ribosome r 
      
      e <- c(1:(r-1))[sM[1:(r-1)] == 0]  
      
      s <- c(1:(r-1))[sM[1:(r-1)] != 0] 
      
      # indexes: m (metabolite), a (all proteins)
      
      m <- 1:(p-1)
      
      # number of transport reactions
      ns <- length(s)
     

      
      source("/app/Kinetics_rxn_shiny.R",local = TRUE)
      
      if (input$runModelBtn | input$runModelBtn3 ) {
        if(input$input2=='IPOPT'){
          source("/app/IPOPT.R",local = TRUE)
        }
        else{
          solver    <- input$input2
          source("/app/GBA_solver_shiny.R",local = TRUE) 
          
        }
      }
      
      opt_state <- matrix(rep(0,(nx+2+p+r+r+r+r)*n_conditions),nrow = n_conditions)
      for (cond in 1:n_conditions) {
        
        rho <- rho_cond[cond]
        
        x  <- x_cond[,cond]
        
        f <- f_opt[cond,]
        
        opt_state[cond,] <- c(conv[cond],mu(f),x,ci(f),tau(ci(f)),v(f),prot(f),f)
      }
      
      colnames(opt_state) <- c("convergence","mu",reactant,paste("tau",reaction),
                               paste("v",reaction),paste("p",reaction),paste("f",reaction))
      
      
      opt_state=as.data.frame(opt_state)
      opt_state2=opt_state
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(opt_state2, file)
        }
      )
      output$downloadData1 <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(opt_state2, file)
        }
      )
      nsucess=sum(length(opt_state$convergence[opt_state$convergence %in% c(0, 4)]))
      if(nsucess==0){
        showModal(shinyalert(
          title = "Error!",
          text = paste("An error has occurred:",
                       "Unfortunately, the numerical optimization was not successful, please try another solver!" ),
          type = "error"
        ))
      }else{
        
        
        
        
        
        totcond=length(opt_state$convergence)
        
       
        
        
        #set an error alert here, if no solution found, it will result in an error.
        p_opt <- opt_state[,paste("p",reaction)]
        phi_opt <- p_opt/c_opt[,p]
        colnames(phi_opt) <- gsub("p", "phi", colnames(phi_opt))
        opt_state=cbind(opt_state,phi_opt)
        colnames(opt_state) <- gsub(" ", "_", colnames(opt_state))
 
         opt_state = opt_state[opt_state$convergence %in% c(0, 4),]
        cx1    <- reactant[grep("x_", reactant)]
        procol <- colnames(opt_state)[grep("p_", colnames(opt_state))]
        phicol <- colnames(opt_state)[grep("phi_", colnames(opt_state))]
        fluxcol <- colnames(opt_state)[grep("v_", colnames(opt_state))]
        metcol <- reactant[-grep("x_", reactant)]
        growth="mu"
        pfs="protein_fractions"
        pcs="protein_concentration"
        flx="reaction_fluxes"
        mtcs="metabolite_concentration"
        
        
          updateSelectInput(session, "x_axis", choices = c(cx1,"Growth Rate" ="mu"))
          updateSelectInput(session, "y_axis", choices = list("Growth rate"=growth,
                                                              "Protein concentration"=procol,
                                                              "Single Protein Fraction"=phicol,
                                                              "Fluxes"=fluxcol,
                                                              "Metabolite concentration"=metcol,
                                                              "Protein Fractions"=pfs,
                                                              "Protein Concentrations"=pcs,
                                                              "Flux of reactions"=flx,
                                                              "Metabolite concentration"=mtcs
                                                              
                                                              
          ))
      
      ##render the results in the table
        output$opt_state=renderReactable(reactable(opt_state,filterable = TRUE,
                                                   resizable = TRUE,class = "table-dark",fullWidth = TRUE,defaultPageSize = 5))
       ##Visualize and update the metabolic pathway
         source("/app/map_visualization_pre.R")
        # visualization
        lb <- rep(0, length(kcatb))
        ub <- rep(0, length(kcatf))
        
        # Update ub and lb based on kcatf and kcatb values
        ub[kcatf > 0] <- rho[1]
        lb[kcatb > 0] <- -1000
        #print(lb)
        ###fluxes
        metabolite_fluxes=opt_state[, reactant][1,]
        reaction_fluxes=opt_state[grep("phi_", colnames(opt_state))][1,]
        # Convert to COBRA model
        cobra_model <- convert_stoichiometric_matrix_to_cobra(Mtotal, reactant, reaction, lb, ub, metabolite_fluxes, reaction_fluxes)
        
        json_data <- toJSON(cobra_model)
        # Send the JSON data to the UI
        session$sendCustomMessage(type = 'jsondata', message = json_data)
     
        
      ##Plots for uploaded model
        v <- reactiveValues(plot = NULL)
      
        
        observe( {
          req(input$x_axis)
            #####
        
            
            if(input$y_axis=='protein_fractions'){
              SJ1 <- gather(opt_state, key = "variable", value = "value", -c(mu,cx1))
              print(SJ1%>%filter(variable %in%  phicol))
              v$plot1 <-SJ1%>%filter(variable %in%  phicol)%>% apex(aes_string(input$x_axis,'value', group='variable'),
                                                                    type = "line")%>%ax_stroke(width = 4) %>%
                ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
                ax_labs(
                  title = paste(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis), "vs", input$y_axis),
                  x = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),
                  y = input$y_axis
                ) %>%
                ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                         
                )%>% 
                ax_colors("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB") %>%
                ax_title(style = list( fontSize=  '20px', cssClass="apexcharts-xaxis-label" ),align = "center"
                )%>%
                
                ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                         axisBorder = list(color = "#666666", height= 3),
                         axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
                
                ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                         axisBorder = list(show = TRUE, color = "#666666",width= 3),
                         axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
                ax_nodata(
                  text = "Sorry no data to visualize",
                  fontSize = "30px"
                ) 
              
              
            } 
            else if(input$y_axis=='protein_concentration'){
              SJ1 <- gather(opt_state, key = "variable", value = "value", -c(mu,cx1))
              print(SJ1%>%filter(variable %in%  procol))
              v$plot3 <-SJ1%>%filter(variable %in%  procol)%>% apex(aes_string(input$x_axis,'value', group='variable'),
                                                                    type = "line")%>%ax_stroke(width = 4) %>%
                ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
                ax_labs(
                  title = paste(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis), "vs", input$y_axis),
                  x = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),
                  y = input$y_axis
                ) %>%
                ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                         
                )%>% 
                ax_colors("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB") %>%
                ax_title(style = list( fontSize=  '20px', cssClass="apexcharts-xaxis-label" ),align = "center"
                )%>%
                
                ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                         axisBorder = list(color = "#666666", height= 3),
                         axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
                
                ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                         axisBorder = list(show = TRUE, color = "#666666",width= 3),
                         axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
                ax_nodata(
                  text = "Sorry no data to visualize",
                  fontSize = "30px"
                ) 
              
              
            }else if(input$y_axis=='reaction_fluxes'){
              SJ1 <- gather(opt_state, key = "variable", value = "value", -c(mu,cx1))
              print(SJ1%>%filter(variable %in%  fluxcol))
              v$plot4 <-SJ1%>%filter(variable %in%  fluxcol)%>% apex(aes_string(input$x_axis,'value', group='variable'),
                                                                     type = "line")%>%ax_stroke(width = 4) %>%
                ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
                ax_labs(
                  title = paste(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis), "vs", input$y_axis),
                  x = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),
                  y = input$y_axis
                ) %>%
                ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                         
                )%>% 
                ax_colors("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB") %>%
                ax_title(style = list( fontSize=  '20px', cssClass="apexcharts-xaxis-label" ),align = "center"
                )%>%
                
                ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                         axisBorder = list(color = "#666666", height= 3),
                         axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
                
                ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                         axisBorder = list(show = TRUE, color = "#666666",width= 3),
                         axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
                ax_nodata(
                  text = "Sorry no data to visualize",
                  fontSize = "30px"
                ) 
              
              
            } else if(input$y_axis=='metabolite_concentration'){
              SJ1 <- gather(opt_state, key = "variable", value = "value", -c(mu,cx1))
              print(SJ1%>%filter(variable %in%  metcol))
              v$plot5 <-SJ1%>%filter(variable %in%  metcol)%>% apex(aes_string(input$x_axis,'value', group='variable'),
                                                                    type = "line")%>%ax_stroke(width = 4) %>%
                ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
                ax_labs(
                  title = paste(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis), "vs", input$y_axis),
                  x = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),
                  y = input$y_axis
                ) %>%
                ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                         
                )%>% 
                ax_colors("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB") %>%
                ax_title(style = list( fontSize=  '20px', cssClass="apexcharts-xaxis-label" ),align = "center"
                )%>%
                
                ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                         axisBorder = list(color = "#666666", height= 3),
                         axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
                
                ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                         axisBorder = list(show = TRUE, color = "#666666",width= 3),
                         axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
                ax_nodata(
                  text = "Sorry no data to visualize",
                  fontSize = "30px"
                ) 
              
              
            }
            else{
              
              ##second graph
              #else{
              v$plot2 <-apex(data = as.data.frame(opt_state),
                             mapping = aes_string(x = input$x_axis, y = input$y_axis),
                             type = "line") %>% 
                ax_stroke(width = 4, colors = "#4477AA") %>%
                ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
                ax_labs(
                  title = paste(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis), "vs", ifelse(input$y_axis == "mu", "Growth rate μ (1/h)", input$y_axis)),
                  x = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),
                  y = ifelse(input$y_axis == "mu", "Growth rate μ (1/h)", input$y_axis)
                ) %>%
                ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                         
                ) %>%
                ax_title(style = list( fontSize=  '20px'),align = "center"
                )%>%
                
                ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                         axisBorder = list(color = "#666666", height= 3),
                         axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
                
                ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                         labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                         axisBorder = list(show = TRUE, color = "#666666",width= 3),
                         axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                         decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
                ax_nodata(
                  text = "Sorry no data to visualize",
                  fontSize = "30px"
                )
              # }
              
              
            }
            
            
          })
        
  
       
        #later(function() {
          output$scatter_plot <- renderApexchart({ v$plot1 })
          output$scatter_plot2 <- renderApexchart({ v$plot2 })
          output$scatter_plot3 <- renderApexchart({ v$plot3 })
          output$scatter_plot4 <- renderApexchart({ v$plot4 })
          output$scatter_plot5 <- renderApexchart({ v$plot5 })
#},1)
        
        output$convergence1 <- renderText({ nsucess})
        output$maxgr <- renderText({sprintf("%.4f", max(opt_state$mu)) })
        output$totpr <- renderText({ sprintf("%.2f",mean(opt_state$Protein)) })
        output$Nrcond <- renderText({ totcond })
        
        ##updating the axes of the plots
        observeEvent(c(input$y_axis,input$x_axis), {
          apexchartProxy("scatter_plot") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = input$y_axis,style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis)), "vs", as.character(input$y_axis)))
            ))
          apexchartProxy("scatter_plot4") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = input$y_axis,style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis)), "vs", as.character(input$y_axis)))
            ))
          apexchartProxy("scatter_plot5") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = ifelse(input$y_axis == "mu", "Growth rate μ (1/h)", input$y_axis),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis)), "vs", as.character(input$y_axis)))
            ))
          apexchartProxy("scatter_plot3") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = input$y_axis,style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis)), "vs", as.character(input$y_axis)))
            ))
          apexchartProxy("scatter_plot2") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = ifelse(input$y_axis == "mu", "Growth rate μ (1/h)", input$y_axis),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis == "mu", "Growth rate μ (1/h)", input$x_axis)), "vs", as.character(ifelse(input$y_axis == "mu", "Growth rate μ (1/h)", input$y_axis))))
            ))
        })
        
      }
        
        
      
      }
  })
  ##Running optimization for create model
  observeEvent(input$runModelBtn2 , { 
    if (input$runModelBtn2 ) {
      
      req(datamanual())
      #print(datamanual())
      
      
      
      
      Mtotal <- datamanual()$Mtotal
         K <- datamanual()$K
      if (any(K[which(Mtotal < 0)]  > 0)) {
        
        K <- as.matrix(K)
        
        } else K <- 0.1*(Mtotal<0)
      #print(K)
      KI <- datamanual()$KI
      KA <- datamanual()$KA
      kcatf <- datamanual()$kcatf
      kcatb <- datamanual()$kcatb
      condition <- datamanual()$condition
      rho_cond <- datamanual()$rho_cond
      x_cond <- datamanual()$x_cond
      # index for external reactants
      reaction <- colnames(datamanual()$Mtotal)
      reactant <- rownames(datamanual()$Mtotal)
      n <- 1:dim(datamanual()$x_cond)[1]
      
      # internal matrix M
      M <- datamanual()$Mtotal[-n,]
      
      # number of external reactants
      nx <- dim(datamanual()$x_cond)[1]
      
      # number of growth conditions
      n_conditions <- dim(datamanual()$x_cond)[2]
      
      # names of internal reactants
      i_reactant <- datamanual()$reactant[-n]
      
      # number of internal reactants
      p <- dim(M)[1]
      
      # number of reactions
      r <- dim(M)[2]
    #  print(r) 
      
      # the sum of each M column 
      sM <- colSums(M)
      
      # delete numerical artifacts
      sM[abs(sM) < 1e-10] <- 0
      
      # indexes for reactions: s (transport), e (enzymatic), and ribosome r 
      
      e <- c(1:(r-1))[sM[1:(r-1)] == 0]  
      
      s <- c(1:(r-1))[sM[1:(r-1)] != 0] 
      
      # indexes: m (metabolite), a (all proteins)
      
      m <- 1:(p-1)
      
      # number of transport reactions
      ns <- length(s)
      
      
      
      source("/app/Kinetics_rxn_shiny.R",local = TRUE)
        if(input$input21=='IPOPT'){
          source("/app/IPOPT.R",local = TRUE)
        }
        else{
          solver    <- input$input21
          source("/app/GBA_solver_shiny.R",local = TRUE) 
          
        }
      
      
      opt_state <- matrix(rep(0,(nx+2+p+r+r+r+r)*n_conditions),nrow = n_conditions)
      for (cond in 1:n_conditions) {
        
        rho <- rho_cond[cond]
        
        x  <- x_cond[,cond]
        
        f <- f_opt[cond,]
        
        opt_state[cond,] <- c(conv[cond],mu(f),x,ci(f),tau(ci(f)),v(f),prot(f),f)
      }
      
      colnames(opt_state) <- c("convergence","mu",reactant,paste("tau",reaction),
                               paste("v",reaction),paste("p",reaction),paste("f",reaction))
      
      
      opt_state=as.data.frame(opt_state)
      opt_state2=opt_state
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(opt_state2, file)
        }
      )
      output$downloadData1 <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(opt_state2, file)
        }
      )
      nsucess=sum(length(opt_state$convergence[opt_state$convergence %in% c(0, 4)]))
      if(nsucess==0){
        showModal(shinyalert(
          title = "Error!",
          text = paste("An error has occurred:",
                       "Unfortunately, the numerical optimization was not successful, please try another solver!" ),
          type = "error"
          
        ))
      }else{
        
        
        
        
        
        totcond=length(opt_state$convergence)
        
        
        
        
        #set an error alert here, if no solution found, it will result in an error.
        p_opt <- opt_state[,paste("p",reaction)]
        phi_opt <- p_opt/c_opt[,p]
        colnames(phi_opt) <- gsub("p", "phi", colnames(phi_opt))
        opt_state=cbind(opt_state,phi_opt)
        colnames(opt_state) <- gsub(" ", "_", colnames(opt_state))
        #opt_state=as.data.frame(opt_state)
        #phi_opt=data.frame(phi_opt)
       # print(colnames(phi_opt))
       # print(colnames(opt_state))
       opt_state = opt_state[opt_state$convergence %in% c(0, 4),]
        cx1    <- reactant[grep("x_", reactant)]
        procol <- colnames(opt_state)[grep("p_", colnames(opt_state))]
        phicol <- colnames(opt_state)[grep("phi_", colnames(opt_state))]
        fluxcol <- colnames(opt_state)[grep("v_", colnames(opt_state))]
        metcol <- reactant[-grep("x_", reactant)]
        growth="mu"
        pfs="protein_fractions"
        pcs="protein_concentration"
        flx="reaction_fluxes"
        mtcs="metabolite_concentration"
        
       
          updateSelectInput(session, "x_axis1", choices = c(cx1,"Growth Rate" ="mu"))
          updateSelectInput(session, "y_axis1", choices = list("Growth rate"=growth,
                                                               "Protein concentration"=procol,
                                                               "Single Protein Fraction"=phicol,
                                                               "Fluxes"=fluxcol,
                                                               "Metabolite concentration"=metcol,
                                                               "Protein Fractions"=pfs,
                                                               "Protein Concentrations"=pcs,
                                                               "Flux of reactions"=flx,
                                                               "Metabolite concentration"=mtcs
                                                               
                                                               
          ))
      
        ##download results of the created model
        output$opt_state=renderReactable(reactable(opt_state,filterable = TRUE,
                                                   resizable = TRUE,class = "table-dark",fullWidth = TRUE,defaultPageSize = 5))
      #pathway visualization of the created model
          source("/app/map_visualization_pre.R")
        # visualization
        lb <- rep(0, length(kcatb))
        ub <- rep(0, length(kcatf))
        
        # Update ub and lb based on kcatf and kcatb values
        ub[kcatf > 0] <- rho[1]
        lb[kcatb > 0] <- -1000
       # print(lb)
        ###fluxes
        metabolite_fluxes=opt_state[, reactant][1,]
        reaction_fluxes=opt_state[grep("phi_", colnames(opt_state))][1,]
        # Convert to COBRA model
        cobra_model <- convert_stoichiometric_matrix_to_cobra(Mtotal, reactant, reaction, lb, ub, metabolite_fluxes, reaction_fluxes)
        
        json_data <- toJSON(cobra_model)
       # print(json_data)
        # Send the JSON data to the UI
        
        
          session$sendCustomMessage(type = 'jsondata1', message = json_data)
       
       ##plots section in the created model
           gg <- reactiveValues(plot = NULL)
        
        observe( {
          req(input$x_axis1)
          #####
          
          if(input$y_axis1=='protein_fractions'){
            SJ1 <- gather(opt_state, key = "variable", value = "value", -c(mu,cx1))
            print(SJ1%>%filter(variable %in%  phicol))
            gg$plot1 <-SJ1%>%filter(variable %in%  phicol)%>% apex(aes_string(input$x_axis1,'value', group='variable'),
                                                                   type = "line")%>%ax_stroke(width = 4) %>%
              ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
              ax_labs(
                title = paste(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1), "vs", input$y_axis1),
                x = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),
                y = input$y_axis1
              ) %>%
              ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                       
              )%>% 
              ax_colors("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB") %>%
              ax_title(style = list( fontSize=  '20px', cssClass="apexcharts-xaxis-label" ),align = "center"
              )%>%
              
              ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                       axisBorder = list(color = "#666666", height= 3),
                       axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
              
              ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                       axisBorder = list(show = TRUE, color = "#666666",width= 3),
                       axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
              ax_nodata(
                text = "Sorry no data to visualize",
                fontSize = "30px"
              ) 
            
            
          } 
          else if(input$y_axis1=='protein_concentration'){
            SJ1 <- gather(opt_state, key = "variable", value = "value", -c(mu,cx1))
            print(SJ1%>%filter(variable %in%  procol))
            gg$plot3 <-SJ1%>%filter(variable %in%  procol)%>% apex(aes_string(input$x_axis1,'value', group='variable'),
                                                                   type = "line")%>%ax_stroke(width = 4) %>%
              ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
              ax_labs(
                title = paste(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1), "vs", input$y_axis1),
                x = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),
                y = input$y_axis1
              ) %>%
              ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                       
              )%>% 
              ax_colors("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB") %>%
              ax_title(style = list( fontSize=  '20px', cssClass="apexcharts-xaxis-label" ),align = "center"
              )%>%
              
              ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                       axisBorder = list(color = "#666666", height= 3),
                       axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
              
              ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                       axisBorder = list(show = TRUE, color = "#666666",width= 3),
                       axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
              ax_nodata(
                text = "Sorry no data to visualize",
                fontSize = "30px"
              ) 
            
            
          }else if(input$y_axis1=='reaction_fluxes'){
            SJ1 <- gather(opt_state, key = "variable", value = "value", -c(mu,cx1))
            print(SJ1%>%filter(variable %in%  fluxcol))
            gg$plot4 <-SJ1%>%filter(variable %in%  fluxcol)%>% apex(aes_string(input$x_axis1,'value', group='variable'),
                                                                    type = "line")%>%ax_stroke(width = 4) %>%
              ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
              ax_labs(
                title = paste(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1), "vs", input$y_axis1),
                x = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),
                y = input$y_axis1
              ) %>%
              ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                       
              )%>% 
              ax_colors("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB") %>%
              ax_title(style = list( fontSize=  '20px', cssClass="apexcharts-xaxis-label" ),align = "center"
              )%>%
              
              ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                       axisBorder = list(color = "#666666", height= 3),
                       axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
              
              ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                       axisBorder = list(show = TRUE, color = "#666666",width= 3),
                       axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
              ax_nodata(
                text = "Sorry no data to visualize",
                fontSize = "30px"
              ) 
            
            
          } else if(input$y_axis1=='metabolite_concentration'){
            SJ1 <- gather(opt_state, key = "variable", value = "value", -c(mu,cx1))
            print(SJ1%>%filter(variable %in%  metcol))
            gg$plot5 <-SJ1%>%filter(variable %in%  metcol)%>% apex(aes_string(input$x_axis1,'value', group='variable'),
                                                                   type = "line")%>%ax_stroke(width = 4) %>%
              ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
              ax_labs(
                title = paste(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1), "vs", input$y_axis1),
                x = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),
                y = input$y_axis1
              ) %>%
              ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                       
              )%>% 
              ax_colors("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB") %>%
              ax_title(style = list( fontSize=  '20px', cssClass="apexcharts-xaxis-label" ),align = "center"
              )%>%
              
              ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                       axisBorder = list(color = "#666666", height= 3),
                       axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
              
              ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                       axisBorder = list(show = TRUE, color = "#666666",width= 3),
                       axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
              ax_nodata(
                text = "Sorry no data to visualize",
                fontSize = "30px"
              ) 
            
            
          }
          else{
            
            ##second graph
            #else{
            gg$plot2 <-apex(data = as.data.frame(opt_state),
                            mapping = aes_string(x = input$x_axis1, y = input$y_axis1),
                            type = "line") %>% 
              ax_stroke(width = 4, colors = "#4477AA") %>%
              ax_grid(yaxis = list(lines = list(show = FALSE)),padding = list(right=20)) %>%
              ax_labs(
                title = paste(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1), "vs", ifelse(input$y_axis1 == "mu", "Growth rate μ (1/h)", input$y_axis1)),
                x = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),
                y = ifelse(input$y_axis1 == "mu", "Growth rate μ (1/h)", input$y_axis1)
              ) %>%
              ax_chart(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)
                       
              ) %>%
              ax_title(style = list( fontSize=  '20px'),align = "center"
              )%>%
              
              ax_xaxis(type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                       axisBorder = list(color = "#666666", height= 3),
                       axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetY= -10)) %>%
              
              ax_yaxis(forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                       labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                       axisBorder = list(show = TRUE, color = "#666666",width= 3),
                       axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                       decimalsInFloat=2,title = list(style=list(fontSize = '20px',colors = "black"), offsetX= -10)) %>%
              ax_nodata(
                text = "Sorry no data to visualize",
                fontSize = "30px"
              )
            # }
            
            
          }
          
          
        })
       
        output$scatter_plot1 <- renderApexchart({ gg$plot1 })
        output$scatter_plot21 <- renderApexchart({ gg$plot2 })
        output$scatter_plot31 <- renderApexchart({ gg$plot3 })
        output$scatter_plot41 <- renderApexchart({ gg$plot4 })
        output$scatter_plot51 <- renderApexchart({ gg$plot5 })
        output$convergence11 <- renderText({ nsucess})
        output$maxgr1 <- renderText({ 
          sprintf("%.4f", max(opt_state$mu))})
        output$totpr1 <- renderText({ sprintf("%.2f",mean(opt_state$Protein)) })
        output$Nrcond1 <- renderText({ totcond })
        
        observe({
          apexchartProxy("scatter_plot1") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = input$y_axis1,style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1)), "vs", as.character(input$y_axis1)))
            ))
          apexchartProxy("scatter_plot41") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = input$y_axis1,style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1)), "vs", as.character(input$y_axis1)))
            ))
          apexchartProxy("scatter_plot51") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = input$y_axis1,style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1)), "vs", as.character(input$y_axis1)))
            ))
          apexchartProxy("scatter_plot31") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = input$y_axis1,style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1)), "vs", as.character(input$y_axis1)))
            ))
          apexchartProxy("scatter_plot21") %>%
            ax_proxy_options(list(#
              
              xaxis = list(
                type = "numeric",tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px',fontWeight= 600, colors = "#666666")),
                axisBorder = list(color = "#666666", height= 3),
                axisTicks=list(color="#666666",height=20,offsetX=0,offsetY=-10),
                decimalsInFloat=2,title = list(text = ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetY= -10)),
              yaxis = list(
                forceNiceScale=TRUE,tickAmount = 5,tickPlacement= 'on',min=0,
                labels = list(style = list(fontSize = '16px', colors = "#666666",fontWeight= 600),align="center"),
                axisBorder = list(show = TRUE, color = "#666666",width= 3),
                axisTicks=list(show=TRUE,color="#666666",width= 20,offsetX=10,offsetY=0),
                decimalsInFloat=2,title = list(text = ifelse(input$y_axis1 == "mu", "Growth rate μ (1/h)", input$y_axis1),style=list(fontSize = '20px',fontWeight= 600,colors = "black"), offsetX= -10)
              ),
              chart=list(animations = list(speed=1500,dynamicAnimation=list(speed=1500)),zoom = list(enabled=FALSE)),
              title=list(text = paste(as.character(ifelse(input$x_axis1 == "mu", "Growth rate μ (1/h)", input$x_axis1)), "vs", as.character(ifelse(input$y_axis1 == "mu", "Growth rate μ (1/h)", input$y_axis1))))
            ))
        })
      }
        
        
        }
    
      
    
      
  

   
      
     
      
    
   
      
      
  })

}

# Run the application 
shinyApp(ui = ui, server = server)



