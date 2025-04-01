convert_stoichiometric_matrix_to_cobra2 <- function(S, metabolites, reactions, lb, ub) {
  # Define the list of all possible cofactors
  #cofactors <- c("ATP", "ADP", "AMP", "NAD", "NADH", "NADP", "NADPH", "FAD", "FADH2", "CoA", "SAM", "THF", "GTP", "GDP", "GMP", "UTP", "UDP", "UMP", "CTP", "CDP", "CMP")
  
  # Convert the cofactors to lowercase for case-insensitive comparison
  #cofactors_lower <- tolower(cofactors)
  


  
  # Create the metabolites list
  metabolites_list <- lapply(1:length(metabolites), function(i) {
    list(
      id = metabolites[i],
      name = "",
      compartment = "",
      notes = list(map_info = list(display_name = metabolites[i]))
    )
  })
  
  # Create the reactions list
  reactions_list <- lapply(1:length(reactions), function(i) {
    reaction_id <- reactions[i]
    metabolites_in_reaction <- S[, i]
    metabolites_list <- setNames(as.list(metabolites_in_reaction[metabolites_in_reaction != 0]), 
                                 metabolites[metabolites_in_reaction != 0])
    #metabolites_in_reaction_lower <- tolower(names(metabolites_list))
    
    # Check for cofactors in the reaction metabolites
    #cofactors_in_reaction <- intersect(metabolites_in_reaction_lower, cofactors_lower)
    
    # Create notes list
    notes_list <- list(map_info = list(reversibility = lb[i] < 0, hidden = FALSE))
    ##cofactors
    # if (length(cofactors_in_reaction) > 0) {
    #   notes_list$map_info$cofactors <- list()
    #   for (cofactor in cofactors_in_reaction) {
    #     cofactor_id <- names(metabolites_list)[metabolites_in_reaction_lower == cofactor]
    #     notes_list$map_info$cofactors[[cofactor_id]] <- list()
    #   }
    # }
    
    list(
      id = reaction_id,
      name = "",
      metabolites = metabolites_list,
      lower_bound = lb[i],
      upper_bound = ub[i],
      gene_reaction_rule = "",
      notes = notes_list
    )
  })
  
  # Create the model template
  model <- list(
    metabolites = metabolites_list,
    reactions = reactions_list,
    genes = list(),
    id = "simple_model",
    compartments = list(),
    notes = list(map_info = list())
  )
  
  # Return the model
  return(model)
}

# # Example usage
# # Define the stoichiometric matrix
# S <- matrix(c(-1, 0, 0, 0, -1, -1, -1, 0, -1, -1,
#               1, -1, 0, 0, 1, 0, 0, -1, -1, 0,
#               0, 1, 0, 0, 0, 1, 0, 1, 0, 1,
#               0, 0, 0, 0, 0, 0, 1, 0, 1, 1,
#               0, 0, 1, -1, 0, 0, 0, 0, 0, -1,
#               0, 0, -1, 1, 0, 0, 0, 0, 1, 0), nrow = 6, byrow = TRUE)
# 
# # Define the metabolites and reactions
# metabolites <- c("A", "B", "C", "ATP", "E", "P")
# reactions <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")
# 
# # Define the bounds and objective coefficients
# lb <- c(0, -1000, 0, 0, 0, 0, 0, -1000, 0, -1000)
# ub <- c(1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000)
# #obj_coef <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# 
# # Define the fluxes for metabolites and reactions
# metabolite_fluxes <- c(0, 0, 0, 0, 0, 0)
# reaction_fluxes <- c(0, 0, 0, 0, 0, 0, 0, 2, 5, -100)
# 
# # Convert to COBRA model
# cobra_model <- convert_stoichiometric_matrix_to_cobra(S, metabolites, reactions, lb, ub, metabolite_fluxes, reaction_fluxes)

