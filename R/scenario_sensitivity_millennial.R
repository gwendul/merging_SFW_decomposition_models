scenario_sensitivity_millennial <- function(
    param_matrix,
    y0 = init_millennial_state(),
    plot_results = TRUE,
    habitat_type = "forest"
) {
  
  # ---------------------------------------------------
  # Load required packages
  # ---------------------------------------------------
  library(deSolve)
  library(rootSolve)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # ---------------------------------------------------
  # Load model & helper code
  # ---------------------------------------------------
  source("R/tree_monomolecular.R")
  source("R/make_tree_forcing.R")
  source("R/millennial_model.R")
  source("R/load_config.R")
  source("R/derive_millennial_parms.R")
  
  # --------------------------------------------------
  # Coerce input to data.frame
  # --------------------------------------------------
  param_df <- as.data.frame(param_matrix, stringsAsFactors = FALSE)
  param_names <- colnames(param_df)
  
  # --------------------------------------------------
  # Load baseline parameters
  # --------------------------------------------------
  parms_base <- load_config("tree_monomolecular")
  parms_base <- modifyList(parms_base, yaml::read_yaml("config/millennial.yml"))
  parms_base <- derive_millennial_parms(parms_base)
  
  # --------------------------------------------------
  # Safety checks
  # --------------------------------------------------
  missing <- setdiff(param_names, names(parms_base))
  if (length(missing) > 0) {
    stop(
      "Parameters not found in MILLENNIAL config: ",
      paste(missing, collapse = ", ")
    )
  }
  
  if (any(!sapply(param_df, is.numeric))) {
    stop("All elements of param_matrix must be numeric parameter values.")
  }
  
  # --------------------------------------------------
  # Prepare storage
  # --------------------------------------------------
  state_names <- names(y0)
  n_scen <- nrow(param_df)
  
  results <- vector("list", n_scen)
  
  cat("Running scenario-based parameter sensitivity\n")
  cat("Number of scenarios:", n_scen, "\n")
  cat("Parameters varied:", paste(param_names, collapse = ", "), "\n")
  
  # --------------------------------------------------
  # Scenario loop
  # --------------------------------------------------
  for (i in seq_len(n_scen)) {
    
    parms_test <- parms_base
    
    # Overwrite parameters with scenario values
    for (p in param_names) {
      parms_test[[p]] <- param_df[[p]][i]
    }
    
    # ---- mass-balance / structural constraints ----
    if ("a_root" %in% param_names) {
      
      flo <- 1 - parms_test$a_root
      
      lprop <- parms_test$a_leaf /
        (parms_test$a_leaf + parms_test$a_wood)
      
      parms_test$a_leaf <- flo * lprop
      parms_test$a_wood <- flo * (1 - lprop)
    }
    
    # Re-derive dependent parameters
    parms_test <- derive_millennial_parms(parms_test)
    
    # Build equilibrium forcing
    parms_test$tree_forcing <-
      if (habitat_type == "forest") {
        make_tree_forcing_equilibrium(parms_test)
      } else {
        make_herb_forcing_equilibrium(parms_test)
      }
    
    # Solve equilibrium
    eq <- rootSolve::stode(
      y     = y0,
      func  = millennial_model,
      parms = parms_test
    )
    
    # Store results
    results[[i]] <- cbind(
      scenario = i,
      param_df[i, , drop = FALSE],
      as.list(eq$y[state_names]),
      Total = sum(eq$y)
    )
    
    cat("  completed scenario", i, "of", n_scen, "\n")
  }
  
  results <- bind_rows(results)
  
  # --------------------------------------------------
  # Optional plotting
  # --------------------------------------------------
  if (plot_results) {
    
    results_long <- results %>%
      pivot_longer(
        cols = all_of(state_names),
        names_to = "state",
        values_to = "value"
      )
    
    print(
      ggplot(results_long,
             aes(x = scenario, y = value)) +
        geom_point() +
        geom_line() +
        facet_wrap(~ state, scales = "free_y") +
        labs(
          x = "Scenario",
          y = "Equilibrium pool size (g C m⁻²)",
          title = "Scenario-based parameter sensitivity (equilibrium)"
        ) +
        theme_minimal()
    )
  }
  
  return(results)
}