
# Local Test Suite for interpret() family functions
# Run this script to verify the logic of:
# - interpret()
# - process_enrichment_input()
# - interpret_hierarchical()
# - interpret_agent()

message("Loading clusterProfiler and dependencies...")
tryCatch({
    devtools::load_all(".")
}, error = function(e) {
    message("Error loading package: ", e$message)
})

# --- Helper Functions for Mocking ---

# Smart Mock for call_llm_fanyi
mock_llm_smart <- function(prompt, model, api_key) {
    # Check prompt content to decide which agent/task is calling
    
    if (grepl("You are 'Agent Cleaner'", prompt)) {
        message("  [Mock] LLM called: Agent Cleaner")
        return(list(
            kept_pathways = c("Pathway 1", "Pathway 2"),
            discarded_pathways = c("Ribosome"),
            reasoning = "Filtered housekeeping."
        ))
    } else if (grepl("You are 'Agent Detective'", prompt)) {
        message("  [Mock] LLM called: Agent Detective")
        return(list(
            key_drivers = c("TF_A", "TF_B"),
            functional_modules = c("Module X"),
            refined_network = list(
                list(source="TF_A", target="Gene1", interaction="activation", reason="db"),
                list(source="TF_B", target="Gene2", interaction="inhibition", reason="lit")
            ),
            network_evidence = "Strong PPI support."
        ))
    } else if (grepl("You are 'Agent Storyteller'", prompt)) {
        message("  [Mock] LLM called: Agent Storyteller")
        return(list(
            overview = "Agent-based interpretation overview.",
            key_mechanisms = "Mechanisms explained.",
            hypothesis = "Hypothesis generated.",
            narrative = "Narrative text."
        ))
    } else if (grepl("characterize the specific biological phenotype", prompt)) {
        message("  [Mock] LLM called: Phenotype Task")
        return(list(
            phenotype = "Mock Phenotype",
            confidence = "Medium",
            reasoning = "Phenotype reasoning.",
            key_processes = c("Process A", "Process B")
        ))
    } else {
        # Default annotation/interpretation task
        message("  [Mock] LLM called: Standard Annotation/Interpretation")
        return(list(
            cell_type = "Mock Cell Type",
            confidence = "High",
            reasoning = "Standard reasoning.",
            markers = c("GeneA", "GeneB"),
            regulatory_drivers = c("TF1", "TF2"),
            overview = "Standard overview.",
            key_mechanisms = list("Mec1" = "Exp1")
        ))
    }
}

mock_llm_empty <- function(prompt, model, api_key) {
    message("  [Mock] LLM returned empty/null.")
    return(NULL)
}

# Inject mock function
inject_mock <- function(mock_func) {
    tryCatch({
        assignInNamespace("call_llm_fanyi", mock_func, ns = "clusterProfiler")
        message("  [Setup] Mock LLM function injected.")
    }, error = function(e) {
        message("  [Setup] Failed to inject mock (namespace locked?): ", e$message)
    })
}

# --- Test Data Creation ---

create_mock_data <- function() {
    # Create a simple data frame acting as enrichment result
    df <- data.frame(
        Cluster = c(rep("Cluster1", 2), rep("Cluster2", 2)),
        ID = c("P1", "P2", "P3", "P4"),
        Description = c("Pathway 1", "Pathway 2", "Pathway 3", "Pathway 4"),
        GeneRatio = c("0.1", "0.1", "0.1", "0.1"),
        p.adjust = c(0.01, 0.02, 0.01, 0.02),
        geneID = c("G1/G2", "G3/G4", "G5/G6", "G7/G8"),
        stringsAsFactors = FALSE
    )
    return(df)
}

# --- Test Cases ---

run_tests <- function() {
    message("\n=== Starting Tests ===\n")
    
    mock_data <- create_mock_data()
    
    # --- Test 1: process_enrichment_input Structure ---
    message("--- Test 1: process_enrichment_input (Unit Test) ---")
    res_input <- process_enrichment_input(mock_data, n_pathways = 5)
    
    if (is.list(res_input) && !is.null(names(res_input))) {
        if (all(names(res_input) %in% c("Cluster1", "Cluster2"))) {
             message("PASS: Output list has correct names.")
        } else {
             message("FAIL: Output list names incorrect. Got: ", paste(names(res_input), collapse=", "))
        }
        
        if (!is.null(res_input$Cluster1$df)) {
             message("PASS: Cluster1 has 'df' component.")
        } else {
             message("FAIL: Cluster1 structure missing 'df'.")
        }
    } else {
        message("FAIL: process_enrichment_input did not return a named list.")
    }

    # --- Test 2: interpret() Standard ---
    message("\n--- Test 2: interpret() Standard Annotation ---")
    inject_mock(mock_llm_smart)
    
    res <- tryCatch({
        interpret(mock_data, task = "annotation")
    }, error = function(e) {
        message("Error running interpret: ", e$message)
        return(NULL)
    })
    
    if (!is.null(res) && inherits(res, "interpretation_list")) {
        message("PASS: Result is interpretation_list")
        if (res[[1]]$cell_type == "Mock Cell Type") {
            message("PASS: Content matches mock")
        } else {
            message("FAIL: Content mismatch")
        }
    } else {
        message("FAIL: Result structure incorrect")
    }

    # --- Test 3: interpret() Phenotype Task ---
    message("\n--- Test 3: interpret() Phenotype Task ---")
    res_pheno <- interpret(mock_data, task = "phenotype")
    
    if (!is.null(res_pheno) && res_pheno[[1]]$phenotype == "Mock Phenotype") {
        message("PASS: Phenotype task handled correctly.")
    } else {
        message("FAIL: Phenotype task failed.")
    }

    # --- Test 4: interpret_hierarchical() ---
    message("\n--- Test 4: interpret_hierarchical() ---")
    # Mock data for hierarchical
    major_data <- data.frame(
        ID = c("MajorP1"), Description = c("Major Pathway"), 
        p.adjust = c(0.01), geneID = c("G1"), 
        Cluster = "MajorClusterA", stringsAsFactors=FALSE
    )
    minor_data <- mock_data # Cluster1, Cluster2
    
    mapping <- list("Cluster1" = "MajorClusterA")
    
    res_hier <- interpret_hierarchical(minor_data, major_data, mapping, task = "cell_type")
    
    if (!is.null(res_hier) && inherits(res_hier, "interpretation_list")) {
        message("PASS: interpret_hierarchical returned list.")
        # Check if Cluster1 was processed (it should be)
        if ("Cluster1" %in% names(res_hier)) {
             message("PASS: Sub-cluster preserved.")
        }
    } else {
        message("FAIL: interpret_hierarchical failed.")
    }

    # --- Test 5: interpret_agent() ---
    message("\n--- Test 5: interpret_agent() ---")
    # interpret_agent returns a list of results (one per cluster)
    res_agent <- interpret_agent(mock_data)
    
    if (!is.null(res_agent) && inherits(res_agent, "interpretation_list")) {
        message("PASS: interpret_agent returned interpretation_list.")
        
        # Check content of Cluster1
        c1 <- res_agent$Cluster1
        if (!is.null(c1$overview) && c1$overview == "Agent-based interpretation overview.") {
            message("PASS: Agent Synthesizer output found.")
        } else {
            message("FAIL: Agent output missing or incorrect.")
        }
        
        if (!is.null(c1$regulatory_drivers) && "TF_A" %in% c1$regulatory_drivers) {
            message("PASS: Agent Detective output found.")
        }
    } else {
        message("FAIL: interpret_agent returned NULL or wrong type.")
        str(res_agent)
    }

    message("\n=== Tests Completed ===")
}

# Run them
run_tests()
