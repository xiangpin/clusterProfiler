##' Interpret enrichment results using Large Language Models (LLM)
##'
##' This function sends the enrichment results (top significant pathways) along with 
##' an optional experimental context to an LLM (e.g., DeepSeek) to generate 
##' a biological interpretation, hypothesis, and narrative suitable for a paper.
##'
##' @title interpret
##' @param x An enrichment result object (e.g., `enrichResult` or `gseaResult`).
##' @param context A string describing the experimental background (e.g., "scRNA-seq of mouse myocardial infarction at day 3").
##' @param n_pathways Number of top significant pathways to include in the analysis. Default is 20.
##' @param model The LLM model to use. Default is "deepseek-chat". Supported models include "deepseek-chat", "glm-4", "qwen-turbo" etc.
##' @param api_key The API key for the LLM. If NULL, it tries to fetch from `getOption('yulab_translate')` based on the model.
##' @return A character string containing the LLM-generated interpretation.
##' @author Guangchuang Yu
##' @export
interpret <- function(x, context = NULL, n_pathways = 20, model = "deepseek-chat", api_key = NULL, task = "interpretation", prior = NULL, add_ppi = FALSE, gene_fold_change = NULL) {
    if (missing(x)) {
        stop("enrichment result 'x' is required.")
    }
    
    # Process input into a list of data frames (one per cluster/group)
    res_list <- process_enrichment_input(x, n_pathways)
    
    if (length(res_list) == 0) {
        return("No significant pathways found to interpret.")
    }
    
    # Process each item
    results <- lapply(names(res_list), function(name) {
        df <- res_list[[name]]
        if (nrow(df) == 0) return(NULL)
        
        # Format pathways for prompt
        # We typically need ID, Description, GeneRatio/NES, p.adjust, geneID
        cols_to_keep <- intersect(c("ID", "Description", "GeneRatio", "NES", "p.adjust", "pvalue", "geneID"), names(df))
        pathway_text <- paste(
            apply(df[, cols_to_keep, drop=FALSE], 1, function(row) {
                paste(names(row), row, sep=": ", collapse=", ")
            }),
            collapse="\n"
        )
        
        # Determine prior for this cluster
        current_prior <- NULL
        if (!is.null(prior)) {
            if (length(prior) == 1 && is.null(names(prior))) {
                 current_prior <- prior
            } else if (name %in% names(prior)) {
                 current_prior <- prior[[name]]
            }
        }
        
        # Determine PPI/Hub Genes info if requested
        ppi_network_text <- NULL
        if (add_ppi) {
            # Extract all unique genes from the top pathways
            all_genes <- unique(unlist(strsplit(df$geneID, "/")))
            if (length(all_genes) > 0) {
                # will add 'organism' slot in compareClusterResult in next release
                #
                #
                #
                # Attempt to get PPI network
                # Note: getPPI requires an enrichment object or vector of proteins. 
                # Here we pass the gene vector. We need taxID. 
                # If x is enrichResult, we can try to extract organism. 
                # For simplicity, we try 'auto' if x is enrichResult, otherwise user might need to ensure species is detectable or we skip.
                
                # Check if we can get taxID from the original object 'x' if it corresponds to this name
                # This is tricky because 'x' might be a list or compareClusterResult.
                # Simplified approach: Use the genes directly.
                
                # We need a safe wrapper for getPPI to avoid crashing
                ppi_res <- tryCatch({
                    # Try to infer taxID from original object if possible, otherwise default behavior
                    # Since process_enrichment_input destroys the original object structure, 
                    # we rely on 'x' if it is a single object, or we skip taxID auto-detection if it's too complex
                    
                    # Heuristic: if x is enrichResult, use it.
                    input_for_ppi <- all_genes
                    current_taxID <- "auto"
                    
                    if (inherits(x, "enrichResult") && !is.list(x)) {
                         # Single result
                         current_taxID <- tryCatch(getTaxID(x@organism), error=function(e) "auto")
                    }
                    
                    # Limit to top 50 genes to avoid huge networks and slow API calls
                    input_for_ppi <- head(all_genes, 50)
                    
                    # Get PPI as igraph
                    g <- getPPI(input_for_ppi, taxID = current_taxID, output = "igraph", network_type = "functional")
                    
                    if (!is.null(g)) {
                         # Convert to edge list string for LLM
                         # Format: GeneA -- GeneB (Score: 0.9)
                         el <- igraph::as_data_frame(g, what="edges")
                         if (nrow(el) > 0) {
                             # Filter low confidence edges if needed, but getPPI already has some defaults.
                             # Take top 50 interactions by score if available, or just first 50 to save tokens
                             if ("score" %in% names(el)) {
                                 el <- el[order(el$score, decreasing = TRUE), ]
                             }
                             el_subset <- head(el, 50)
                             
                             # Construct text representation
                             edges_text <- apply(el_subset, 1, function(row) {
                                 score_info <- ""
                                 if ("score" %in% names(row)) score_info <- paste0(" (Score: ", row["score"], ")")
                                 paste0(row["from"], " -- ", row["to"], score_info)
                             })
                             paste(edges_text, collapse = "\n")
                         } else {
                             NULL
                         }
                    } else {
                        NULL
                    }
                }, error = function(e) {
                    # warning("Failed to retrieve PPI info: ", e$message)
                    NULL
                })
                
                if (!is.null(ppi_res)) {
                    ppi_network_text <- ppi_res
                }
            }
        }
        
        # Determine Fold Change info if provided
        fc_text <- NULL
        if (!is.null(gene_fold_change)) {
             # gene_fold_change should be a named vector of logFC
             # We filter for genes present in the pathways
             all_genes <- unique(unlist(strsplit(df$geneID, "/")))
             common_genes <- intersect(all_genes, names(gene_fold_change))
             
             if (length(common_genes) > 0) {
                 # Sort by absolute FC to show most regulated genes
                 fc_subset <- gene_fold_change[common_genes]
                 fc_subset <- fc_subset[order(abs(fc_subset), decreasing = TRUE)]
                 
                 # Take top 20
                 top_fc <- head(fc_subset, 20)
                 fc_text <- paste(names(top_fc), round(top_fc, 2), sep = ":", collapse = ", ")
             }
        }

        # Construct Prompt based on task
        if (task == "annotation" || task == "cell_type") {
            prompt <- construct_annotation_prompt(pathway_text, context, name, current_prior, ppi_network_text, fc_text)
        } else if (task == "phenotype" || task == "phenotyping") {
            prompt <- construct_phenotype_prompt(pathway_text, context, name, ppi_network_text, fc_text)
        } else {
            prompt <- construct_interpretation_prompt(pathway_text, context, ppi_network_text, fc_text)
        }
        
        # Call LLM via fanyi
        res <- call_llm_fanyi(prompt, model, api_key)
        
        # If result is a list (JSON parsed), add cluster name
        if (is.list(res)) {
            res$cluster <- name
            
            # Post-process refined_network if present to be an igraph object
            if (!is.null(res$refined_network)) {
                # refined_network is a list of lists/dataframes from JSON
                # Convert to dataframe
                rn_df <- tryCatch({
                    # It might be a list of lists or a dataframe already depending on jsonlite parsing
                    if (is.data.frame(res$refined_network)) {
                        res$refined_network
                    } else {
                        do.call(rbind, lapply(res$refined_network, as.data.frame))
                    }
                }, error = function(e) NULL)
                
                if (!is.null(rn_df) && nrow(rn_df) > 0) {
                     # Create igraph object
                     # Columns expected: source, target, interaction, reason
                     # Map source/target to from/to for igraph
                     colnames(rn_df)[colnames(rn_df) == "source"] <- "from"
                     colnames(rn_df)[colnames(rn_df) == "target"] <- "to"
                     
                     # Ensure we have from and to
                     if ("from" %in% names(rn_df) && "to" %in% names(rn_df)) {
                         res$network <- igraph::graph_from_data_frame(rn_df, directed = FALSE) # Assuming undirected for simplicity or directed if interaction implies
                     }
                }
            }
        }
        return(res)
    })
    
    names(results) <- names(res_list)
    
    # Return structure
    if (length(results) == 1 && names(results)[1] == "Default") {
        return(results[[1]])
    } else {
        class(results) <- c("interpretation_list", "list")
        return(results)
    }
}

#' Interpret enrichment results using a hierarchical strategy (Major -> Minor clusters)
#'
#' @title interpret_hierarchical
#' @param x_minor Enrichment result for sub-clusters (e.g., compareClusterResult or list of enrichResult).
#' @param x_major Enrichment result for major clusters.
#' @param mapping A named vector mapping sub-cluster IDs (names in x_minor) to major cluster IDs (names in x_major).
#' @param model LLM model.
#' @param api_key API key.
#' @param task Task type, default is "cell_type".
#' @return A list of interpretation results.
#' @author Guangchuang Yu
#' @export
interpret_hierarchical <- function(x_minor, x_major, mapping, model = "deepseek-chat", api_key = NULL, task = "cell_type") {
    
    # 1. Interpret Major Clusters
    message("Step 1: Interpreting Major Clusters to establish lineage context...")
    res_major <- interpret(x_major, context = NULL, model = model, api_key = api_key, task = "cell_type")
    
    # 2. Interpret Sub-clusters with Context
    message("Step 2: Interpreting Sub-clusters using hierarchical constraints...")
    
    # Use internal helper to process x_minor into list of dataframes
    res_list_minor <- process_enrichment_input(x_minor, n_pathways = 20)
    
    results <- lapply(names(res_list_minor), function(name) {
        # name is the sub-cluster ID
        
        # Determine Major Context
        specific_context <- NULL
        if (name %in% names(mapping)) {
            major_id <- mapping[[name]]
            
            # Extract major result
            major_info <- NULL
            if (inherits(res_major, "interpretation_list") && major_id %in% names(res_major)) {
                major_info <- res_major[[major_id]]
            } else if (inherits(res_major, "interpretation")) {
                # Handle case where res_major might be a single result (if only 1 major cluster)
                # Check if it matches major_id or is just default
                if (!is.null(res_major$cluster) && res_major$cluster == major_id) {
                     major_info <- res_major
                } else if (is.null(res_major$cluster)) {
                     # Assume it's the only one
                     major_info <- res_major
                }
            }
            
            if (!is.null(major_info) && !is.null(major_info$cell_type)) {
                 major_label <- major_info$cell_type
                 specific_context <- paste0("Hierarchical Constraint: This cluster is a confirmed subcluster of the '", major_label, "' lineage (identified in major cluster analysis). Please focus on distinguishing the specific subtype or state within this lineage.")
            }
        }
        
        if (is.null(specific_context)) {
            warning(paste("No major lineage context found for sub-cluster:", name))
        }
        
        # Call interpret for this single cluster
        # We pass the dataframe directly
        res <- interpret(res_list_minor[[name]], context = specific_context, model = model, api_key = api_key, task = task)
        
        # Ensure cluster name is preserved
        if (is.list(res)) res$cluster <- name
        
        return(res)
    })
    
    names(results) <- names(res_list_minor)
    class(results) <- c("interpretation_list", "list")
    return(results)
}

process_enrichment_input <- function(x, n_pathways) {
    # Helper to convert object to data frame
    get_df <- function(obj) {
        if (inherits(obj, "compareClusterResult") || inherits(obj, "enrichResult") || inherits(obj, "gseaResult")) {
            return(as.data.frame(obj))
        } else if (is.data.frame(obj)) {
            return(obj)
        }
        stop("Unsupported input type. Expected enrichResult, compareClusterResult, gseaResult, or data.frame.")
    }
    
    # Helper to get top N
    get_top_n <- function(df, n) {
        if (nrow(df) == 0) return(df)
        if ("p.adjust" %in% names(df)) {
            df <- df[order(df$p.adjust), ]
        } else if ("pvalue" %in% names(df)) {
            df <- df[order(df$pvalue), ]
        }
        head(df, n)
    }
    
    # Check if input is a list of enrichment objects (Mixed Database Strategy)
    if (is.list(x) && !inherits(x, "enrichResult") && !inherits(x, "gseaResult") && !inherits(x, "compareClusterResult") && !is.data.frame(x)) {
        # Convert all elements to data frames
        dfs <- lapply(x, get_df)
        
        # Check if they look like compareCluster results (have 'Cluster' column)
        has_cluster <- all(sapply(dfs, function(d) "Cluster" %in% names(d)))
        
        combined_df <- do.call(rbind, dfs)
        
        if (has_cluster) {
            # Split by Cluster and get top N for each cluster
            df_list <- split(combined_df, combined_df$Cluster)
            return(lapply(df_list, function(d) get_top_n(d, n_pathways)))
        } else {
            # Assume single group (e.g. list of enrichResult for same sample)
            return(list(Default = get_top_n(combined_df, n_pathways)))
        }
    } else {
        # Single object
        df <- get_df(x)
        if ("Cluster" %in% names(df)) {
            # compareClusterResult
            df_list <- split(df, df$Cluster)
            return(lapply(df_list, function(d) get_top_n(d, n_pathways)))
        } else {
            # enrichResult / gseaResult
            return(list(Default = get_top_n(df, n_pathways)))
        }
    }
}

construct_interpretation_prompt <- function(pathways, context, ppi_network = NULL, fold_change = NULL) {
    base_prompt <- "You are an expert biologist and bioinformatics researcher. I have performed functional enrichment analyses using multiple databases (e.g., KEGG, Reactome, GO, ChEA/Transcription Factors, Disease Ontologies)."
    
    if (!is.null(context) && context != "") {
        base_prompt <- paste0(base_prompt, "\n\nExperimental Context:\n", context)
    }
    
    base_prompt <- paste0(base_prompt, "\n\nTop Enriched Terms (Mixed Sources):\n", pathways)
    
    if (!is.null(ppi_network)) {
        base_prompt <- paste0(base_prompt, "\n\nPPI Network (Edge List from STRING):\n", ppi_network)
    }
    
    if (!is.null(fold_change)) {
        base_prompt <- paste0(base_prompt, "\n\nTop Regulated Genes (Log2 Fold Change):\n", fold_change)
    }
    
    base_prompt <- paste0(base_prompt, "\n\nPlease use a **Chain-of-Thought** approach to analyze these results before generating the final report. Follow these reasoning steps:
1. **Source Deconvolution**: Identify the nature of the enriched terms. Distinguish between:
    - **Biological Processes/Pathways** (e.g., 'Cell Cycle', 'TCR Signaling') -> WHAT is happening.
    - **Upstream Regulators/TFs** (e.g., 'E2F1', 'NFKB1 target genes') -> WHO is driving it.
    - **Phenotypes/Diseases** (e.g., 'Inflammation') -> WHAT is the outcome.
2. **Gene-Level Analysis**: Identify shared key genes and unique functional modules across the pathways.
3. **Causal Integration**: Construct a regulatory narrative connecting the Regulators (TFs) to the Processes (Pathways). For example: 'Enrichment of E2F1 targets explains the observed upregulation of Cell Cycle pathways'.
4. **Contextual Mapping**: Map these themes to the provided experimental context (e.g., tissue, treatment, timepoint) to distinguish expected vs. unexpected findings.
5. **Network Analysis**: Use the provided PPI connections to identify **functional modules** (e.g., a receptor-ligand pair or a protein complex). Identify key hubs that drive the biological process.
6. **Network Refinement**: Prune the provided PPI network to retain only the most biologically relevant edges that explain the context.

Based on this deep analysis, please provide a comprehensive biological interpretation formatted as a JSON object with the following keys:
- overview: A high-level summary of the key biological processes identified.
- regulatory_drivers: Identify key transcription factors or master regulators found in the enrichment list (or inferred from hubs) and explain their potential role in driving the observed pathways.
- key_mechanisms: Explain the underlying biological mechanisms, grouping pathways into major themes.
- crosstalk: Discuss potential interactions and regulatory networks between these pathways.
- hypothesis: Formulate a coherent biological hypothesis connecting the pathways ('What') to the biological meaning ('So What').
- narrative: Write a cohesive paragraph suitable for the 'Results' or 'Discussion' section of a scientific paper.
- refined_network: A list of objects representing the refined core regulatory network. Each object should have 'source', 'target', 'interaction' (e.g., binding, activation, inhibition), and 'reason' (why this edge is relevant).
- network_evidence: Describe specific protein complexes or signaling modules found in the network that support your conclusion.

**GROUNDING INSTRUCTION:**
Every biological claim or interpretation MUST be supported by specific evidence from the provided enrichment results.
- When mentioning a biological process, cite the supporting pathway name(s) in parentheses.
- Example: 'The cluster shows signs of proliferation (supported by: Cell cycle, DNA replication)'.
- Do not make claims that cannot be directly inferred from the provided list.

Please be scientifically rigorous, citing standard biological knowledge where appropriate, and avoid hallucinations.

Ensure the response is a valid JSON object. Do not include any markdown formatting (like ```json).")
    
    return(base_prompt)
}

construct_annotation_prompt <- function(pathways, context, cluster_id, prior = NULL, ppi_network = NULL, fold_change = NULL) {
    base_prompt <- paste0("You are an expert cell biologist. I have a cell cluster (", cluster_id, ") from a single-cell RNA-seq experiment.")
    
    if (!is.null(context) && context != "") {
        base_prompt <- paste0(base_prompt, "\n\nExperimental Context:\n", context)
    }
    
    if (!is.null(prior) && prior != "") {
        base_prompt <- paste0(base_prompt, "\n\nPreliminary Annotation (from automated tool):\n", prior)
        base_prompt <- paste0(base_prompt, "\n\nEnriched Terms (Mixed Sources: Pathways/TFs):\n", pathways)
        
        if (!is.null(ppi_network)) {
            base_prompt <- paste0(base_prompt, "\n\nPPI Network (Edge List from STRING):\n", ppi_network)
        }
        
        if (!is.null(fold_change)) {
            base_prompt <- paste0(base_prompt, "\n\nTop Regulated Genes (Log2 Fold Change):\n", fold_change)
        }
        
        base_prompt <- paste0(base_prompt, "\n\nTask:
Please validate and refine the preliminary annotation based on the enrichment evidence. Use the following logic:
1. **Source Deconvolution**: Identify enriched Transcription Factors (TFs) vs. Biological Pathways.
2. **Validation**: Do the pathways support the preliminary label?
3. **Refinement**: Can you assign a more specific subtype or functional state? (e.g., 'T cells' -> 'Th17' if 'STAT3' and 'IL17 signaling' are enriched).
4. **Correction**: If the evidence strongly contradicts the preliminary label, propose the correct cell type.
5. **Network Analysis**: Use the provided PPI connections to identify **functional modules** (e.g., a receptor-ligand pair or a protein complex).
6. **Integration**: Use these identified modules and TFs to solidify your cell type classification.

**GROUNDING INSTRUCTION:**
- Every conclusion must be backed by specific pathways or markers from the input list.
- In the 'reasoning' field, explicitly mention which terms support your decision.
- Example: 'Refined to Exhausted T cells because of the presence of PD-1 signaling and LAG3 expression.'

Provide the result as a JSON object with the following keys:
- cell_type: The final identified cell type label (refined or corrected).
- refinement_status: 'Confirmed', 'Refined', or 'Corrected'.
- confidence: 'High', 'Medium', or 'Low'.
- reasoning: Explain why you confirmed, refined, or corrected the label, citing specific evidence (Pathways, TFs, Markers).
- regulatory_drivers: Identify key TFs or master regulators that define this cell state.
- markers: A list of key markers or pathways from the input that support this decision.
- refined_network: A list of objects representing the refined core regulatory network. Each object should have 'source', 'target', 'interaction', and 'reason'.
- network_evidence: Describe specific protein complexes or signaling modules found in the network that support your conclusion.

Ensure the response is a valid JSON object. Do not include any markdown formatting (like ```json).")
    } else {
        base_prompt <- paste0(base_prompt, "\n\nEnriched Terms (Mixed Sources: Pathways/TFs):\n", pathways)
        
        if (!is.null(ppi_network)) {
            base_prompt <- paste0(base_prompt, "\n\nPPI Network (Edge List from STRING):\n", ppi_network)
        }
        
        if (!is.null(fold_change)) {
            base_prompt <- paste0(base_prompt, "\n\nTop Regulated Genes (Log2 Fold Change):\n", fold_change)
        }
        
        base_prompt <- paste0(base_prompt, "\n\nBased on these enrichment results, identify the cell type of this cluster. 
Use the following logic:
1. **Source Deconvolution**: Distinguish between Cell Type Markers, Biological Pathways, and Upstream TFs.
2. Analyze the specific marker genes or cell type signatures present in the terms.
3. Analyze the functional pathways (e.g., KEGG, Reactome) to infer cell state or function.
4. Combine these evidences to assign a specific Cell Type Label and, if possible, a functional state (e.g., 'CD8+ T cells - Exhausted').
5. **Network Analysis**: Use the provided PPI connections to identify **functional modules**.
6. **Integration**: Use these identified modules and TFs to solidify your cell type classification.

**GROUNDING INSTRUCTION:**
- Do not guess based on weak evidence.
- Explicitly cite the markers, pathways, or TFs that led to your classification in the 'reasoning' section.

Provide the result as a JSON object with the following keys:
- cell_type: The identified cell type label.
- confidence: 'High', 'Medium', or 'Low'.
- reasoning: A brief explanation of why this cell type was assigned, citing specific markers or pathways.
- regulatory_drivers: Identify key TFs or master regulators that define this cell type/state.
- markers: A list of key markers or pathways from the input that support this decision.
- refined_network: A list of objects representing the refined core regulatory network. Each object should have 'source', 'target', 'interaction', and 'reason'.
- network_evidence: Describe specific protein complexes or signaling modules found in the network that support your conclusion.

Ensure the response is a valid JSON object. Do not include any markdown formatting (like ```json).")
    }
    
    return(base_prompt)
}

construct_phenotype_prompt <- function(pathways, context, group_id, ppi_network = NULL, fold_change = NULL) {
    base_prompt <- paste0("You are an expert biologist. I have a list of enriched pathways/terms for a biological group (", group_id, "). The enrichment may include results from multiple databases (e.g., Pathways, TFs, Ontologies).")
    
    if (!is.null(context) && context != "") {
        base_prompt <- paste0(base_prompt, "\n\nExperimental Context:\n", context)
    }
    
    base_prompt <- paste0(base_prompt, "\n\nEnriched Terms:\n", pathways)
    
    if (!is.null(ppi_network)) {
        base_prompt <- paste0(base_prompt, "\n\nPPI Network (Edge List from STRING):\n", ppi_network)
    }
    
    if (!is.null(fold_change)) {
        base_prompt <- paste0(base_prompt, "\n\nTop Regulated Genes (Log2 Fold Change):\n", fold_change)
    }
    
    base_prompt <- paste0(base_prompt, "\n\nBased on these enrichment results, characterize the specific biological phenotype or functional state of this group.
Use the following logic:
1. **Source Deconvolution**: Separate observed processes (Pathways) from drivers (TFs).
2. Synthesize the enriched terms to identify the dominant biological theme (e.g., Inflammation, Cell Cycle, Metabolism, Stress Response).
3. Be specific about the direction or nature of the state (e.g., 'M1 Macrophage Polarization' is better than just 'Immune response'; 'G2/M Arrest' is better than 'Cell Cycle').
4. Assign a concise Phenotype Label.
5. **Network Analysis**: Use the provided PPI connections to identify **functional modules** and key drivers.

**GROUNDING INSTRUCTION:**
- Every phenotype claim must be supported by evidence.
- In the 'reasoning' section, list the specific pathways or TFs that define this phenotype.

Provide the result as a JSON object with the following keys:
- phenotype: A concise label for the biological phenotype/state.
- confidence: 'High', 'Medium', or 'Low'.
- reasoning: Explanation of how the terms support this phenotype.
- regulatory_drivers: Identify key TFs or master regulators that drive this phenotype.
- key_processes: A list of key pathways/terms that define this phenotype.
- refined_network: A list of objects representing the refined core regulatory network. Each object should have 'source', 'target', 'interaction', and 'reason'.
- network_evidence: Describe specific protein complexes or signaling modules found in the network that support your conclusion.

Ensure the response is a valid JSON object. Do not include any markdown formatting (like ```json).")
    
    return(base_prompt)
}

#' @importFrom jsonlite toJSON fromJSON
call_llm_fanyi <- function(prompt, model, api_key) {
    if (!requireNamespace("fanyi", quietly = TRUE)) {
        stop("Package 'fanyi' is required for interpret(). Please install it.")
    }
    
    # Call fanyi::chat_request
    res_content <- tryCatch({
        fanyi::chat_request(x = prompt, model = model, api_key = api_key, max_tokens = 4096)
    }, error = function(e) {
        stop("Failed to call fanyi::chat_request. Error: ", e$message)
    })
    
    # Try to parse JSON response if the prompt asked for JSON
    tryCatch({
        # Clean up potential markdown code blocks like ```json ... ```
        json_str <- res_content
        if (grepl("```json", json_str)) {
            json_str <- sub(".*?```json\\s*", "", json_str)
            json_str <- sub("\\s*```.*", "", json_str)
        } else if (grepl("```", json_str)) {
            json_str <- sub(".*?```\\s*", "", json_str)
            json_str <- sub("\\s*```.*", "", json_str)
        }
        
        if (!requireNamespace("jsonlite", quietly = TRUE)) {
             stop("Package 'jsonlite' is required.")
        }
        
        parsed_res <- jsonlite::fromJSON(json_str)
        class(parsed_res) <- c("interpretation", class(parsed_res))
        return(parsed_res)
    }, error = function(e) {
        warning("Failed to parse JSON response from LLM. Returning raw text. Error: ", e$message)
        return(res_content)
    })
}

#' @method print interpretation
#' @export
print.interpretation <- function(x, ...) {
    # Check if it is an annotation result
    if (!is.null(x$cell_type)) {
        cat("## Cell Type Annotation\n\n")
        if (!is.null(x$cluster)) cat(sprintf("### Cluster: %s\n\n", x$cluster))
        
        cat(sprintf("**Cell Type:** %s\n", x$cell_type))
        cat(sprintf("**Confidence:** %s\n", x$confidence))
        cat("\n**Reasoning:**\n", x$reasoning, "\n")
        
        if (!is.null(x$markers)) {
            cat("\n**Supporting Markers/Pathways:**\n")
            if (is.list(x$markers) || length(x$markers) > 1) {
                cat(paste("-", unlist(x$markers), collapse="\n"), "\n")
            } else {
                cat(x$markers, "\n")
            }
        }
        cat("\n")
        return(invisible(x))
    }

    # Check if it is a phenotyping result
    if (!is.null(x$phenotype)) {
        cat("## Phenotype Characterization\n\n")
        if (!is.null(x$cluster)) cat(sprintf("### Group/Cluster: %s\n\n", x$cluster))
        
        cat(sprintf("**Phenotype:** %s\n", x$phenotype))
        cat(sprintf("**Confidence:** %s\n", x$confidence))
        cat("\n**Reasoning:**\n", x$reasoning, "\n")
        
        if (!is.null(x$key_processes)) {
            cat("\n**Key Processes:**\n")
            if (is.list(x$key_processes) || length(x$key_processes) > 1) {
                cat(paste("-", unlist(x$key_processes), collapse="\n"), "\n")
            } else {
                cat(x$key_processes, "\n")
            }
        }
        cat("\n")
        return(invisible(x))
    }

    cat("## Interpretation Result\n\n")
    if (!is.null(x$cluster)) cat(sprintf("### Cluster: %s\n\n", x$cluster))
    
    if (!is.null(x$overview)) {
        cat("### 1. Overview\n")
        cat(x$overview, "\n\n")
    }
    
    if (!is.null(x$regulatory_drivers)) {
        cat("### 2. Regulatory Drivers (TFs/Hubs)\n")
        if (is.list(x$regulatory_drivers) || length(x$regulatory_drivers) > 1) {
             # If list or vector
             drivers <- unlist(x$regulatory_drivers)
             cat(paste("-", drivers, collapse="\n"), "\n\n")
        } else {
             cat(x$regulatory_drivers, "\n\n")
        }
    }
    
    if (!is.null(x$key_mechanisms)) {
        cat("### 3. Key Mechanisms\n")
        if (is.list(x$key_mechanisms)) {
            for (mechanism_name in names(x$key_mechanisms)) {
                cat(sprintf("#### %s\n", mechanism_name))
                mechanism <- x$key_mechanisms[[mechanism_name]]
                
                # Check if mechanism is a list (structured) or just a character string
                if (is.list(mechanism)) {
                    if (!is.null(mechanism$explanation)) {
                        cat(mechanism$explanation, "\n")
                    }
                    if (!is.null(mechanism$pathways)) {
                        cat("**Pathways:** ", paste(mechanism$pathways, collapse = ", "), "\n")
                    }
                    if (!is.null(mechanism$genes)) {
                        cat("**Key Genes:** ", paste(head(mechanism$genes, 10), collapse = ", "), ifelse(length(mechanism$genes) > 10, "...", ""), "\n")
                    }
                } else {
                    # mechanism is likely a simple character string description
                    cat(mechanism, "\n")
                }
                cat("\n")
            }
        } else {
            cat(x$key_mechanisms, "\n\n")
        }
    }
    
    if (!is.null(x$crosstalk)) {
        cat("### 4. Crosstalk & Interactions\n")
        cat(x$crosstalk, "\n\n")
    }
    
    if (!is.null(x$hypothesis)) {
        cat("### 5. Hypothesis\n")
        if (is.list(x$hypothesis)) {
            if (!is.null(x$hypothesis$what)) {
                cat("**Observation (What):** ", x$hypothesis$what, "\n\n")
            }
            if (!is.null(x$hypothesis$so_what)) {
                cat("**Implication (So What):** ", x$hypothesis$so_what, "\n\n")
            }
        } else {
            cat(x$hypothesis, "\n\n")
        }
    }
    
    if (!is.null(x$narrative)) {
        cat("### 6. Narrative Draft\n")
        cat(x$narrative, "\n\n")
    }
    
    if (!is.null(x$network)) {
        cat("### 7. Refined Regulatory Network\n")
        # Simple ASCII visualization of the network
        # Edge list with interaction type
        el <- igraph::as_data_frame(x$network, what = "edges")
        if (nrow(el) > 0) {
            cat("Key Interactions:\n")
            for (i in 1:nrow(el)) {
                interaction_type <- ifelse("interaction" %in% names(el), paste0(" (", el[i, "interaction"], ")"), "")
                reason <- ifelse("reason" %in% names(el), paste0(" - ", el[i, "reason"]), "")
                cat(sprintf("  %s -- %s%s%s\n", el[i, "from"], el[i, "to"], interaction_type, reason))
            }
            cat("\n")
        }
        
        if (!is.null(x$network_evidence)) {
             cat("**Network Evidence:**\n")
             cat(x$network_evidence, "\n\n")
        }
    }
    
    invisible(x)
}

#' @method print interpretation_list
#' @export
print.interpretation_list <- function(x, ...) {
    cat("# Enrichment Interpretation / Annotation Report\n\n")
    for (i in seq_along(x)) {
        print(x[[i]])
        cat("---\n\n")
    }
    invisible(x)
}
