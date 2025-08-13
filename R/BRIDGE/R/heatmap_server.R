#' @export
heatmap_server <- function(id, rv) {
  
  moduleServer(id, function(input, output, session) {
    shiny::observe({
      lapply(rv$table_names, function(tbl_name) {
        output[[paste0("raw_ht_", tbl_name)]] <- shiny::renderPlot({
          ComplexHeatmap::Heatmap(rv$ht_matrix[[tbl_name]], show_row_dend = FALSE, show_row_names = FALSE)
        })
      })
    })
  })
}

#' @export
raw_heatmap_server <- function(input, output, session, rv) {
    shiny::observe({
    lapply(rv$table_names, function(tbl_name) {

        raw_heatmap_plot <- ExtendedTask$new(function(raw_input){
            promises::future_promise({
                raw_data <- raw_input$raw_data
                timepoint_cols <- raw_input$timepoint_cols
                datatype <- raw_input$datatype

                clean <- as.matrix(raw_data[timepoint_cols])
                rn <- switch(datatype,
                proteomics        = raw_data$Gene_Name,
                rnaseq            = raw_data$Gene_Name,
                phosphoproteomics = raw_data$pepG
                )
                rownames(clean) <- rn
                list(clean_matrix = clean)
            })
        })

        shiny::observeEvent(input[[paste0("compute_raw_ht_", tbl_name)]], {
            raw_input <- list(
                raw_data = isolate(rv$tables[[tbl_name]]),
                timepoint_cols = isolate(rv$time_cols[[tbl_name]]), 
                datatype = isolate(rv$datatype[[tbl_name]])
            )

            raw_heatmap_plot$invoke(raw_input)
        })
          
        output[[paste0("raw_ht_", tbl_name)]] <- shiny::renderPlot({
            res <- raw_heatmap_plot$result()
            shiny::req(res)
            ComplexHeatmap::Heatmap(res$clean_matrix, show_row_dend = FALSE, show_row_names = FALSE)
        })       
    })
  })
}


#' Server module for DEP heatmap + volcano
#' Expects:
#'   rv$tables[[tbl]], rv$time_cols[[tbl]], rv$datatype[[tbl]]
#'   rv$ht_matrix[[tbl]], rv$optimal_k_individual[[tbl]]
#' @export
dep_heatmap_server <- function(input, output, session, rv, cache, db_path) {

    # Track which tables have handlers registered
    registered <- reactiveVal(character())

    # Safe single numeric reader
    .read_num <- function(x) {
        if (is.null(x) || length(x) != 1) return(NA_real_)
        y <- suppressWarnings(as.numeric(x))
        if (length(y) != 1) return(NA_real_)
        y
    }

    # returns TRUE if x is an env/func we must drop
    .is_bad <- function(x) is.environment(x) || is.function(x)

    # strip envs/funs in a generic object
    sanitize_for_cache <- function(x) {
    if (.is_bad(x)) return(NULL)
    if (is.list(x)) {
        y <- lapply(x, sanitize_for_cache)
        # drop NULLs while preserving names/order
        if (!is.null(names(y))) {
            y <- y[!vapply(y, is.null, logical(1))]
        } else {
           keep <- !vapply(y, is.null, logical(1))
          y <- y[keep]
        }
        return(y)
    }
    x
    }

    # sanitize a DataFrame (rowData/colData): list-columns & metadata()
    sanitize_DataFrame <- function(df) {
        # 1) list-columns: drop env/func elements
        for (nm in colnames(df)) {
            col <- df[[nm]]
            if (is.list(col)) {
            col <- lapply(col, function(el) if (.is_bad(el)) NULL else el)
            # keep length (replace NULLs by NA_character_ to keep rectangular)
            nulls <- vapply(col, is.null, logical(1))
            if (any(nulls)) {
                # pick a sensible NA type based on the first non-NULL element
                first_good <- which(!nulls)[1]
                if (is.na(first_good)) {
                # all NULL: make it logical NA column
                col <- rep(NA, nrow(df))
                } else {
                proto <- col[[first_good]]
                if (is.numeric(proto))      col[nulls] <- list(NA_real_)
                else if (is.integer(proto)) col[nulls] <- list(NA_integer_)
                else if (is.logical(proto)) col[nulls] <- list(NA)
                else                        col[nulls] <- list(NA_character_)
                col <- I(col)  # keep as list-column
                }
            } else {
                col <- I(col)
            }
            df[[nm]] <- col
            }
        }
        # 2) metadata() on the DataFrame itself
        md <- S4Vectors::metadata(df)
        md <- sanitize_for_cache(md)
        S4Vectors::metadata(df) <- md
        df
    }

    # full SE sanitizer: top-level metadata + rowData/colData + (optional) rowRanges metadata
    sanitize_SE <- function(se) {
        # A) top-level metadata
        md <- S4Vectors::metadata(se)
        md <- sanitize_for_cache(md)
        S4Vectors::metadata(se) <- md

        # B) rowData / colData
        rd <- SummarizedExperiment::rowData(se)
        cd <- SummarizedExperiment::colData(se)
        rd <- sanitize_DataFrame(rd)
        cd <- sanitize_DataFrame(cd)
        SummarizedExperiment::rowData(se) <- rd
        SummarizedExperiment::colData(se) <- cd

        # C) (optional) rowRanges metadata if present and supported
        if ("rowRanges" %in% slotNames(se)) {
            rr <- try(SummarizedExperiment::rowRanges(se), silent = TRUE)
            if (!inherits(rr, "try-error") && !is.null(rr)) {
            mdr <- try(S4Vectors::metadata(rr), silent = TRUE)
            if (!inherits(mdr, "try-error")) {
                mdr <- sanitize_for_cache(mdr)
                try(S4Vectors::metadata(rr) <- mdr, silent = TRUE)
                try(SummarizedExperiment::rowRanges(se) <- rr, silent = TRUE)
            }
            }
        }

        se
    }

    register_for_table <- function(tbl_name) {
        #message("[", tbl_name, "] registering handlers")

        # ----- One-time UI init for gene search choices -----
        observeEvent(rv$tables[[tbl_name]], ignoreInit = FALSE, once = TRUE, {
            updateSelectizeInput(
            session, paste0("volcano_search_", tbl_name),
            choices = rv$tables[[tbl_name]]$Gene_Name, server = TRUE
            )
            updateSelectizeInput(
            session, paste0("search_gene_", tbl_name),
            choices = rv$tables[[tbl_name]]$Gene_Name, server = TRUE
            )
        })

        # ----- Build a cache key for the DEP base per table -----
        dep_key <- reactive({
            cols_key <- paste(isolate(rv$time_cols[[tbl_name]]), collapse = "_")
            paste(tbl_name, cols_key, isolate(rv$datatype[[tbl_name]]), sep = "_")
        })

        # ----- Compute / load DEP base and update contrasts -----
        observeEvent(
            list(rv$tables[[tbl_name]], rv$time_cols[[tbl_name]], rv$datatype[[tbl_name]]),
            {
                withCallingHandlers({
                    key <- dep_key()
                    dt  <- isolate(rv$datatype[[tbl_name]])
                    dat <- isolate(rv$tables[[tbl_name]])
                    req(!is.null(dat), !is.null(dt))

                    if (cache$exists(key)) {
                        #message("[", tbl_name, "] load DEP base: ", key)
                        dep_out <- cache$get(key)
                        #message("[", tbl_name, "] DEP loaded from cache")
                    } else {
                        #message("[", tbl_name, "] compute DEP base: ", key, " (", dt, ")")
                        dep_out <- switch(dt,
                            proteomics        = dep2_proteomics(dat, tbl_name, trimws(rv$time_cols[[tbl_name]])),
                            phosphoproteomics = dep2_phosphoproteomics(dat, tbl_name, trimws(rv$time_cols[[tbl_name]])),
                            rnaseq            = dep2_rnaseq(dat, tbl_name),
                            stop("[", tbl_name, "] unknown datatype: ", dt)
                        )

                        # sanitize metadata (drop envs/functions) before caching
                        dep_out <- sanitize_SE(dep_out)

                        # Preflight serialization: this mimics what storr/rds would do
                        ok <- TRUE
                        tryCatch({
                            invisible(serialize(dep_out, NULL))  # in-memory
                        }, error = function(e) {
                            ok <<- FALSE
                            #message("[", tbl_name, "] serialize(preflight) failed: ", conditionMessage(e))
                            # help identify the culprit:
                            md  <- S4Vectors::metadata(dep_out)
                            rdm <- S4Vectors::metadata(SummarizedExperiment::rowData(dep_out))
                            cdm <- S4Vectors::metadata(SummarizedExperiment::colData(dep_out))
                            # list bad fields (top level only) in each metadata
                            list_bad <- function(m) {
                                if (is.null(m)) return(character())
                                nms <- names(m); if (is.null(nms)) nms <- as.character(seq_along(m))
                                nms[vapply(m, .is_bad, logical(1))]
                            }
                            bad_top  <- list_bad(md)
                            bad_rd   <- list_bad(rdm)
                            bad_cd   <- list_bad(cdm)
                            #if (length(bad_top)) message("  top-level metadata bad: ", paste(bad_top, collapse=", "))
                            #if (length(bad_rd))  message("  rowData metadata bad: ",   paste(bad_rd,  collapse=", "))
                            #if (length(bad_cd))  message("  colData metadata bad: ",   paste(bad_cd,  collapse=", "))
                        })
                        if (!ok) stop("DEP object still contains non-serializable fields; see logs.") #else message("Preflight checks successful")

                        # Continue
                        tryCatch({
                        cache$set(key, dep_out)
                        }, error = function(e) {
                            #message("[", tbl_name, "] cache.set failed: ", conditionMessage(e))
                            md <- S4Vectors::metadata(dep_out)  # <- use S4Vectors
                            bad <- names(Filter(function(v) is.environment(v) || is.function(v), md))
                            if (length(bad)) {
                                #message("[", tbl_name, "] metadata contains env/function in: ", paste(bad, collapse = ", "))
                            }
                            stop(e)
                        })
                    }

                    # store only the key (tiny), not the big object
                    rv$dep_key[[tbl_name]] <- key

                    # derive valid contrasts and push them to the UI
                    rd_names <- colnames(SummarizedExperiment::rowData(dep_out))
                    sig_cols <- grep("_significant$", rd_names, value = TRUE)
                    choices  <- sub("_significant$", "", sig_cols)
                    rv$contrasts[[tbl_name]] <- choices

                    shinyWidgets::updateVirtualSelect(
                        session,
                        paste0("comparison_volcano_", tbl_name),
                        choices  = choices,
                        selected = if (length(choices)) choices[1] else NULL
                    )

                    # compact logging (avoid str() which returns NULL and is huge)
                    dims <- tryCatch(dim(SummarizedExperiment::assay(dep_out)), error = function(...) NULL)
                    #message("[", tbl_name, "] DEP ready ",
                    #        if (!is.null(dims)) paste0("(", paste(dims, collapse = "x"), ") ") else "",
                    #        "with ", length(choices), " contrasts")

                            # check what is causing errors
                    # store only the key (tiny), not the big object
                    step <- "set dep_key"
                    tryCatch({
                        rv$dep_key[[tbl_name]] <- key
                        #message("[", tbl_name, "] OK: ", step)
                    }, error = function(e) {
                        #message("[", tbl_name, "] FAIL: ", step, " -> ", conditionMessage(e)); stop(e)
                    })

                    # derive valid contrasts
                    rd_names <- colnames(SummarizedExperiment::rowData(dep_out))
                    sig_cols <- grep("_significant$", rd_names, value = TRUE)
                    choices  <- sub("_significant$", "", sig_cols)

                    #step <- "set contrasts"
                    #tryCatch({
                    #    # force plain character vector
                    #    choices <- as.character(choices)
                    #    rv$contrasts[[tbl_name]] <- choices
                    #    #message("[", tbl_name, "] OK: ", step, " (n=", length(choices), ")")
                    #}, error = function(e) {
                    #    #message("[", tbl_name, "] FAIL: ", step, " -> ", conditionMessage(e)); stop(e)
                    #})

                    #step <- "update virtualSelect"
                    #tryCatch({
                    #    shinyWidgets::updateVirtualSelect(
                    #        session,
                    #        paste0("comparison_volcano_", tbl_name),
                    #        choices  = choices,
                    #        selected = if (length(choices)) choices[1] else NULL
                    #    )
                    #    #message("[", tbl_name, "] OK: ", step)
                    #}, error = function(e) {
                    #    #message("[", tbl_name, "] FAIL: ", step, " -> ", conditionMessage(e)); stop(e)
                    #})

                    #step <- "final #message"
                    #tryCatch({
                    #    dims <- tryCatch(dim(SummarizedExperiment::assay(dep_out)), error = function(...) NULL)
                    #    #message("[", tbl_name, "] DEP ready ",
                    #            if (!is.null(dims)) paste0("(", paste(dims, collapse = "x"), ") ") else "",
                    #            "with ", length(choices), " contrasts")
                    #}, error = function(e) {
                    #    #message("[", tbl_name, "] FAIL: ", step, " -> ", conditionMessage(e)); stop(e)
                    #})
                    invisible(NULL)
                }, error = function(e) {
                    cat("\n=== HARD TRACE ===\n")
                    cs <- sys.calls()
                    for (i in seq_along(cs)) cat(i, ":", deparse(cs[[i]], nlines = 1L), "\n")
                    cat("=== /HARD TRACE ===\n\n")
                    stop(e)
                })
            },
            ignoreInit = FALSE
        )
        
        

        # =====================================================
        # VOLCANO — compute in future; return df + sig table
        # =====================================================
        volcano_res <- reactiveVal(NULL)

        observeEvent(input[[paste0("compute_volcano_", tbl_name)]], {
            #message("[", tbl_name, "] compute_volcano")

            key <- isolate(rv$dep_key[[tbl_name]]); req(key)
            read_num <- function(x){ 
                y <- suppressWarnings(as.numeric(x)); if (length(y)!=1 || !is.finite(y)) NA_real_ else y 
            }

            contrs <- isolate(rv$contrasts[[tbl_name]])
            if (is.null(contrs) || !length(contrs)) {
                showNotification("No contrasts available yet.", type = "error"); return()
            }

            contrast <- isolate(input[[paste0("comparison_volcano_", tbl_name)]])
            if (is.null(contrast) || !nzchar(contrast)) contrast <- contrs[1]

            pcut_raw <- read_num(isolate(input[[paste0("volcano_pcutoff_",  tbl_name)]]))
            fccut    <- read_num(isolate(input[[paste0("volcano_fccutoff_", tbl_name)]]))
            dtype    <- isolate(rv$datatype[[tbl_name]])
            req(is.finite(pcut_raw), is.finite(fccut))

            alpha <- pcut_raw; if (!is.finite(alpha) || alpha <= 0) alpha <- .Machine$double.xmin; if (alpha > 1) alpha <- 1

            promises::future_promise({
                # --- fresh, plain DB connection on worker ---
                con_w <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
                on.exit(try(DBI::dbDisconnect(con_w), silent = TRUE), add = TRUE)
                cache_w <- storr::storr_dbi(tbl_data = "storr_data", tbl_keys = "storr_keys", con = con_w)

                dep_out <- cache_w$get(key)
                dep_cut <- DEP2::add_rejections(dep_out, alpha = alpha, lfc = fccut)

                lfc_col  <- paste0(contrast, "_diff")
                pval_col <- paste0(contrast, "_p.adj")
                rd <- SummarizedExperiment::rowData(dep_cut)

                if (dtype == "proteomics") {
                    df <- data.frame(gene_names = stringr::str_to_title(rd$Gene_Name),
                                    pval = rd[[pval_col]], log2FC = rd[[lfc_col]], stringsAsFactors = FALSE)
                } else if (dtype == "phosphoproteomics") {
                    df <- data.frame(peptide = rd$pepG,
                                    pval = rd[[pval_col]], log2FC = rd[[lfc_col]], stringsAsFactors = FALSE)
                } else {
                    df <- data.frame(gene_ID = rownames(rd),
                                    pval = rd[[pval_col]], log2FC = rd[[lfc_col]], stringsAsFactors = FALSE)
                }
                df$pval <- pmax(df$pval, .Machine$double.xmin)

                sig_se    <- DEP2::get_signicant(dep_cut, contrast)
                sig_table <- DEP2::get_results(sig_se)

                list(df = df, table = sig_table, dtype = dtype)
            },
            seed     = TRUE,
            globals  = list(db_path = db_path, key = key, contrast = contrast, alpha = alpha, fccut = fccut, dtype = dtype),
            packages = c("DBI","RSQLite","storr","DEP2","SummarizedExperiment","stringr")
            ) %...>% volcano_res %...!% (function(e) {
                showNotification(conditionMessage(e), type = "error")
                #message("[", tbl_name, "] volcano error: ", conditionMessage(e))
            })
        }, ignoreInit = TRUE)

        # Native plotly render (fast on main thread)
        output[[paste0("volcano_", tbl_name)]] <- plotly::renderPlotly({
            res <- volcano_res(); req(res)
            df  <- res$df
            df$neglog10p <- -log10(df$pval)

            # Optional downsample if huge (keeps UI snappy)
            if (nrow(df) > 50000) {
                idx <- sample.int(nrow(df), 50000)
                df  <- df[idx, , drop = FALSE]
            }

            # Labels + highlighting
            text_col <- if ("gene_names" %in% names(df)) df$gene_names else if ("peptide" %in% names(df)) df$peptide else df$gene_ID
            hl <- input[[paste0("volcano_search_", tbl_name)]] |>
                    stringr::str_split(",", simplify = TRUE) |>
                    stringr::str_trim() |>
                    stringr::str_to_title()
            size <- ifelse(text_col %in% hl, 9, 5)

            plotly::plot_ly(
                df,
                x = ~log2FC,
                y = ~neglog10p,
                type = "scatter",
                mode = "markers",
                text = ~text_col,
                marker = list(size = size),
                hovertemplate = paste0(
                "%{text}<br>",
                "log2FC: %{x:.3f}<br>",
                "-log10(p): %{y:.3f}<extra></extra>"
                )
            ) %>% plotly::layout(
                xaxis = list(title = "log2FC"),
                yaxis = list(title = "-log10(p)")
            )
        }) %>% bindEvent(
        input[[paste0("compute_volcano_", tbl_name)]],
        input[[paste0("volcano_search_", tbl_name)]],
        input[[paste0("volcano_pcutoff_", tbl_name)]],
        input[[paste0("volcano_fccutoff_", tbl_name)]],
        ignoreInit = TRUE
        )

        # Sig table uses the one returned by the future (no work on main)
        output[[paste0("volcano_sig_table_", tbl_name)]] <- DT::renderDT({
            res <- volcano_res(); req(res)
            DT::datatable(
                res$table,
                extensions = "Buttons",
                options = list(
                scrollX = TRUE, pageLength = 10, dom = "Bfrtip",
                buttons = c("copy","csv","excel","pdf","print")
                )
            )
        })

        # =====================================================
        # HEATMAP: compute + DRAW in future; stream PNG to UI
        # =====================================================
        heat_res <- reactiveVal(NULL)

        observeEvent(input[[paste0("recompute_heatmap_", tbl_name)]], {
            #message("[", tbl_name, "] heatmap recompute")

            key <- isolate(rv$dep_key[[tbl_name]]); req(key)

            read_num <- function(x){ y <- suppressWarnings(as.numeric(x)); if (length(y)!=1 || !is.finite(y)) NA_real_ else y }
            pcut_raw  <- read_num(isolate(input[[paste0("heat_pcutoff_",  tbl_name)]]))
            fccut_raw <- read_num(isolate(input[[paste0("heat_fccutoff_", tbl_name)]]))
            clustering <- isTRUE(isolate(input[[paste0("clustering_",   tbl_name)]]))
            num_k      <- isolate(input[[paste0("num_clusters_", tbl_name)]])
            if (!is.finite(pcut_raw) || !is.finite(fccut_raw)) {
                showNotification("Heatmap cutoffs invalid.", type = "error"); return()
            }
            alpha <- 10^(-pcut_raw); if (!is.finite(alpha) || alpha <= 0) alpha <- .Machine$double.xmin; if (alpha > 1) alpha <- 1

            promises::future_promise(
                {
                # --- open a fresh, plain DB connection on the worker ---
                con_w <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
                on.exit(try(DBI::dbDisconnect(con_w), silent = TRUE), add = TRUE)
                cache_w <- storr::storr_dbi(tbl_data = "storr_data", tbl_keys = "storr_keys", con = con_w)

                dep_out <- cache_w$get(key)

                # Apply cutoffs & get significant set
                dep_cut <- DEP2::add_rejections(dep_out, alpha = alpha, lfc = fccut_raw)
                dep_sig <- DEP2::get_signicant(dep_cut)

                if (is.null(dep_sig) || NROW(SummarizedExperiment::rowData(dep_sig)) == 0) {
                    tf <- tempfile(fileext = ".png")
                    grDevices::png(tf, width = 1200, height = 900, res = 144)
                    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
                    plot.new(); title("No significant features for current cutoffs")
                    return(list(
                    img = list(src = tf, width = 1200, height = 900),
                    optimal_k = NA_integer_, pcut = alpha, fccut = fccut_raw, df = data.frame()
                    ))
                }

                sig_mat <- SummarizedExperiment::assay(dep_sig)

                # Compute k on significant matrix (sandbox graphics to avoid device leakage)
                opt_k <- NA_integer_
                if (!is.null(dim(sig_mat)) && nrow(sig_mat) >= 2 && ncol(sig_mat) >= 2) {
                    tfk <- tempfile(fileext = ".png")
                    grDevices::png(tfk, width = 1200, height = 900, res = 144)
                    op <- graphics::par(no.readonly = TRUE)
                    on.exit({ try(graphics::par(op), silent = TRUE); try(grDevices::dev.off(), silent = TRUE); try(unlink(tfk), silent = TRUE) }, add = TRUE)
                    elbow <- NbClust::NbClust(as.matrix(sig_mat), distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")
                    opt_k <- as.numeric(names(sort(table(elbow$Best.nc[1, ]), decreasing = TRUE)[1]))
                }

                # Selected genes table
                sig_df <- DEP2::get_results(dep_sig)

                # Draw heatmap to PNG off-main
                tf <- tempfile(fileext = ".png")
                grDevices::png(tf, width = 1200, height = 900, res = 144)
                on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
                if (isTRUE(clustering) && is.finite(num_k) && !is.na(num_k)) {
                    DEP2::plot_heatmap(dep_cut, kmeans = TRUE, k = num_k)
                } else {
                    DEP2::plot_heatmap(dep_cut)
                }

                list(img = list(src = tf, width = 1200, height = 900),
                    optimal_k = opt_k, pcut = alpha, fccut = fccut_raw, df = sig_df)
            },
            seed    = TRUE,
            globals = list(db_path = db_path, key = key, clustering = clustering, num_k = num_k, alpha = alpha, fccut_raw = fccut_raw),
            packages = c("DBI","RSQLite","storr","DEP2","SummarizedExperiment","NbClust","grDevices","graphics")
            ) %...>% heat_res %...!% (function(e) {
                showNotification(conditionMessage(e), type = "error")
                #message("[", tbl_name, "] heatmap error: ", conditionMessage(e))
            })
        }, ignoreInit = TRUE)

        # show the pre-rendered PNG
        output[[paste0("ht_", tbl_name)]] <- renderImage({
        res <- heat_res(); req(res)
        list(src = res$img$src,
            width = res$img$width,
            height = res$img$height,
            alt = "DEP heatmap")
        }, deleteFile = TRUE)

        # now uses df directly from the future
        output[[paste0("ht_sig", tbl_name)]] <- DT::renderDT({
        res <- heat_res(); req(res)
        DT::datatable(
            res$df,
            extensions = "Buttons",
            options = list(scrollX = TRUE, pageLength = 10, dom = "Bfrtip",
                        buttons = c("copy","csv","excel","pdf","print"))
        )
        })

        output[[paste0("optimal_k", tbl_name)]] <- renderUI({
        res <- heat_res(); req(res)
        if (is.null(rv$optimal_k_individual[[tbl_name]]) && !is.na(res$optimal_k)) {
            rv$optimal_k_individual[[tbl_name]] <- res$optimal_k
        }
        htmltools::tagList(
            tags$ul(tags$li(sprintf(
            "The optimal k for this table (significant subset) is: %s",
            ifelse(is.na(res$optimal_k), "N/A (too few rows)", res$optimal_k)
            )))
        )
        })
    }

    # Self-registering: runs now and on changes to rv$table_names
    observe({
        tbls <- rv$table_names
        if (is.null(tbls) || !length(tbls)) return()
        todo <- setdiff(tbls, registered())
        if (!length(todo)) return()
        for (t in todo) register_for_table(t)
        registered(union(registered(), todo))
        #message("dep_heatmap_server: registered = ", paste(registered(), collapse = ", "))
        flush.console()
    })
}