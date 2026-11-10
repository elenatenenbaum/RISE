# ------------------------------------------------------------
# run_condition_model()
#
# Fits a linear mixed model and returns emmeans-based estimates
# as a formatted table (or raw data frame if return_data = TRUE).
#
# Supports:
#   - Single condition variable (e.g. Age)
#   - Multiple additive condition variables
#   - Two-way interaction (interaction = TRUE)
#
# For Age-only models, automatically prepends:
#   - Overall marginal mean (collapsing across Age)
#   - Age 12 - Age 6 contrast
#
# For interaction models, produces per-level blocks containing:
#   - Overall (marginalised over Age)
#   - Age 12 - Age 6 contrast
#   - Cell means by Age
# Followed by the omnibus interaction contrast.
#
# All p-values and CIs use adjust = "none" throughout.
#
# Inputs:
#   data           : data frame containing metric and condition vars
#   condition_var  : character vector of predictor variable name(s)
#   metric         : name of the outcome variable (character)
#   chance_level   : null hypothesis value for emmeans tests (default 0.5)
#   caption_label  : label for the kable caption (default = metric)
#   random_effects : character vector of random effect terms
#                    (default "(1 | response_uuid)")
#   interaction    : whether to fit a two-way interaction (default FALSE)
#                    requires exactly 2 condition variables;
#                    first is assumed to be Age, second the grouping factor
#   return_data    : if TRUE, return raw data frame instead of kable (default TRUE)
#
# Returns:
#   If return_data = TRUE : a data frame with columns
#                           Condition, Estimate (95% CI), P-value
#   If return_data = FALSE: a kableExtra HTML table with packed row groups
# ------------------------------------------------------------

run_condition_model <- function(
    data,
    condition_var,
    metric,
    chance_level = 0.5,
    caption_label = metric,
    random_effects = c("(1 | response_uuid)"),
    interaction = FALSE,
    return_data = TRUE
) {
  
  # ----------------------------------------------------------
  # Build formula
  # ----------------------------------------------------------
  
  if (interaction) {
    if (length(condition_var) != 2) {
      stop("interaction = TRUE requires exactly 2 condition variables.")
    }
    condition_term <- paste(condition_var, collapse = " * ")
  } else {
    condition_term <- paste(condition_var, collapse = " + ")
  }
  
  formula_text <- paste(
    metric, "~",
    condition_term, "+",
    paste(random_effects, collapse = " + ")
  )
  
  mod <- lmerTest::lmer(
    formula = as.formula(formula_text),
    data    = data,
    REML    = TRUE
  )
  
  # ----------------------------------------------------------
  # Compute emmeans grid across all condition variables
  # ----------------------------------------------------------
  
  emm_formula <- as.formula(
    paste0("~ ", paste(condition_var, collapse = " * "))
  )
  
  emm <- emmeans::emmeans(mod, emm_formula)
  
  # ----------------------------------------------------------
  # INTERACTION branch: Age x Stimulus_type (or any 2-var interaction)
  # Produces per-stimulus-type block of: Overall, Age 12 - Age 6, Age 6, Age 12
  # ----------------------------------------------------------
  
  if (interaction) {
    
    age_var  <- condition_var[1]   # first var assumed to be Age
    stim_var <- condition_var[2]   # second var is the grouping factor
    
    age_levels <- levels(factor(data[[age_var]]))
    stim_levels <- sort(unique(as.character(data[[stim_var]])))
    
    has_age_contrast <- all(c("6", "12") %in% age_levels)
    
    fmt <- function(x) sprintf("%.2f", x)
    
    fmt_pval <- function(p) {
      ifelse(p < .001, "< 0.001", as.character(round(p, 3)))
    }
    
    blocks <- lapply(stim_levels, function(stim) {
      
      # --- Cell means for this stimulus type ---
      emm_stim <- emmeans::emmeans(
        mod,
        as.formula(paste0("~ ", age_var, " | ", stim_var)),
        at = setNames(list(stim), stim_var)
      )
      
      cell_tests <- summary(emmeans::test(emm_stim, null = chance_level, adjust = "none")) %>%
        as.data.frame()
      cell_cis   <- as.data.frame(confint(emm_stim, adjust = "none")) %>%
        dplyr::select(lower.CL, upper.CL)
      
      cell_tab <- bind_cols(cell_tests, cell_cis) %>%
        mutate(
          Condition = paste0(stim, " | Age ", .data[[age_var]]),
          `Estimate (95% CI)` = paste0(
            fmt(emmean), " (", fmt(lower.CL), ", ", fmt(upper.CL), ")"
          ),
          `P-value` = fmt_pval(p.value)
        ) %>%
        dplyr::select(Condition, `Estimate (95% CI)`, `P-value`)
      
      # --- Overall for this stimulus type (marginalise over Age) ---
      emm_overall_stim <- emmeans::emmeans(
        mod,
        as.formula(paste0("~ ", stim_var)),
        at = setNames(list(stim), stim_var)
      )
      
      ov_test <- summary(emmeans::test(emm_overall_stim, null = chance_level, adjust = "none")) %>%
        as.data.frame()
      ov_ci   <- as.data.frame(confint(emm_overall_stim, adjust = "none")) %>%
        dplyr::select(lower.CL, upper.CL)
      
      overall_row <- data.frame(
        Condition           = paste0(stim, " | Overall"),
        `Estimate (95% CI)` = paste0(
          fmt(ov_test$emmean), " (", fmt(ov_ci$lower.CL), ", ", fmt(ov_ci$upper.CL), ")"
        ),
        `P-value`           = fmt_pval(ov_test$p.value),
        check.names = FALSE
      )
      
      # --- Age 12 - Age 6 contrast within this stimulus type ---
      if (has_age_contrast) {
        
        contrast_vec <- setNames(
          ifelse(age_levels == "12", 1, ifelse(age_levels == "6", -1, 0)),
          age_levels
        )
        
        contrast_result <- emmeans::contrast(
          emm_stim,
          method = list("Age 12 - Age 6" = contrast_vec), adjust = "none"
        )
        
        ct_test <- summary(contrast_result) %>% as.data.frame()
        ct_ci   <- as.data.frame(confint(contrast_result, adjust = "none")) %>%
          dplyr::select(any_of(c("lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL")))
        colnames(ct_ci) <- c("lower.CL", "upper.CL")
        
        contrast_row <- data.frame(
          Condition           = paste0(stim, " | Age 12 - Age 6"),
          `Estimate (95% CI)` = paste0(
            fmt(ct_test$estimate), " (", fmt(ct_ci$lower.CL), ", ", fmt(ct_ci$upper.CL), ")"
          ),
          `P-value`           = fmt_pval(ct_test$p.value),
          check.names = FALSE
        )
        
      } else {
        contrast_row <- NULL
      }
      
      # --- Interaction contrast (Age x Stimulus_type) ---
      # Computed once at the end outside the loop; placeholder here
      
      # Stack: Overall, Age 12 - Age 6, then age cell means (sorted)
      bind_rows(overall_row, cell_tab, contrast_row)
    })
    
    tab <- do.call(bind_rows, blocks)
    
    # --- Append the Age x group interaction contrast ---
    # revpairwise on Age -> "12 - 6" (matches the per-group "Age 12 - Age 6" rows);
    # pairwise on the grouping factor -> "Study 1 - Study 2"
    interaction_contrast <- emmeans::contrast(
      emm,
      interaction = c("revpairwise", "pairwise"),   # order follows condition_var: Age first, group second
      adjust = "none"
    )
    
    ix_test <- summary(interaction_contrast) %>% as.data.frame()
    ix_ci   <- as.data.frame(confint(interaction_contrast, adjust = "none")) %>%
      dplyr::select(any_of(c("lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL")))
    colnames(ix_ci) <- c("lower.CL", "upper.CL")
    
    age_col  <- grep("^Age_",  names(ix_test), value = TRUE)[1]   # e.g. "Age_revpairwise"
    grp_col  <- grep("pairwise$", names(ix_test), value = TRUE) |> setdiff(age_col)
    
    ix_tab <- bind_cols(ix_test, ix_ci) %>%
      mutate(
        Condition = paste0("Age ", gsub(" - ", " - Age ", .data[[age_col]]),
                           " | ", .data[[grp_col]]),
        `Estimate (95% CI)` = paste0(fmt(estimate), " (", fmt(lower.CL), ", ", fmt(upper.CL), ")"),
        `P-value` = fmt_pval(p.value)
      ) %>%
      dplyr::select(Condition, `Estimate (95% CI)`, `P-value`)
    
    tab <- bind_rows(tab, ix_tab)
    
    if (return_data) return(tab)
    
    # kable with grouping
    # rows per stim block, taken from the blocks actually built above
    # (Overall + Age contrast if present + one row per age level)
    block_n        <- vapply(blocks, nrow, integer(1))
    block_end      <- cumsum(block_n)
    block_start    <- block_end - block_n + 1
    n_stim         <- length(stim_levels)
    n_interaction  <- nrow(ix_tab)
    
    out <- tab %>%
      knitr::kable(
        format  = "html",
        caption = paste("Mixed model estimates for", caption_label)
      ) %>%
      kableExtra::kable_styling("striped")
    
    for (i in seq_along(stim_levels)) {
      out <- kableExtra::pack_rows(out, stim_levels[i], block_start[i], block_end[i])
    }
    
    out <- kableExtra::pack_rows(
      out, "Interaction contrast",
      sum(block_n) + 1,
      nrow(tab)
    )
    
    return(out)
  }
  
  # ----------------------------------------------------------
  # NON-INTERACTION branch
  # Single or additive multi-variable model:
  # returns cell means, with Overall + Age contrast prepended
  # when condition_var == "Age"
  # ----------------------------------------------------------
  
  tests <- summary(emmeans::test(emm, null = chance_level, adjust = "none"))
  cis   <- as.data.frame(confint(emm, adjust = "none"))
  
  tests <- tests %>%
    bind_cols(cis %>% dplyr::select(lower.CL, upper.CL))
  
  fmt      <- function(x) sprintf("%.2f", x)
  fmt_pval <- function(p) ifelse(p < .001, "< 0.001", as.character(round(p, 3)))
  
  tab <- tests %>%
    as.data.frame() %>%
    mutate(
      Condition = as.character(.data[[condition_var]]),
      `Estimate (95% CI)` = paste0(
        fmt(emmean), " (", fmt(lower.CL), ", ", fmt(upper.CL), ")"
      ),
      `P-value` = fmt_pval(p.value)
    ) %>%
    dplyr::select(Condition, `Estimate (95% CI)`, `P-value`)
  
  # When condition_var is "Age", prepend an Overall marginal mean
  # and an Age 12 - Age 6 contrast above the cell means
  age_overall_row  <- NULL
  age_contrast_row <- NULL
  age_prefix_tab <- NULL
  
  if (length(condition_var) == 1 && condition_var == "Age") {
    
    emm_overall  <- emmeans::emmeans(mod, ~ 1)
    overall_test <- summary(emmeans::test(emm_overall, null = chance_level, adjust = "none")) %>% as.data.frame()
    overall_ci   <- as.data.frame(confint(emm_overall, adjust = "none")) %>% dplyr::select(lower.CL, upper.CL)
    
    overall_row <- data.frame(
      Condition           = "Overall",
      `Estimate (95% CI)` = paste0(
        fmt(overall_test$emmean), " (",
        fmt(overall_ci$lower.CL), ", ",
        fmt(overall_ci$upper.CL), ")"
      ),
      `P-value`     = fmt_pval(overall_test$p.value),
      check.names   = FALSE
    )
    
    age_levels_d <- levels(factor(data[["Age"]]))
    
    if (all(c("6", "12") %in% age_levels_d)) {
      
      contrast_result <- emmeans::contrast(
        emm,
        method = list("Age 12 - Age 6" = setNames(
          ifelse(age_levels_d == "12", 1, ifelse(age_levels_d == "6", -1, 0)),
          age_levels_d
        )),
        adjust = "none"   
      )
      
      ct_test <- summary(contrast_result) %>% as.data.frame()
      ct_ci   <- as.data.frame(confint(contrast_result, adjust = "none")) %>%
        dplyr::select(any_of(c("lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL")))
      colnames(ct_ci) <- c("lower.CL", "upper.CL")
      
      contrast_row <- data.frame(
        Condition           = "Age 12 - Age 6",
        `Estimate (95% CI)` = paste0(
          fmt(ct_test$estimate), " (", fmt(ct_ci$lower.CL), ", ", fmt(ct_ci$upper.CL), ")"
        ),
        `P-value`   = fmt_pval(ct_test$p.value),
        check.names = FALSE
      )
      
    } else {
      warning("Age levels '6' and/or '12' not found; Age 12 - Age 6 contrast skipped.")
      contrast_row <- NULL
    }
    
    age_overall_row  <- overall_row
    age_contrast_row <- contrast_row
  }
  
  n_cellmean  <- nrow(tab)
  n_overall   <- if (is.null(age_overall_row))  0L else nrow(age_overall_row)
  n_contrast  <- if (is.null(age_contrast_row)) 0L else nrow(age_contrast_row)
  
  # Order: Overall, cell means by Age, then Age 12 - Age 6 last
  tab <- bind_rows(age_overall_row, tab, age_contrast_row)
  
  if (return_data) return(tab)
  
  out <- tab %>%
    knitr::kable(
      format  = "html",
      caption = paste("Mixed model estimates for", caption_label)
    ) %>%
    kableExtra::kable_styling("striped")
  
  if (n_overall > 0) {
    out <- out %>%
      kableExtra::pack_rows("Overall", 1, n_overall) %>%
      kableExtra::pack_rows("Cell means (by Age)", n_overall + 1, n_overall + n_cellmean)
    if (n_contrast > 0) {
      out <- out %>%
        kableExtra::pack_rows("Contrast", n_overall + n_cellmean + 1, n_overall + n_cellmean + n_contrast)
    }
  }
  
  out
}
