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
    
    age_levels  <- sort(unique(as.character(data[[age_var]])))
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
      bind_rows(overall_row, contrast_row, cell_tab)
    })
    
    tab <- do.call(bind_rows, blocks)
    
    # --- Append the Age x Stimulus_type interaction contrast ---
    interaction_contrast <- emmeans::contrast(emm, interaction = "pairwise", adjust = "none")
    
    ix_test <- summary(interaction_contrast) %>% as.data.frame()
    ix_ci   <- as.data.frame(confint(interaction_contrast, adjust = "none")) %>%
      dplyr::select(any_of(c("lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL")))
    colnames(ix_ci) <- c("lower.CL", "upper.CL")
    
    ix_tab <- bind_cols(ix_test, ix_ci) %>%
      mutate(
        Condition = apply(
          dplyr::select(., ends_with("pairwise")),
          1, paste, collapse = " | "
        ),
        `Estimate (95% CI)` = paste0(
          fmt(estimate), " (", fmt(lower.CL), ", ", fmt(upper.CL), ")"
        ),
        `P-value` = fmt_pval(p.value)
      ) %>%
      dplyr::select(Condition, `Estimate (95% CI)`, `P-value`)
    
    tab <- bind_rows(tab, ix_tab)
    
    if (return_data) return(tab)
    
    # kable with grouping
    n_per_stim     <- 2 + length(age_levels) + if (has_age_contrast) 1 else 0
    # rows per stim block: Overall + Age contrast (if present) + one per age level
    n_stim         <- length(stim_levels)
    n_interaction  <- nrow(ix_tab)
    
    out <- tab %>%
      knitr::kable(
        format  = "html",
        caption = paste("Mixed model estimates for", caption_label)
      ) %>%
      kableExtra::kable_styling("striped")
    
    for (i in seq_along(stim_levels)) {
      start <- (i - 1) * n_per_stim + 1
      end   <- i * n_per_stim
      out   <- kableExtra::pack_rows(out, stim_levels[i], start, end)
    }
    
    out <- kableExtra::pack_rows(
      out, "Interaction contrast",
      n_stim * n_per_stim + 1,
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
    
    age_prefix_tab <- bind_rows(overall_row, contrast_row)
  }
  
  n_prefix   <- nrow(age_prefix_tab)
  n_cellmean <- nrow(tab)
  
  if (!is.null(age_prefix_tab)) {
    tab <- bind_rows(age_prefix_tab, tab)
  }
  
  if (return_data) return(tab)
  
  out <- tab %>%
    knitr::kable(
      format  = "html",
      caption = paste("Mixed model estimates for", caption_label)
    ) %>%
    kableExtra::kable_styling("striped")
  
  if (!is.null(age_prefix_tab)) {
    out <- out %>%
      kableExtra::pack_rows("Summary", 1, n_prefix) %>%
      kableExtra::pack_rows("Cell means (by Age)", n_prefix + 1, n_prefix + n_cellmean)
  }
  
  out
}
