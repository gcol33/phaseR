#' Parse formula for random effects
#'
#' @param formula A formula potentially containing `(1 | group)` terms
#'
#' @return A list with fixed formula, RE terms, and group names
#' @keywords internal
#'
parse_formula_re <- function(formula) {

  formula_str <- deparse(formula, width.cutoff = 500)

  # Detect (1 | group) pattern
  re_pattern <- "\\(\\s*1\\s*\\|\\s*([a-zA-Z_][a-zA-Z0-9_.]*)\\s*\\)"
  re_matches <- gregexpr(re_pattern, formula_str, perl = TRUE)

  if (re_matches[[1]][1] == -1) {
    # No random effects
    return(list(
      fixed_formula = formula,
      has_re = FALSE,
      re_groups = character(0)
    ))
  }

  # Extract group names
  re_terms <- regmatches(formula_str, re_matches)[[1]]
  group_names <- gsub(re_pattern, "\\1", re_terms, perl = TRUE)

  # Remove RE terms from formula to get fixed-effects formula
  fixed_str <- gsub("\\+?\\s*\\(\\s*1\\s*\\|\\s*[^)]+\\)", "", formula_str)
  fixed_str <- gsub("~\\s*\\+", "~", fixed_str)  # Clean up leading +
  fixed_str <- trimws(fixed_str)

  if (fixed_str == "~" || fixed_str == "") {
    fixed_str <- "~ 1"
  }

  # Handle two-sided formulas
  if (grepl("^[^~]+~", fixed_str)) {
    # Has response - make sure it's preserved
    parts <- strsplit(fixed_str, "~")[[1]]
    if (length(parts) == 2 && trimws(parts[2]) == "") {
      fixed_str <- paste0(trimws(parts[1]), " ~ 1")
    }
  }

  list(
    fixed_formula = as.formula(fixed_str),
    has_re = TRUE,
    re_groups = unique(group_names)
  )
}


#' Build random effect index structure
#'
#' @param data Data frame
#' @param group_var Name of grouping variable
#'
#' @return List with group indices and counts
#' @keywords internal
#'
build_re_index <- function(data, group_var) {

  if (!group_var %in% names(data)) {
    stop(sprintf("Grouping variable '%s' not found in data", group_var),
         call. = FALSE)
  }

  groups <- as.factor(data[[group_var]])
  group_idx <- as.integer(groups)
  n_groups <- nlevels(groups)
  group_levels <- levels(groups)

  list(
    idx = group_idx,
    n_groups = n_groups,
    levels = group_levels
  )
}
