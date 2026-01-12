#' Parse formula for random effects
#'
#' Parses formulas with one or more random intercept terms.
#' Supports: `y ~ x + (1|site) + (1|year)`
#'
#' @param formula A formula potentially containing `(1 | group)` terms
#'
#' @return A list with:
#'   - `fixed_formula`: Formula without RE terms
#'   - `has_re`: Logical, TRUE if any RE terms found
#'   - `re_terms`: List of RE specifications (one per term)
#'   - `n_re_terms`: Number of RE terms
#'   - `re_groups`: Character vector of group names (for backward compat)
#'
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
      re_terms = list(),
      n_re_terms = 0L,
      re_groups = character(0)
    ))
  }

  # Extract group names
  re_terms_str <- regmatches(formula_str, re_matches)[[1]]
  group_names <- gsub(re_pattern, "\\1", re_terms_str, perl = TRUE)

  # Build structured RE term list
  re_terms <- lapply(seq_along(group_names), function(i) {
    list(
      group = group_names[i],
      type = "intercept",
      term_idx = i
    )
  })
  names(re_terms) <- group_names

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
    re_terms = re_terms,
    n_re_terms = length(re_terms),
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


#' Build multiple random effect structures
#'
#' @param data Data frame
#' @param re_terms List of RE term specifications from parse_formula_re()
#'
#' @return List with:
#'   - `n_re_terms`: Number of RE terms
#'   - `re_list`: List of RE structures (one per term)
#'   - `total_re_params`: Total number of RE parameters
#'
#' @keywords internal
#'
build_multi_re <- function(data, re_terms) {

 if (length(re_terms) == 0) {
    return(list(
      n_re_terms = 0L,
      re_list = list(),
      total_re_params = 0L
    ))
  }

  re_list <- lapply(re_terms, function(term) {
    re_info <- build_re_index(data, term$group)
    re_info$group_var <- term$group
    re_info$type <- term$type
    re_info
  })
  names(re_list) <- names(re_terms)

  # Total RE params: sum of n_groups + one sigma per term
  total_re_params <- sum(vapply(re_list, function(x) x$n_groups + 1L, integer(1)))

  list(
    n_re_terms = length(re_list),
    re_list = re_list,
    total_re_params = total_re_params
  )
}
