#' Parse formula for random effects
#'
#' Parses formulas with random intercept and/or slope terms.
#' Supports:
#' - `y ~ x + (1|site)` - random intercept
#' - `y ~ x + (1 + x|site)` - correlated random intercept and slope
#' - `y ~ x + (1 + x||site)` - uncorrelated random intercept and slope
#' - `y ~ x + (0 + x|site)` or `(x|site)` - random slope only
#' - `y ~ x + (1|site/plot)` - nested random effects (expands to (1|site) + (1|site:plot))
#'
#' @param formula A formula potentially containing random effect terms
#'
#' @return A list with:
#'   - `fixed_formula`: Formula without RE terms
#'   - `has_re`: Logical, TRUE if any RE terms found
#'   - `re_terms`: List of RE specifications (one per term)
#'   - `n_re_terms`: Number of RE terms
#'   - `re_groups`: Character vector of group names (for backward compat)
#'   - `has_slopes`: Logical, TRUE if any RE term has slopes
#'   - `has_nested`: Logical, TRUE if nested RE terms were expanded
#'   - `nested_interactions`: Character vector of interaction terms to create
#'
#' @keywords internal
#'
parse_formula_re <- function(formula) {

  formula_str <- deparse(formula, width.cutoff = 500)

  # First, expand any nested RE terms: (1|a/b) -> (1|a) + (1|a:b)
  nested_result <- expand_nested_re(formula_str)
  formula_str <- nested_result$formula_str
  has_nested <- nested_result$has_nested
  nested_interactions <- nested_result$interactions

  # Pattern for RE terms: (... | group) or (... || group)
  # Captures: full match, content before bar, bar type (| or ||), group name
  # Group name can include colons for interaction terms like site:plot
  re_pattern <- "\\(\\s*([^|]+?)\\s*(\\|\\|?)\\s*([a-zA-Z_][a-zA-Z0-9_.:]*(?::[a-zA-Z_][a-zA-Z0-9_.]*)*)\\s*\\)"
  re_matches <- gregexpr(re_pattern, formula_str, perl = TRUE)

  if (re_matches[[1]][1] == -1) {
    # No random effects
    return(list(
      fixed_formula = formula,
      has_re = FALSE,
      re_terms = list(),
      n_re_terms = 0L,
      re_groups = character(0),
      has_slopes = FALSE,
      has_nested = FALSE,
      nested_interactions = character(0)
    ))
  }

  # Extract all RE term strings
  re_terms_str <- regmatches(formula_str, re_matches)[[1]]

  # Parse each RE term
  re_terms <- list()
  group_names <- character(0)
  has_slopes <- FALSE

  for (i in seq_along(re_terms_str)) {
    term_str <- re_terms_str[i]
    parsed <- parse_re_term(term_str)
    group_names <- c(group_names, parsed$group)
    parsed$term_idx <- i
    re_terms[[parsed$group]] <- parsed
    if (parsed$type != "intercept") {
      has_slopes <- TRUE
    }
  }

  # Remove RE terms from formula to get fixed-effects formula
  fixed_str <- gsub("\\+?\\s*\\([^)]*\\|\\|?[^)]+\\)", "", formula_str)
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
    re_groups = unique(group_names),
    has_slopes = has_slopes,
    has_nested = has_nested,
    nested_interactions = nested_interactions
  )
}


#' Expand nested RE syntax
#'
#' Expands `(1|a/b)` to `(1|a) + (1|a:b)` and `(1|a/b/c)` to
#' `(1|a) + (1|a:b) + (1|a:b:c)`.
#'
#' @param formula_str Formula string
#' @return List with expanded formula_str, has_nested flag, and interactions
#' @keywords internal
#'
expand_nested_re <- function(formula_str) {

  # Pattern for nested RE: (content | group1/group2/...)
  # The group part contains forward slashes
  nested_pattern <- "\\(\\s*([^|]+?)\\s*\\|\\s*([a-zA-Z_][a-zA-Z0-9_.]*(?:/[a-zA-Z_][a-zA-Z0-9_.]*)*)\\s*\\)"

  matches <- gregexpr(nested_pattern, formula_str, perl = TRUE)

  if (matches[[1]][1] == -1) {
    return(list(
      formula_str = formula_str,
      has_nested = FALSE,
      interactions = character(0)
    ))
  }

  # Extract matched terms
  matched_terms <- regmatches(formula_str, matches)[[1]]
  interactions <- character(0)
  has_nested <- FALSE

  for (term in matched_terms) {
    # Parse the term
    match <- regmatches(term, regexec(nested_pattern, term, perl = TRUE))[[1]]
    content <- match[2]  # e.g., "1" or "1 + x"
    groups <- match[3]   # e.g., "site/plot" or "site/plot/subplot"

    # Check if there's actually nesting (contains /)
    if (!grepl("/", groups)) {
      next
    }

    has_nested <- TRUE

    # Split groups
    group_parts <- strsplit(groups, "/")[[1]]

    # Build expanded terms
    expanded_terms <- character(0)
    cumulative_group <- ""

    for (i in seq_along(group_parts)) {
      if (i == 1) {
        cumulative_group <- group_parts[i]
      } else {
        cumulative_group <- paste(cumulative_group, group_parts[i], sep = ":")
        interactions <- c(interactions, cumulative_group)
      }
      expanded_terms <- c(expanded_terms, sprintf("(%s|%s)", content, cumulative_group))
    }

    # Replace the nested term with expanded terms
    replacement <- paste(expanded_terms, collapse = " + ")
    formula_str <- sub(term, replacement, formula_str, fixed = TRUE)
  }

  list(
    formula_str = formula_str,
    has_nested = has_nested,
    interactions = unique(interactions)
  )
}


#' Parse a single RE term string
#'
#' @param term_str String like "(1 + x | group)" or "(1 + x || group)"
#' @return List with group, type, vars, correlated
#' @keywords internal
#'
parse_re_term <- function(term_str) {

  # Extract parts: content | group or content || group
  # Group name can include colons for interaction terms like site:plot
  pattern <- "\\(\\s*(.+?)\\s*(\\|\\|?)\\s*([a-zA-Z_][a-zA-Z0-9_.:]*(?::[a-zA-Z_][a-zA-Z0-9_.]*)*)\\s*\\)"
  match <- regmatches(term_str, regexec(pattern, term_str, perl = TRUE))[[1]]

  if (length(match) < 4) {
    stop("Invalid random effect term: ", term_str, call. = FALSE)
  }

  content <- trimws(match[2])
  bar_type <- match[3]
  group <- match[4]

  correlated <- (bar_type == "|")

  # Parse content to extract variables
  # Possible forms:
  # "1" -> intercept only

  # "1 + x" or "1 + x + z" -> intercept + slopes
  # "0 + x" or "x" -> slope only (no intercept)

  # Split on + and clean up
  parts <- strsplit(content, "\\+")[[1]]
  parts <- trimws(parts)
  parts <- parts[parts != ""]

  has_intercept <- FALSE
  slope_vars <- character(0)

  for (p in parts) {
    if (p == "1") {
      has_intercept <- TRUE
    } else if (p == "0") {
      # Explicit no intercept
      has_intercept <- FALSE
    } else {
      # It's a variable name
      slope_vars <- c(slope_vars, p)
    }
  }

  # Determine type
  if (has_intercept && length(slope_vars) == 0) {
    type <- "intercept"
  } else if (has_intercept && length(slope_vars) > 0) {
    type <- "intercept_slope"
  } else if (!has_intercept && length(slope_vars) > 0) {
    type <- "slope"
  } else {
    stop("Invalid random effect specification: ", term_str, call. = FALSE)
  }

  list(
    group = group,
    type = type,
    has_intercept = has_intercept,
    slope_vars = slope_vars,
    correlated = correlated,
    n_coefs = as.integer(has_intercept) + length(slope_vars)
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
#'   - `has_slopes`: Whether any term has slopes
#'
#' @keywords internal
#'
build_multi_re <- function(data, re_terms) {

  if (length(re_terms) == 0) {
    return(list(
      n_re_terms = 0L,
      re_list = list(),
      total_re_params = 0L,
      has_slopes = FALSE
    ))
  }

  has_slopes <- FALSE
  re_list <- lapply(re_terms, function(term) {
    re_info <- build_re_index(data, term$group)
    re_info$group_var <- term$group
    re_info$type <- term$type
    re_info$has_intercept <- term$has_intercept
    re_info$slope_vars <- term$slope_vars
    re_info$correlated <- term$correlated
    re_info$n_coefs <- term$n_coefs

    # Build Z matrix if there are slopes
    if (term$type != "intercept") {
      has_slopes <<- TRUE
      re_info$Z <- build_re_z_matrix(data, term)
    } else {
      re_info$Z <- NULL
    }

    re_info
  })
  names(re_list) <- names(re_terms)

  # Calculate total RE parameters
  # For each term:
  #   - n_groups * n_coefs raw RE values
  #   - Variance params: correlated uses Cholesky (n_coefs*(n_coefs+1)/2)
  #                      uncorrelated uses diagonal (n_coefs)
  total_re_params <- sum(vapply(re_list, function(x) {
    n_raw <- x$n_groups * x$n_coefs
    if (x$n_coefs == 1) {
      # Just intercept: 1 log_sigma
      n_var <- 1L
    } else if (x$correlated) {
      # Correlated: Cholesky lower triangle
      n_var <- x$n_coefs * (x$n_coefs + 1L) / 2L
    } else {
      # Uncorrelated: diagonal variances
      n_var <- x$n_coefs
    }
    as.integer(n_raw + n_var)
  }, integer(1)))

  list(
    n_re_terms = length(re_list),
    re_list = re_list,
    total_re_params = total_re_params,
    has_slopes = has_slopes
  )
}


#' Build Z matrix (RE design matrix) for a term with slopes
#'
#' @param data Data frame
#' @param term RE term specification
#'
#' @return Matrix with n_obs rows and n_coefs columns
#' @keywords internal
#'
build_re_z_matrix <- function(data, term) {

  n_obs <- nrow(data)
  n_coefs <- term$n_coefs

  Z <- matrix(0, nrow = n_obs, ncol = n_coefs)
  col_idx <- 1

  # Intercept column (all 1s)
  if (term$has_intercept) {
    Z[, col_idx] <- 1
    col_idx <- col_idx + 1
  }

  # Slope columns
  for (var in term$slope_vars) {
    if (!var %in% names(data)) {
      stop(sprintf("Slope variable '%s' not found in data", var), call. = FALSE)
    }
    Z[, col_idx] <- data[[var]]
    col_idx <- col_idx + 1
  }

  Z
}
