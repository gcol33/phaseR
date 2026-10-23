# Create interaction variables for nested RE

Creates columns like "site:plot" from "site" and "plot" columns.

## Usage

``` r
create_nested_interactions(data, interactions)
```

## Arguments

- data:

  Data frame

- interactions:

  Character vector of interaction terms (e.g., "site:plot")

## Value

Modified data frame with new columns
