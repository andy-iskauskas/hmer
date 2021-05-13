# Implausibility colour palettes
redgreen <- c('#00FF00', '#18FF00', '#31FF00', '#49FF00', '#62FF00', '#7AFF00', '#93FF00', '#ABFF00', '#C4FF00', '#DDFF00',
              '#E0E200', '#E4C600', '#E8AA00', '#EC8D00', '#EF7100', '#F35500', '#F73800', '#FB1C00', '#FF0000', '#FF0000')
colourblind <- c('#1aff1a', '#2af219', '#3ae618', '#4ada17', '#5acd16', '#6bc115', '#7bb514', '#8ba813', '#9b9c12', '#ac9011',
                 '#a2831d', '#98762a', '#8e6936', '#845d43', '#7a504f', '#70435c', '#663768', '#5c2a75', '#521d81', '#48118e')
redgreencont <- list(low = '#00FF00', mid = '#DDFF00', high = '#FF0000')
colourblindcont <- list(low = '#1aff1a', mid = '#ac9011', high = '#48118e')

# Scales inputs: important since the emulators should take inputs purely in [-1,1]
scale_input <- function(x, r, forward = TRUE) {
  centers <- purrr::map(r, ~(.x[2]+.x[1])/2)
  scales <- purrr::map(r, ~(.x[2]-.x[1])/2)
  if (is.null(names(x))) {
    centers <- unlist(centers, use.names = F)
    scales <- unlist(scales, use.names = F)
  }
  if (forward)
    output <- (x-centers)/scales
  else
    output <- x * scales + centers
  if (!"data.frame" %in% class(output))
    return(data.frame(output))
  return(output)
}

# Helper to convert functions to names
function_to_names <- function(f, var_names, function_form = TRUE) {
  f_str <- deparse(body(f))
  subbed <- stringr::str_replace_all(gsub("x\\[\\[(\\d*)\\]\\](^\\d)?", "var_names[\\1]\\2", f_str), "var_names\\[\\d*\\]", function(x) eval(parse(text = x)))
  if (function_form) {
    subbed <- gsub("\\s\\*\\s", ":", subbed)
  }
  else {
    subbed <- gsub("\\s\\*\\s", "\\*", subbed)
  }
  return(subbed)
}

# Evaluate multiple functions over points
eval_funcs <- function(funcs, points, ...) {
  pointsdim <- (length(dim(points)) != 0)
  manyfuncs <- (typeof(funcs) != "closure")
  if (manyfuncs && pointsdim)
    return(apply(points, 1, function(x) purrr::map_dbl(funcs, purrr::exec, x, ...)))
  if (manyfuncs)
    return(purrr::map_dbl(funcs, purrr::exec, points, ...))
  if (pointsdim) {
    return(tryCatch(apply(points, 1, funcs, ...),
                    error = function(cond1) {
                      tryCatch(purrr::exec(funcs, points, ...),
                               error = function(cond2) {
                                 print(cond1, cond2)
                               })
                    }))
  }
  return(purrr::exec(funcs, points, ...))
}
