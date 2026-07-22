delete_males <- function(data) {
  data <- data |> dplyr::filter(gender == 0)
}