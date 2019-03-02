#' hello_country

#' @param countries string vertices countaining countries.
#' @examples
#' hello_country(countries)
#' @export
hello_country <- function(countries) {
  for (i in countries){
    print(i)
  }
}
