#' Prices of 6,283 cars listed on Autotrader
#'
#' A dataset containing the prices and other attributes of over 6000 cars in the Minneapolis area
#' @format A data frame with 6283 rows and 11 variables:
#' \describe{
#'   \item{price}{price, in US dollars}
#'   \item{Car_Info}{Raw description from website}
#'   \item{Link}{hyperlink to listing (must be appended to https://www.autotrader.com/)}
#'   \item{Make}{Car manufacturer}
#'   \item{Year}{Year car manufactured}
#'   \item{Location}{Location of listing}
#'   \item{Radius}{Radius chosen for search}
#'   \item{mileage}{mileage on vehicle}
#'   \item{status}{used/certified}
#'   \item{model}{make and model, separated by space}
#'   \item{yearsold}{2017 - autotrader$Year}
#' }
#' @source
#' \url{https://www.autotrader.com/}
#' @source
#' Suman Kundu, Yurii S. Aulchenko and A. Cecile J.W. Janssens (2020). PredictABEL:Assessment of Risk Prediction Models. R package version 1.2-4.
#' https://CRAN.R-project.org/package=PredictABEL
"autotrader"
