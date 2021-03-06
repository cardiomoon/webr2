#' Prices of 6,283 cars listed on Autotrader
#'
#' A dataset containing the prices and other attributes of over 6000 cars in the Minneapolis area
#' @format A data frame with 6283 rows and 11 variables:
#' \describe{
#'   \item{Car_Info}{Raw description from website}
#'   \item{Make}{Car manufacturer}
#'   \item{model}{make and model, separated by space}
#'   \item{Year}{Year car manufactured}
#'   \item{price}{price, in US dollars}
#'   \item{mileage}{mileage on vehicle}
#'   \item{yearsold}{2017 - autotrader$Year}
#'   \item{status}{used/certified}
#' }
#' @source
#' \url{https://www.autotrader.com/}
#' @source
#' Suman Kundu, Yurii S. Aulchenko and A. Cecile J.W. Janssens (2020). PredictABEL:Assessment of Risk Prediction Models. R package version 1.2-4.
#' https://CRAN.R-project.org/package=PredictABEL
"autotrader"
