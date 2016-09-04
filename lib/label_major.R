#' Label the major locations in Alaska
#'
#' @param font The font to print the text in.
#' @param colour The colour of the font.
#' @param size The font size to print the text in.
#'
label_major <- function(font, colour = "black", size) {
  text(x = 200000, y = 1800000, "Alaska",
    font = font, col = colour, cex = size)
  text(x = -1000000, y = 1300000, "Bering \n Sea",
    font = font, col = colour, cex = size)
  text(x = -1200000, y = 250000, "Aleutian Islands",
    font = font, col = colour, cex = size)
  text(x = 600000, y = 500000, "Gulf of \n Alaska",
    font = font, col = colour, cex = size)
}
