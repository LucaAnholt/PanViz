#' multi_hex_col_mix
#'
#' @description This is a helper function that merges any vector of hex colours
#' @param col_vector - vector of hex colours
#'
#' @return - a single mixed hex color from inputted hex codes
#'
multi_hex_col_mix <- function(col_vector){
  rgb <- as.vector(grDevices::col2rgb(col_vector))
  rgb <- as.data.frame(matrix(rgb, ncol = 3,  byrow = TRUE), stringsAsFactors = FALSE)
  mix <- colorspace::mixcolor(0.5, colorspace::RGB(R = rgb$V1[1], G = rgb$V2[1], B = rgb$V3[1]), colorspace::RGB(R = rgb$V1[2], G = rgb$V2[2], B = rgb$V3[2]))
  count <- as.numeric(seq(1, length(rgb$V1))[-c(1,2)])
  for(i in count){
    mix <- colorspace::mixcolor(0.5, mix, colorspace::RGB(R = rgb$V1[i], G = rgb$V2[i], B = rgb$V3[i]))
  }
  mix <- c(mix@coords[1], mix@coords[2], mix@coords[3])
  mix <- unique(grDevices::rgb(red = mix[1], green = mix[2], blue = mix[3], maxColorValue = 255))
  return(mix)
}

