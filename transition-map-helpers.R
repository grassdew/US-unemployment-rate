#
# Helper function to output frame index with leading zeros. Used for filenames.
#
getFrameIndex <- function(x) {
  if (x == 0) {
    frame_index <- "00000"
  } else if (x < 10) {
    frame_index <- paste("0000", x, sep="")
  } else if (x < 100) {
    frame_index <- paste("000", x, sep="")
  } else if (x < 1000) {
    frame_index <- paste("00", x, sep="")
  } else if (x < 10000) {
    frame_index <- paste("0", x, sep="")
  } else {
    frame_index <- x
  }
  
  return(frame_index)
}