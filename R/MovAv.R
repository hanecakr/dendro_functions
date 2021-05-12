#### Moving average ----
# vectorized function
# moving average with default window-width of 11 and values assigned to central position. 
# moving average near edges of x are computed over a shrinking window.

# x vector of numerical values
# w width of window
# when w is even, include one more value from future
# align = "center" average assigned to centre of window
# align = "left" average of current value and next (w-1) values
# align = "right" average of current value and previous (w-1) values
# edges = "nofill" fill with NA when window does not include w values
# edges = "fill" average of reducing number of values at the edges of the vector

MovAv <- function(x,
                  w = 11,
                  align = "center",
                  edges = "fill") {

        if (align == "center") {
                before <- floor((w - 1) / 2)
                after  <- ceiling((w - 1) / 2)
        } else if (align == "right") {
                before <- w - 1
                after  <- 0
        } else if (align == "left") {
                before <- 0
                after  <- w - 1
        } else {
                print("'align' should be 'center', 'left' or 'right'")
        }
        
        run_mean <- matrix(NA, nrow = length(x), ncol = 1)
        n <- length(x)
        if (edges == "fill") {
                for (i in 1:n) {
                        if (is.na(x[i]) == FALSE) {
                                run_mean[i] <- mean(x[max(0, (i - before)):(i + after)], na.rm = TRUE)
                        }
                        else if (is.na(x[i]) == TRUE) {
                                run_mean[i] <- NA
                        }
                        else{
                                print("Error in dataframe. Missing values?")
                        }
                }
                
        } else if (edges == "nofill") {
                for (i in 1:n) {
                        # when x start with non-NA
                        if (i - before <= 0) {
                                run_mean[i] <- NA
                        } else {
                                run_mean[i] <-
                                        mean(x[max(0, (i - before)):(i + after)], na.rm = FALSE)
                        }
                }
                
        } else {
                print("Take care of the edges! 'Edges' should be 'nofill' or 'fill'.")
        }
        
        return(run_mean)
}
