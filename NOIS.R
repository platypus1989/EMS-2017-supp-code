### The Negatively-Oriented Interval Score
### (NOIS) from Gneiting, Raftery 2007
###
### Methods for centered prediction intervals,
### and upper 95% intervals for positive data.
###
### Hunter R. Merrill
### 2016.07.20

NOIS = function(Y.obs, #N-vector of predicted values
                Y.int, #Nx2 matrix, 1st col lower, 2nd col upper, can be NA for 1-sided intervals
                sided = "upper", #one of 'center', 'upper', what type of CI
                alpha = 0.05){ #alpha for the interval ({1-alpha}100% for one-sided,  {1-alpha/2}100% for 'center')
  
  if(sided == "center"){
    
    out = mean( (Y.int[ ,2] - Y.int[ ,1]) + 
                  (2 / alpha) * (Y.int[ ,1] - Y.obs) * (Y.int[ ,1] > Y.obs) + 
                  (2 / alpha) * (Y.obs - Y.int[ ,2]) * (Y.int[ ,2] < Y.obs) )
    
  } else if(sided == "upper"){
    
    out = mean(Y.int[ ,2] + (1 / alpha) * (Y.obs - Y.int[ ,2]) * (Y.int[ ,2] < Y.obs) )
    
  } else stop("'sided' must be one of 'center', 'upper'")
  
  return(out)
}



