
rho_test_perm = function(x, y, R = 1000, verbose = F, return_dist = T){
  
  rho_minus_ones = NA
  
  if(return_dist){
    rho_minus_ones = rho(x, y, return_dist = T)
  }
  
  rho_minus_one_mean = rho(x, y)
  
  rho_minus_one_perm = rep(0, R)
  
  for(i in 1:R){
    
    rho_minus_one_perm[i] = rho(sample(x), sample(y))
    
    if(verbose) cat("\r", round(100*i/R, 2), "%")
    
  }
  
  p_value = mean(rho_minus_one_perm > rho_minus_one_mean)
  
  return(list(rho_minus_one_mean = rho_minus_one_mean,
              rho_minus_one_sample = rho_minus_ones,
              p_value = p_value)
         )
  
}