loglik <- function(data, param){
  
  d_tx <- data$d_tx
  e_tx <- data$e_tx
  m_tx <- param$m_tx
  
  lambda <- e_tx * m_tx
  sum(d_tx * log(lambda) - lambda)
}

loglik_param <- function(age, param){
  
  a_x <- param[1:A]
  
  m_tx <- a_x
  
  m_tx
}



