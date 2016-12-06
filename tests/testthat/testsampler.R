library(associationr)

context("testing sampler")

test_that("sample.phenotype", {
  n = 10
  L = 500
  G = matrix(sample(c(0,1),n*L,replace = TRUE),n,L)
  res = sample.phenotype(G,0.2,100)

  #res$outlier
  #res$X

})
