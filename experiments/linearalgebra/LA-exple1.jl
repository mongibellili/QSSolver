using LinearAlgebra

m=[1; 2; 3]
#display(m);println()
#display(diagm(m));println()
m2=rand(2,2)
#display(m2);println()
display(eigen(m2))
#display(diag(m2))