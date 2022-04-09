#array3 = Array{Int64}(undef, 2)
#display(array3)
tx = Vector{Array{Float64}}(undef, 2)
tx[1]=Array{Float64}[]
tx[2]=Array{Float64}[]



#=
push!(tx[1],11)
push!(tx[1],12)
push!(tx[2],21)
push!(tx[2],22)
=#
display(tx)