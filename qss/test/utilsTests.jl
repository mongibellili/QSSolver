#array3 = Array{Int64}(undef, 2)
#display(array3)
#=
arr=Array{Float64}(undef, 2)
tx = Vector{Array{Float64}}(undef, 2)
tx[1]=arr
tx[2]=Array{Float64}[]
#push!(tx[1],11)

arr[1]=11
push!(tx[1],12)
push!(tx[2],21)
push!(tx[2],22)
display(tx)
=#
#=
quantum = Vector{Float64}(undef, 2)
tx = Vector{Array{Float64}}(undef, 2)
#quantum[1]=2 # is legal because if type float vector will allocate 0.0  in first element
#tx[1]=2 # is illegal...access index 1 which does not exist
quantum[1]=-2
display(abs(quantum[1]))
#display(tx)
=#
#=
jacobian=[-2.5 20.0 1.0; 0.0 7.0 1.0]
display(jacobian);println()
display(ndims(jacobian));println()
display(size(jacobian));println()
display(size(jacobian,1));println()
=#
#jacobian=[-2.5 20.0; 0.0 7.0]
#jacobian[1,4]=1.0
#display(jacobian)
#depMatrix=Array{Float64, 2} 
#depMatrix=zeros((2, 3))
#display(depMatrix)
#depMatrix[1,1]=1.0
#=quantum = Vector{Float64}(undef, 2)
quantum[1]=2
display(quantum[1])
quantum[1]=3
display(quantum[1])
=#
x= 1>2 ? 1 : 2
display(x)
@printf "%.0f %.1f %f" 0.5 0.025 -0.0078125
