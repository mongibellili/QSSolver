using Plots;gr()
x = [0,1,2,3,4,5,6]; y = [3,4,4,4,2,2,0]; # These are the plotting data
#display(plot(x,y))
#display(plot!(y,x))
#readline()

tx = Vector{Array{Float64}}(undef, 2)
tx[1]=Array{Float64}[]
tx[2]=Array{Float64}[]
push!(tx[1],11)
push!(tx[1],12)
push!(tx[2],21)
push!(tx[2],22)
push!(tx[1],13)
push!(tx[1],14)
push!(tx[2],23)
push!(tx[2],24)
t = Vector{Array{Float64}}(undef, 2)
t[1]=Array{Float64}[]
t[2]=Array{Float64}[]
push!(t[1],0)
push!(t[1],1.5)
push!(t[2],0)
push!(t[2],2)
push!(t[1],2.3)
push!(t[1],3)
push!(t[2],2.1)
push!(t[2],2.6)
println(length(t[1]))
println(length(tx[1]))
#display(plot(x,y))
#readline()
for i=1:2
   # display(plot!(t[i],tx[i]))
end
readline()