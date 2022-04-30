using BenchmarkTools


# create vectors of 500 random points each
function createvectors()
vector_points = [[0.185, 0.917], [0.695, 0.167]]
vector_points
end
function createtuples()
tuple_points  = ((0.185, 0.917), (0.695, 0.167))
tuple_points 
end
#= @show vector_points
@show tuple_points =#

function test(dep) 
    #dep=[[2],[1, 2]]
    s=0
    index=1
    for j=1:20
        index=3-index
        for i = 1:length(dep[index])
            s+=dep[index][i]
            #s+=1
        end
    end 
    s
end
#= function test1(dep)
    #dep=((2),(1, 2))
    for j=1:20
        index=rand(1:2)
        for i = 1:length(dep[index])
            #j=dep[index][i]

        end
    end    
end =#
 dep=((2,0),(1, 2))
#dep=[[2,2],[1, 2]]
#display(typeof(dep))
@btime test(dep)

#display(testTuples());println()
#display(test())
#= function testVectors()
    s=0
    for j=1:length(vector_points)
    for i = 1:length(vector_points[j])
        s += vector_points[1][i]
    end
    end
    s
end =#
#display(testVectors());println()
#@time testVectors()






#= 
rand_vector_point() = [rand(), rand()] # note the []
rand_tuple_point()  = (rand(), rand()) # note the ()

# create vectors of 500 random points each
vector_points = [rand_vector_point() for _ in 1:50]
tuple_points  = [rand_tuple_point()  for _ in 1:50]

# define a simple function calculating pairwise differences
#= function difference_matrix(points)
    [p1 .- p2 for p1 in points, p2 in points]
end =#
function difference_matrix(points)
    s=0
    for j=1:length(points)
    for i = 1:length(points[j])
        s += points[1][i]
    end
    end
end
# run each version once, just to get compilation out of the way
difference_matrix(vector_points)
difference_matrix(tuple_points)

println("Vector version:")
@btime difference_matrix(vector_points)
println("Tuple version:")
@btime difference_matrix(tuple_points) =#