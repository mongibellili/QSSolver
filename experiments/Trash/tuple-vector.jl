using BenchmarkTools
rand_vector_point() = [rand(), rand()] # note the []
rand_tuple_point()  = (rand(), rand()) # note the ()

# create vectors of 500 random points each
vector_points = [rand_vector_point() for _ in 1:2]
tuple_points  = [rand_tuple_point()  for _ in 1:2]
@show vector_points
#@show tuple_points
# define a simple function calculating pairwise differences
function difference_matrix(points)
    #[p1 .- p2 for p1 in points, p2 in points]
    [p1  for p1 in points]
end

# run each version once, just to get compilation out of the way
@show difference_matrix(vector_points)
println()
#@show difference_matrix(tuple_points)


#@btime difference_matrix(vector_points)

#@btime difference_matrix(tuple_points)