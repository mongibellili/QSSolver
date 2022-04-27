using StaticArrays
using BenchmarkTools

macro computeSome_thing(states, x,index ,j)
    esc(quote for $j in 1:$states
        $x[$index]+=$x[$j]*0.5 end
    end)
end
function computeSome_thing(states, x,index )
    for j = 1:states
      x[index] += x[j] * 0.5 
    end
end

function usefunction(v2::MVector{2,Float64})
    for i =1:1e+1
        computeSome_thing(2,v2,2)
    end

end
function usemacro(v2::MVector{2,Float64})
  j=1
    for i =1:1e+1
        @computeSome_thing(2,v2,2,j)
    end

end
j=1
v2 = MVector{2,Float64}(1.0,2.0)
#@computeSome_thing(2,v2,2,j)
#display(v2)


#= display(usemacro(v2,j));println()
display(v2) =#
#= display(usefunction(v2))
display(v2) =#


@btime usemacro(v2)
@btime usefunction(v2)


#= @btime computeSome_thing(2,v2,2)
@btime @computeSome_thing(2,v2,2,j)
 =#
