
struct Sol
    savedTimes::Vector{Float64} 
    savedVars::Vector{Array{Taylor0{Float64}}}
end
function plotSol(sol::Sol)
    numPoints=length(sol.savedTimes)
    numVars=length(sol.savedVars)
    for k=1:numVars
      temp = []
      for i = 1:numPoints #each point is a taylor
          push!(temp, sol.savedVars[k][i].coeffs[1])
      end
     display(plot!(sol.savedTimes, temp,label="x$k")) 

    #=  if k%2==0#k=2
      display(plot!(sol.savedTimes,temp,seriestype = :scatter,label="x$k",xlims=(20,50),ylims=(0.698,0.702)))
     end =#

    # display(plot!(sol.savedTimes, temp,seriestype = :scatter,label="x$k",xlims=(10,12),ylims=(-0.5,0.5)))  
    #------------ twoVarSys1----------------
    #-----medium zoom---
    #display(plot!(sol.savedTimes, temp,label="x$k",xlims=(20,50),ylims=(-0.00003,0.00003))) 
#------------ twoVarSys13----------------
#mliqss1 adequate
#= if k%2==0
display(plot!(sol.savedTimes,label="x$k",xlims=(10,50),ylims=(-0.02,0.1))) 
else
  display(plot!(sol.savedTimes,label="x$k",line=(:dot,4),xlims=(10,50),ylims=(-0.02,0.1))) 
end =#
#liqss2 adequate
#= if k%2==0
  display(plot!(sol.savedTimes,label="x$k",xlims=(10,50),ylims=(-0.02,1.1))) #1.1 for small delta, 3.1 for a larger delta
  else
    display(plot!(sol.savedTimes,label="x$k",line=(:dot,4),xlims=(10,50),ylims=(-0.02,1.1))) 
  end =#
 # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(10.0,10.0000001),ylims=(0.0318803,0.0318810))) 
       
    
    # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(10.0,10.0000001),ylims=(0.0318803,0.0318810))) 
  # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(8.8,11.2),ylims=(-0.5,0.5))) 
    # display(plot!(sol.savedTimes, temp,label="x$k",xlims=(10.0,10.0000001),ylims=(0.0318803,0.0318810))) 
      
    end
      println("press enter to exit")
      readline()
end
function plotSol(savedTimes::Vector{Float64} ,  savedVars::Vector{Array{Taylor0{Float64}}})
    numPoints=length(savedTimes)
    numVars=length(savedVars)
    for k=1:numVars
      temp = []
      for i = 1:numPoints #each point is a taylor
          push!(temp, savedVars[k][i].coeffs[1])
      end
      display(plot!(savedTimes, temp,label="x$k"))
    end
      println("press enter to exit")
      readline()
end
function plotSol(savedTimes::Vector{Float64} , savedVar::Vector{Taylor0{Float64}})
      numPoints=length(savedTimes)
      temp = []
      for i = 1:numPoints
          push!(temp, savedVar[i].coeffs[1])
      end
      display(plot!(savedTimes, temp))
      println("press enter to exit")
      readline()
end
function evaluateSol(sol::Sol,index::Int,t::Float64)
        for i=1:length(sol[1])#savedTimes
            if sol[1][i]>t # i is closest lower point
                return sol[2][index][i-1](t-sol[1][i-1])#taylor evaluation after small elapsed with the point before (i-1)
            end
        end
end
(sol::Sol)(index::Int,t::Float64) = evaluateSol(sol,index,t)

  
function getindex(s::Sol, i::Int64)
      if i==1
         return s.savedTimes
      elseif i==2
         return s.savedVars
      else
         error("sol has 2 attributes: time and states")
      end
end