
##mother save## 
#= function save_Sol(sol::Sol{T,O},xvars::Int...)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=plot()
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if xvars[1]!=-1
    for k in xvars
      temp = []
      for i = 1:numPoints 
          push!(temp, sol.savedVars[k][i].coeffs[1])
      end
      p1=plot!(p1,sol.savedTimes, temp,marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])")
    end
  else
    for k=1:numVars
      temp = []
      for i = 1:numPoints 
          push!(temp, sol.savedVars[k][i].coeffs[1])
      end
      p1=plot!(p1,sol.savedTimes, temp,marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])")
    end
  end
  p1=plot!(p1, title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount) ")
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
end =#
function save_Sol(sol::Sol{T,O},xvars::Int...) where{T,O}# no labels
  #numPoints=length(sol.savedTimes[1])
  #numVars=length(sol.savedVars)
  p1=plot()
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if xvars[1]!=-1
    for k in xvars
      p1=plot!(p1,sol.savedTimes[k], sol.savedVars[k],marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])")
    end
  else
    for k=1:T
      p1=plot!(p1,sol.savedTimes[k], sol.savedVars[k],marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])")
    end
  end
  p1=plot!(p1, title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount) ",legend = true)
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
end
function save_Sol(sol::Sol{T,O}) where{T,O}
  save_Sol(sol,-1)
end


#if you want to put a note and lims in the graph
function save_Sol(sol::Sol{T,O},note::String,xvars::Int...;xlims::Tuple{Float64, Float64},ylims::Tuple{Float64, Float64}) where{T,O}
  p1=plot()
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  if xvars[1]!=-1
    for k in xvars
      p1=plot!(p1,sol.savedTimes[k], sol.savedVars[k],marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])")
    end
  else
    for k=1:T
      p1=plot!(p1,sol.savedTimes[k], sol.savedVars[k],marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])")
    end
  end
  if xlims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount) \n $note", xlims=xlims ,ylims=ylims)
  else
    p1=plot!(p1, title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount) \n $note")
  end
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(note)_$(timestamp).png")
end


function save_Sol(sol::Sol{T,O},note::String;xlims::Tuple{Float64, Float64},ylims::Tuple{Float64, Float64}) where{T,O}
  save_Sol(sol,note,-1;xlims=xlims,ylims=ylims)
end


function getPlot(sol::Sol{T,O}) where{T,O}
  numPoints=length(sol.savedTimes)
  #numVars=length(sol.savedVars)
  p1=plot()#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  for k=1:T
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    p1=plot!(p1,sol.savedTimes, temp,marker=(:circle),markersize=2,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.simulStepCount) ",label="x$k) steps ")
    #p1=plot!(p1,sol.savedTimes, temp,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.simulStepCount)",label="x$k"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
    #= if xlims!=(0.0,0.0) 
    p1=plot!( xlims=xlims ,ylims=ylims#= (-2.04e-1,2.06e-1) =#)
    end =#
  end
  p1
end

function getPlot(sol::Sol{T,O},index::Int) where{T,O}
  numPoints=length(sol.savedTimes)
  #numVars=length(sol.savedVars)
  #p1=plot()#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  #mydate=now()
  #timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  #for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[index][i].coeffs[1])
    end
    p1=plot!(sol.savedTimes, temp,marker=(:circle),markersize=2,#= title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.simulStepCount) ", =#label="x$index ")
    #p1=plot!(p1,sol.savedTimes, temp,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.simulStepCount)",label="x$k"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
    #= if xlims!=(0.0,0.0) 
    p1=plot!( xlims=xlims ,ylims=ylims#= (-2.04e-1,2.06e-1) =#)
    end =#
 # end
 # p1
end

#= function stack_save_Sol(sol::Sol{T,O}) where{T,O}
  numPoints=length(sol.savedTimes)
  #numVars=length(sol.savedVars)
  p1=nothing#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  for k=1:T
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    p1=plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.absQ)",label="_$(sol.algName)_x$k "#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
  end
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")  
end =#


function plotSol(sol::Sol{T,O}) where{T,O}
  numPoints=length(sol.savedTimes)
  #numVars=length(sol.savedVars)
  p1=plot()
  for k=1:T
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    display(plot!(p1,sol.savedTimes, temp,marker=(:circle),markersize=3,title="$(sol.algName)  _$(sol.absQ)_$(sol.simulStepCount)",label="x$k"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#))
  end
    println("press enter to exit")
    readline() 
end
#= 
function plot_save_SolVars(sol::Sol{T,O}) where{T,O}
  numPoints=length(sol.savedTimes)
 # numVars=length(sol.savedVars)
  p1=nothing;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  for k=1:T
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    if k%2==0#k=2
    p1=plot(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
    
    savefig(p1, "plot_$(sol.sysName)_x$(k)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
    #savefig(p1, "liqx$k.png")
    else
    p2=plot(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k"#= ,xlims=(877.5,877.8),ylims=(20.19475,20.19480) =#)
    savefig(p2, "plot_$(sol.sysName)_x$(k)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
    end
    
      
  end
 
  
    #=  display(plot(p1, p2, layout = 2))
    println("press enter to exit")
    readline()  =#
   
end



  
function stackplotSol(sol::Sol{T,O}) where{T,O}
  numPoints=length(sol.savedTimes)
 # numVars=length(sol.savedVars)
 # p1=nothing;p2=nothing
  for k=1:T
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end   
    display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.absQ)",label="$(sol.algName)_x$k"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#))
  end

    # display(plot(p1, p2, layout = 2))
    println("press enter to exit")
    readline() 
   
end

function plot_save_Sol(sol::Sol{T,O}) where{T,O}
numPoints=length(sol.savedTimes)
#numVars=length(sol.savedVars)
p1=plot()#;p2=nothing
#timestamp=string(round(now(),Minute(1)))# remove colon file name complains
mydate=now()
timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
for k=1:T
  temp = []
  for i = 1:numPoints #each point is a taylor
      push!(temp, sol.savedVars[k][i].coeffs[1])
  end
  p1=plot!(p1,sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
end
savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
display(p1)
println("press enter to exit")
readline()  
end
function stackplot_save_Sol(sol::Sol{T,O}) where{T,O}
numPoints=length(sol.savedTimes)
#numVars=length(sol.savedVars)
p1=nothing#;p2=nothing
#timestamp=string(round(now(),Minute(1)))# remove colon file name complains
mydate=now()
timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
for k=1:T
  temp = []
  for i = 1:numPoints #each point is a taylor
      push!(temp, sol.savedVars[k][i].coeffs[1])
  end
  p1=plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.absQ)",label="_$(sol.algName)_x$(k)"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
end
savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
display(p1)
println("press enter to exit")
readline()    
end


function plotSol_Der1(sol::Sol{T,O}) where{T,O}
  numPoints=length(sol.savedTimes)
  #numVars=length(sol.savedVars)
  for k=1:T
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
   display(plot!(sol.savedTimes, temp,label="x$k")) 
   temp2 = []
    for i = 1:numPoints #each point is a taylor
        push!(temp2, sol.savedVars[k][i].coeffs[2])
    end
   display(plot!(sol.savedTimes, temp2,title="$(sol.algName)",label="dx$k")) 
    
  end
    println("press enter to exit")
    readline()
end
function plotSol(savedTimes::Vector{Float64} ,  savedVars::Vector{Array{Taylor0}})
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
function plotSol(savedTimes::Vector{Float64} , savedVar::Vector{Taylor0})
      numPoints=length(savedTimes)
      temp = []
      for i = 1:numPoints
          push!(temp, savedVar[i].coeffs[1])
      end
      display(plot!(savedTimes, temp))
      println("press enter to exit")
      readline()
end =#