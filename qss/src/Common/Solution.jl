
struct Sol{T}
    savedTimes::Vector{Float64} 
    savedVars::Vector{Array{Taylor0{Float64}}}
    algName::String
    sysName::String
    numSteps ::MVector{T,Int}
    absQ::Float64
end
function getError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  sumTrueSqr=0.0
  sumDiffSqr=0.0
  relerror=0.0
  if index<=numVars
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      sumDiffSqr+=(sol.savedVars[index][i].coeffs[1]-ft)*(sol.savedVars[index][i].coeffs[1]-ft)
      sumTrueSqr+=ft*ft
    end
    relerror=sqrt(sumDiffSqr/sumTrueSqr)
  else
    error("the system contains only $numVars variables!")
  end
  return relerror
end

function plotRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  if index<=numVars
    temp = []
    tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      if ft>1e-12 || ft < -1e-12
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)/ft))
        push!(tempt,sol.savedTimes[i])
      end
    end
  # display(plot!(tempt, temp,title="RelError:(S-T)/T for x$index",label="$(sol.algName)")) 
  display(plot!(tempt, temp,title="RelError:(S-T)/T for x$index",label="$(sol.algName)",xlims=(80,200),ylims=(0.0,0.0002)) )
   
  else
    error("the system contains only $numVars variables!")
  end
  println("press enter to exit")
  readline()
end
function plotAbsoluteError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  if index<=numVars
    temp = []
   # tempt = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      #if ft>1e-2 || ft < -1e-2
        push!(temp, abs((sol.savedVars[index][i].coeffs[1]-ft)))
       # push!(tempt,sol.savedTimes[i])
     # end
    end
   #display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)")) 
   display(plot!(sol.savedTimes, temp,title="AbsError:(S-T) for x$index",label="$(sol.algName)",xlims=(80,200),ylims=(0.0,0.0002))) 
   
  else
    error("the system contains only $numVars variables!")
  end
  println("press enter to exit")
  readline()
end
function plotCumulativeSquaredRelativeError(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  sumTrueSqr=0.0
  sumDiffSqr=0.0
  
  if index<=numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      sumDiffSqr+=(sol.savedVars[index][i].coeffs[1]-ft)*(sol.savedVars[index][i].coeffs[1]-ft)
      sumTrueSqr+=ft*ft
        push!(temp, sqrt(sumDiffSqr/sumTrueSqr))
    end
   display(plot!(sol.savedTimes, temp,title="Error:sqrt(∑(S-T)^2/∑T^2) for x$index",label="$(sol.algName)")) 
  else
    error("the system contains only $numVars variables!")
  end
  println("press enter to exit")
  readline()
end
function plotMSE(sol::Sol,index::Int,f::Function)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
 # sumTrueSqr=0.0
  sumDiffSqr=0.0
  
  if index<=numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
      ft=f(sol.savedTimes[i])
      sumDiffSqr+=(sol.savedVars[index][i].coeffs[1]-ft)*(sol.savedVars[index][i].coeffs[1]-ft)
     # sumTrueSqr+=ft*ft
        push!(temp, (sumDiffSqr/i))
    end
   display(plot!(sol.savedTimes, temp,title="Error:(∑(S-T)^2/i) for x$index",label="$(sol.algName)")) 
  else
    error("the system contains only $numVars variables!")
  end
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


function plotSol(sol::Sol)
    numPoints=length(sol.savedTimes)
    numVars=length(sol.savedVars)
   # p1=nothing;p2=nothing
    for k=1:numVars
      temp = []
      for i = 1:numPoints #each point is a taylor
          push!(temp, sol.savedVars[k][i].coeffs[1])
      end
      if k%2==0#k=2
      plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
     # savefig(p1, "x2_$(sol.algName)_.png")
      #savefig(p1, "liqx$k.png")
      else
      plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k"#= ,xlims=(877.5,877.8),ylims=(20.19475,20.19480) =#)
     # savefig(p2, "x1_$(sol.algName)_.png")
      end
      
   #  display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k"#= ,xlims=(435,530),ylims=(0,1) =#))
   #display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k",xlims=(62.2,62.5),ylims=(9.8,10.35)))  #both medium zoom
  # display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k",xlims=(63.2,63.4),ylims=(9.9,10.3))) #mliqss2 zoom
   #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(435,800),ylims=(-0.02,0.1)))
   #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(800,2000),ylims=(-0.02,0.1)))#+1e-19
  # display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(877.5,877.8),ylims=(-0.01,0.03)))#++1e-19
 #=   println("press enter to exit")
   readline() =#
   #x1
 #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(918,920),ylims=(20.1889,20.1890))) #+eps=1e-6
 #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(700,2000),ylims=(20.17,20.21)))#++1e-19 
 #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(842,850),ylims=(20.1873,20.1874)))#++1e-19
 #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(847.35,847.65),ylims=(20.1873,20.1874)))#+mliqss2+1e-19
 #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(877.5,877.8),ylims=(20.19475,20.19480)))#+liqss2+1e-19
   #  display(plot!(sol.savedTimes, temp,title="System1-qss2",label="x$k",xlims=(0,0.6),ylims=(-0.4,1))) 
    #=  if k%2==1#k=1 sys5 x1
      display(plot!(sol.savedTimes,temp#= ,line=(:dot,3) =#,marker=(:circle)#= marker=([:cross :d]) =#,title="$(sol.algName)",label="$(sol.algName): x1",xlims=(80,200),ylims=(0.7495,0.7505))) # system5 against true solution
     end =#
    #=  if k%2==0#k=2 sys5 x2
      display(plot!(sol.savedTimes,temp#= ,line=(:dot,3) =#,marker=(:circle)#= marker=([:cross :d]) =#,title="$(sol.algName)",label="$(sol.algName): x2",xlims=(80,200),ylims=(3.999,4.0020))) # system5 against true solution
     end =#
    #=  if k%2==0#k=2 sys3 x2
      display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =##= line=(:dot,3), =#marker=(:circle),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(7,50),ylims=(0.6998,0.7003)))
     # display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =#line=(:dot,3),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(2,20),ylims=(0.6991,0.701)))
     # display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =#line=(:dot,3),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(2,15),ylims=(0.698,0.81)))
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
   
    
      # display(plot(p1, p2, layout = 2))
      println("press enter to exit")
      readline() 
     
end
function plot_save_SolVars(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=nothing;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate)))
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
    if k%2==0#k=2
    p1=plot(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
    
    savefig(p1, "plot_$(sol.sysName)_x$(k)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
    #savefig(p1, "liqx$k.png")
    else
    p2=plot(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(20.19475,20.19480) =#)
    savefig(p2, "plot_$(sol.sysName)_x$(k)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
    end
    
      #  display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k"#= ,xlims=(435,530),ylims=(0,1) =#))
      #display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k",xlims=(62.2,62.5),ylims=(9.8,10.35)))  #both medium zoom
      # display(plot!(sol.savedTimes,marker=(:circle), temp,title="$(sol.algName)",label="$(sol.algName):-x$k",xlims=(63.2,63.4),ylims=(9.9,10.3))) #mliqss2 zoom
      #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(435,800),ylims=(-0.02,0.1)))
      #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(800,2000),ylims=(-0.02,0.1)))#+1e-19
      # display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(877.5,877.8),ylims=(-0.01,0.03)))#++1e-19
      #=   println("press enter to exit")
      readline() =#
      #x1
      #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(918,920),ylims=(20.1889,20.1890))) #+eps=1e-6
      #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(700,2000),ylims=(20.17,20.21)))#++1e-19 
      #  display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(842,850),ylims=(20.1873,20.1874)))#++1e-19
      #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(847.35,847.65),ylims=(20.1873,20.1874)))#+mliqss2+1e-19
      #display(plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.algName)",label="x$k",xlims=(877.5,877.8),ylims=(20.19475,20.19480)))#+liqss2+1e-19
      #  display(plot!(sol.savedTimes, temp,title="System1-qss2",label="x$k",xlims=(0,0.6),ylims=(-0.4,1))) 
        #=  if k%2==1#k=1 sys5 x1
          display(plot!(sol.savedTimes,temp#= ,line=(:dot,3) =#,marker=(:circle)#= marker=([:cross :d]) =#,title="$(sol.algName)",label="$(sol.algName): x1",xlims=(80,200),ylims=(0.7495,0.7505))) # system5 against true solution
        end =#
        #=  if k%2==0#k=2 sys5 x2
          display(plot!(sol.savedTimes,temp#= ,line=(:dot,3) =#,marker=(:circle)#= marker=([:cross :d]) =#,title="$(sol.algName)",label="$(sol.algName): x2",xlims=(80,200),ylims=(3.999,4.0020))) # system5 against true solution
        end =#
        #=  if k%2==0#k=2 sys3 x2
          display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =##= line=(:dot,3), =#marker=(:circle),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(7,50),ylims=(0.6998,0.7003)))
        # display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =#line=(:dot,3),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(2,20),ylims=(0.6991,0.701)))
        # display(plot!(sol.savedTimes,temp,#= seriestype = :scatter =#line=(:dot,3),title="$(sol.algName)",label="$(sol.algName): x2",xlims=(2,15),ylims=(0.698,0.81)))
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
 
  
    #=  display(plot(p1, p2, layout = 2))
    println("press enter to exit")
    readline()  =#
   
end
function plot_save_Sol(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=nothing#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate)))
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
  
    p1=plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
    
    
    #savefig(p1, "liqx$k.png")
    
    
    
  end
 
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
     display(p1)
    println("press enter to exit")
    readline() 
   
end
function save_Sol(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  p1=nothing#;p2=nothing
  #timestamp=string(round(now(),Minute(1)))# remove colon file name complains
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate)))
  for k=1:numVars
    temp = []
    for i = 1:numPoints #each point is a taylor
        push!(temp, sol.savedVars[k][i].coeffs[1])
    end
  
    p1=plot!(sol.savedTimes, temp,marker=(:circle),title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)",label="x$k:$(sol.numSteps[k]) steps"#= ,xlims=(877.5,877.8),ylims=(-0.01,0.03) =#)
    
    
    #savefig(p1, "liqx$k.png")
    
    
    
  end
 
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(timestamp).png")
    #=  display(p1)
    println("press enter to exit")
    readline()  =#
   
end
function plotSol_Der1(sol::Sol)
  numPoints=length(sol.savedTimes)
  numVars=length(sol.savedVars)
  for k=1:numVars
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


  
function getindex(s::Sol, i::Int64)
      if i==1
         return s.savedTimes
      elseif i==2
         return s.savedVars
      else
         error("sol has 2 attributes: time and states")
      end
end