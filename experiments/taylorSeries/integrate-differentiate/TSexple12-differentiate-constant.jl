module modqss
    using TaylorSeries
    using StaticArrays
    using TimerOutputs
   use_show_default(true)
   
   function qss_Integrate(initCond::Vector{Float64},f::SVector{N, Function},order::Int)where{N}
       # reset_timer!()
       #println("integration")
      Debug=false
      # Debug=true
        ft=15.0
        absQ=1e-6
        relQ=1e-3
        #dep= @SVector[[0,2],[1,2]]
        dep= @SVector[[0]]
        states=length(initCond)
        x=Taylor1{Float64}[]
        q=Taylor1{Float64}[]
        t=Taylor1(order)
        tt=Taylor1([0.0,0.0],2)
        #t=Taylor1([1e-120,1,0],order-1) # this is bad cuz 1/t explodes
      #  initTime=5.562684646268003e-308 ########**********sqrt and xreal exponents complain about expansion around zero....this is bad cuz solution explodes (tested)
      #  initTime= 1e-152 # only a not too small number works
        initTime=0.0
        interpolationOrder=1 # 1 means savedVars = x + t derX
        arr = []
        for i = 1:states
             push!(arr, [])
        end
        savedVars = SVector{states,Array{Taylor1{Float64}}}(tuple(arr...))    
        savedTimes = Array{Float64}([initTime])
        tx = @MVector zeros(states)
        tq = @MVector zeros(states)
        nextStateTime = @MVector zeros(states)
        nextInputTime = @MVector zeros(states)
        quantum = @MVector zeros(states)
        for i=1:states
            push!(x,Taylor1(zeros(order+1),order)+initCond[i])
            push!(q,Taylor1(zeros(order),order-1))
            push!(savedVars[i],Taylor1(zeros(interpolationOrder+1),interpolationOrder))
        end
        if Debug println("initialsavedVars  = ",savedVars) end
         for k=1:order # compute initial derivatives for x and q
           # println("k= ",k)
            for i=1:states  
                q[i].coeffs[k] =x[i].coeffs[k]
            end
          # println("k $k q= ",q) 
           for i=1:states 
            x[i].coeffs[k+1] =((differentiate(tt+f[i](q,t+initTime),k-1)).coeffs[1] )/factorial(k) # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
           end
        end     
        if Debug       
        println("initialX with derivatives= ",x)
        println("initialq with derivatives= ",q)
        end
       
        for i = 1:states
            tx[i] = initTime
            tq[i] = initTime
            for j=1:interpolationOrder
                savedVars[i][1].coeffs[j]=x[i].coeffs[j]
            end
            quantum[i] = relQ * abs(x[i].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)
            
            if quantum[i] < absQ
                quantum[i] = absQ
            end
            
            computeNextTime(Val(2), i, initTime, nextStateTime, x, quantum)
            computeNextInputTime(Val(2), i, initTime,0.0,0.0+t,f[i], nextInputTime, x,q, quantum)
        end

        if Debug 
          println("initial quantum= ",quantum)
          println("intial nextstateTime= ",nextStateTime) 
          println("intial nextInputTime= ",nextInputTime) 
        end
        simt = initTime
        count=0
 #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
        while simt < ft && count < 1555
           count+=1

           sch = updateScheduler(nextStateTime,nextInputTime)
           simt= sch[2]
            index = sch[1]

            #  simt, index = findmin(nextStateTime)
            if Debug 
              println("loop$count : nextStateTime rightafter findmin = ",nextStateTime)
              println("loop$count : nextInputTime rightafter findmin = ",nextInputTime)
              println("loop$count :  ",sch[3])
             end
            if Debug  println("simtime= ",simt," index=  ",index) end           
      # println(x[index](t))
              elapsed = simt - tx[index]
    #  if sch[3] ==:ST_STATE

           
             #shift taylor
            for k=1:order #evaluate x at new time, and derivatives...remember coeffs store the over facto also              
                 x[index].coeffs[k] =(((differentiate(x[index],k-1)(t+elapsed))).coeffs[1])/factorial(k-1)
            end
           # x[index].coeffs[1]=x[index](t+elapsed).coeffs[1]
           # x[index].coeffs[2]=((differentiate(x[index]))(t+elapsed)).coeffs[1]
          
           
            if Debug println("loop $count :x after evaluatex = ",x) end
            quantum[index] = relQ * abs(x[index].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
            if quantum[index] < absQ
                quantum[index] = absQ
            end
            if Debug  println("loop $count :quantum = ",quantum) end
            tx[index]=simt
            for k=1:order 
                q[index].coeffs[k] =x[index].coeffs[k] #updateQ
            end
            if Debug println("loop $count :q after updateQ = ",q) end
            tq[index]=simt #tq needed for higher orders
            if sch[3] ==:ST_STATE
                        computeNextTime(Val(2), index, simt, nextStateTime, x, quantum) #
                          # computeNextInputTime(Val(2), index, simt,elapsed,elapsed+t,f[index], nextInputTime, x,q, quantum)
                          if Debug println("loop $count : nextStateTime after computenexttime= ",nextStateTime) end
                                for i = 1:length(dep[index])
                              
                                j = dep[index][i] # this line before the next line or vice versa gave the same bench results
                                if j != 0   
                                  # println("dependency loop $j= ")              
                                  elapsed = simt - tx[j]                 
                                  if elapsed > 0
                                    #@timeit "integrate state" 
                                    x[j].coeffs[1]=x[j](t+elapsed).coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
                                    q[j].coeffs[1] =x[j].coeffs[1]
                                    tx[j] = simt

                                  end
                                  if Debug 
                                    println("loop $count : x before derivative = ",x) 
                                  println("loop $count : q before derivative = ",q) 
                                  end
                                  #@timeit "compute deriv" 
                                  computeDerivative(Val(2),j, f[j], x, q, simt,t)
                                  #@timeit "Recompute next" 
                                  if Debug 
                                    println("loop $count : x after derivative = ",x) 
                                    println("loop $count : q after derivative = ",q) 
                                  end
                                  reComputeNextTime(Val(2),j, simt, nextStateTime, x, q, quantum)
                                  if Debug println("loop $count : nextStateTime after reComputenexttime= ",nextStateTime) end
                              
                                end#end if j!=0
                          end#end for
            else  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself
                        # computeNextTime(Val(2), index, simt, nextStateTime, x, quantum) #
                       # computeDerivative(Val(2),index, f[index], x, q, simt+t)
                        computeNextInputTime(Val(2), index, simt,elapsed,t,f[index], nextInputTime, x,q, quantum)
                          if Debug println("under inputchange loop $count : nextInputTime after computenextInputtime= ",nextInputTime) end
                          computeDerivative(Val(2),index, f[index], x, q, simt,t)
                            
                            #@timeit "Recompute next" 
                            if Debug 
                              println("under inputchange loop $count : x after derivative = ",x) 
                              println("under input change loop $count : q after derivative = ",q) # q not updated here
                            end
                           # reComputeNextTime(Val(2),index, simt, nextStateTime, x, q, quantum)
                            for i = 1:length(dep[index])
                              
                                j = dep[index][i] # this line before the next line or vice versa gave the same bench results
                                if j != 0   
                                # println("dependency loop $j= ")              
                                  elapsed = simt - tx[j]                 
                                  if elapsed > 0
                                    #@timeit "integrate state" 
                                    x[j].coeffs[1]=x[j](t+elapsed).coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
                                    q[j].coeffs[1] =x[j].coeffs[1]
                                    tx[j] = simt
                
                                  end
                                  if Debug 
                                    println("under inputchange loop $count : x before derivative = ",x) 
                                  println("under inputchange loop $count : q before derivative = ",q) 
                                  end
                                  #@timeit "compute deriv" 
                                  computeDerivative(Val(2),j, f[j], x, q, simt,t)
                                  #@timeit "Recompute next" 
                                  if Debug 
                                    println("under inputchange loop $count : x after derivative = ",x) 
                                    println("under inputchange loop $count : q after derivative = ",q) 
                                  end
                                  reComputeNextTime(Val(2),j, simt, nextStateTime, x, q, quantum)
                                  if Debug println("under inputchange loop $count : nextStateTime after reComputenexttime= ",nextStateTime) end
                              
                                end#end if j!=0
                            end#end for
                
            end

#=           else#input 
                  computeNextInputTime(Val(2), index, simt,elapsed,elapsed+t,f[index], nextInputTime, x,q, quantum)
                 #=  for k=1:order #evaluate x at new time, and derivatives...remember coeffs store the over facto also              
                    x[index].coeffs[k] =(((differentiate(x[index],k-1)(t+elapsed))).coeffs[1])/factorial(k-1)
                  end =#
                  computeDerivative(Val(2),index, f[index], x, q, simt+t)
          end#end if state/input =#
            for k = 1:states
                temp=Taylor1(interpolationOrder)
                for j=1:interpolationOrder
                    temp.coeffs[j]=x[k].coeffs[j]
                end
               # println("temp= ",temp)
                push!(savedVars[k], temp)
            end
            push!(savedTimes, simt)
        end#end while
        (savedTimes, savedVars)
    end#end integrate

    function computeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}}, quantum::MVector{T,Float64})where{T}
      absDeltaT=1e-12
        if (x[i].coeffs[3]) != 0
            #= println("currentTime = ",currentTime)
            println("quantum[$i] = ",quantum[i])
            println("(x[$i].coeffs[3])*2 = ",(x[i].coeffs[3])*2) =#
            #nextTime[i] = currentTime + sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
            tempTime=max(sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))),absDeltaT)
            if tempTime!=absDeltaT #normal
                nextTime[i] = currentTime + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
            else#usual sqrt(quant/der) is very small
              println("(x[$i].coeffs[3])*2 = ",(x[i].coeffs[3])*2)
              x[i].coeffs[3]=sign(x[i].coeffs[3])*(abs(quantum[i])/(absDeltaT*absDeltaT))/2# adjust second derivative if it is too high
              nextTime[i] = currentTime + tempTime
              println("corrected (x[$i].coeffs[3])*2 = ",(x[i].coeffs[3])*2)
            end
           #println("schedule state at= ",nextTime[i])
        else
          nextTime[i] = Inf
        end
    end
    function computeNextInputTime(::Val{2}, i::Int, currentTime::Float64,elapsed::Float64, t::Taylor1{Float64} ,f::Function,nextInputTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}}, quantum::MVector{T,Float64})where{T}
      df=0.0
      tt=Taylor1([0.0,0.0],2)
      oldDerDerX=((x[i].coeffs[3])*2)
      newDerDerX=((differentiate(tt+f(q,currentTime+t),1)).coeffs[1] )# do not put /factorial(2) cuz here we actually need derder not store the coeff
      println("-------------------------------------------",f(q,t))
        if elapsed > 0
          println("elapsed= ",elapsed) 
          println("oldDerDerX= ",oldDerDerX) 
          println("newDerDerX= ",newDerDerX) 
          df=(newDerDerX-oldDerDerX)/(elapsed*2)
          #df=newDerDerX
          println("df= ",df) 
        else
          df= quantum[i]*1e12
        end
          
       if df!=0
        nextInputTime[i]=currentTime+(abs(quantum[i] / df))^(1/3)
        
       else
        nextInputTime[i] = Inf
        end
#=         oldDerX=((x[i].coeffs[2]))
        newDerX=((f(q,t)).coeffs[1] )
          if elapsed > 0
            df=(newDerX-oldDerX)/elapsed
            println("df= ",df) 
          else
            df= quantum[i]*1e6
          end
            
         if df!=0
          nextInputTime[i]=currentTime+(abs(quantum[i] / df))^(1/2)
         else
          nextInputTime[i] = Inf
          end =#
    end
    function computeDerivative( ::Val{2} ,j::Int, f::Function ,x::Vector{Taylor1{Float64}} , q::Vector{Taylor1{Float64}},simt::Float64, t::Taylor1{Float64}  )
        #x[j].coeffs[2] =((differentiate(f(q,t),0)).coeffs[1] )/factorial(1) 
        tt=Taylor1([0.0,0.0],2)
        x[j].coeffs[2] =(tt+f(q,t)).coeffs[1]
        x[j].coeffs[3] =((differentiate(tt+f(q,t+simt),1)).coeffs[1] )/factorial(2)
        #println("x[j].coeffs[2]= ",x[j].coeffs[2])
    end
    function reComputeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}}, quantum::MVector{T,Float64})where{T}
        coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],-(x[index].coeffs[3])*2]
       #=  println("coef[1]= ",coef[1])
        println("coef[2]= ",coef[2])
        println("coef[3]= ",coef[3]) =#
        time1 = currentTime + minPosRoot(coef, Val(2))
       # println("time1= ",time1)
        coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
       # println("coef[1]= ",coef[1])
        time2 = currentTime + minPosRoot(coef, Val(2))
      #  println("time2= ",time2)
        nextTime[index] = time1 < time2 ? time1 : time2
    end
    function minPosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
        mpr=-1 #coef1=c, coef2=b, coef3=a
        if coeff[3] == 0 #|| (10000 * abs(coeff[3])) < abs(coeff[2])
            if coeff[2] == 0
              mpr = Inf
            else 
              mpr = -coeff[1] / coeff[2]
            end
            if mpr < 0
              mpr = Inf
            end
        else 
           #double disc;
            disc = coeff[2] * coeff[2] - 4 * coeff[3] * coeff[1]#b^2-4ac
            if disc < 0 # no real roots
              mpr = Inf
            else 
              #double sd, r1;
              sd = sqrt(disc);
              r1 = (-coeff[2] + sd) / (2 * coeff[3]);
              if r1 > 0 
                mpr = r1;
              else 
                mpr = Inf;
              end
              r1 = (-coeff[2] - sd) / (2 * coeff[3]);
              if ((r1 > 0) && (r1 < mpr)) 
                mpr = r1;
              end
            end
            
        end
        return mpr
    end

    function updateScheduler(nextStateTime::MVector{T,Float64},nextInputTime :: MVector{T,Float64} )where{T}   
    
      # which is faster? finding the minimum or implementing a priority queue
      minStateTime=Inf
     # minState_index=0  # what if all nextstateTime= Inf ...especially at begining????? min_index stays 0!!!
      minState_index=1
      minInputTime=Inf
      minInput_index=0
      ST_STATE=1
      ST_INPUT=2
      for i=1:T
          if nextStateTime[i]<minStateTime
              minStateTime=nextStateTime[i]
              minState_index=i
          end
      end
      for i=1:T
          if nextInputTime[i] < minInputTime
              minInputTime=nextInputTime[i]
              minInput_index=i
          end
      end
  
  
  
  #=     if minState_index==0 
          println(" static system! all derivatives are null!")
          return (1,minTime) # later throw exception or maybe draw horizontal lines at initial conditions
      end  =#
      if minInputTime<minStateTime
         # println("an event N",minEvent_index, "about to occur! at time= ",minEventTime)
          return (minInput_index,minInputTime,:ST_INPUT)
      else
          return (minState_index,minStateTime,:ST_STATE)
      end
  
  
      
  end
  







end#end module

#-----------user space----------------
using TaylorSeries
using StaticArrays
using BenchmarkTools
using Plots;gr()

initCond=[1.0] #[x1zero,x2zero]
order=2
#Debug=true
Debug=false
function f1(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    #q[2]*30*exp(t) #cos= 1.0 - 0.5 t² + 𝒪(t³) for order 2
    #q[1]*cos(10*t)
    #cos(10*t)#alone causes base error LoadError: DomainError with Inf:sin(x) is only defined for finite x.

    #q[1]+cos(10*t)
     1 # error : if the user enter a number (constant), under the hood it should be converted to taylor constant
end
#= function f2(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    #-q[1]-q[2]
   # (1 - q[1]^2) * q[2] - q[1]  #Van der Pol
  # 1/(t+1) + sin(t)*sqrt(t)  #DomainError with Taylor1{Float64}([0.0, 1.0], 1): First non-vanishing Taylor1 coefficient must correspond
                              #to an **even power** in order to expand `sqrt` around 0
   #  abs(q[1])  #DomainError with Taylor1{Float64}([0.0, 0.0], 1):
                 # The 0th order Taylor1 coefficient must be non-zero
                  # (abs(x) is not differentiable at x=0).
                   -q[1]-q[2] +sqrt(t)#+(t)^3.2
                #  1/(3-t)
#=                 temp=0
                if 0<(1-t).coeffs[1] < 1e-6 # very small and positif, we increase it a bit
                  temp=1/(1-t+1e-5)
                elseif 0>(1-t).coeffs[1] > -1e-6 # very small and negative, we decrease it a bit
                  temp=1/(1-t-1e-5)
                else
                  temp=1/(1-t)
                end
                temp =#

end =#
 jacobian=SVector{1,Function}(f1) 
 sol=modqss.qss_Integrate(initCond,jacobian,order)
    #=println(sol[1])
    println("x1")
    println(sol[2][1])
    println("x2")
    println(sol[2][2]) =#


if !Debug
temp1=[]

for i=1:length(sol[2][1])
    push!(temp1,sol[2][1][i].coeffs[1])
   
end
display(plot!(sol[1],temp1))

end


#= using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
       # du[1] = u[1]*cos(10*t) #*30*exp(t)
       du[1]= cos(10*t)
       # du[2]=(1 - u[1]^2) * u[2] - u[1] 
      # du[2]=1/(t+1) + sin(t)*sqrt(t)
     # du[2]=sqrt(t)
    
    end
    tspan = (0.0,15)
    u0 = [0.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end
odeDiffEquPackage() =#
println("done") 
readline()
#@btime modqss.qss_Integrate(initCond,jacobian,order)

