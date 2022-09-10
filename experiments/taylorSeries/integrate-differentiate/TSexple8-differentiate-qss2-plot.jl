module modqss
    using TaylorSeries
    using StaticArrays
    using TimerOutputs
   use_show_default(true)
   
   function qss_Integrate(initCond::Vector{Float64},f::SVector{N, Function},order::Int)where{N}
       # reset_timer!()
       Debug=false
        ft=5.0
        absQ=1e-6
        relQ=1e-3
        dep= @SVector[[0,2],[1,2]]
        states=length(initCond)
        x=Taylor1{Float64}[]
        q=Taylor1{Float64}[]
        t=Taylor1(order-1)
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
            x[i].coeffs[k+1] =((differentiate(f[i](q,t),k-1)).coeffs[1] )/factorial(k) # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
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
        end

        if Debug 
          println("initial quantum= ",quantum)
          println("intial nextstateTime= ",nextStateTime) 
        end
        simt = initTime
        count=0
 #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
        while simt < ft #&& count < 3
            #count+=1
            simt, index = findmin(nextStateTime)
            if Debug println("loop$count : nextStateTime rightafter findmin = ",nextStateTime) end
            if Debug  println("simtime= ",simt," index=  ",index) end           
           # println(x[index](t))
           elapsed = simt - tx[index]
             #shift taylor
            for k=1:order #evaluate x at new time, and derivatives...remember coeffs store the over facto also              
                 x[index].coeffs[k] =(((differentiate(x[index],k-1)(t+elapsed))).coeffs[1])#/factorial(k-1)  # some waste here...since result of differentiate could be used in next double differentiate
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
            computeNextTime(Val(2), index, simt, nextStateTime, x, quantum) #
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
                  computeDerivative(Val(2),j, f[j], x, q, t)
                  #@timeit "Recompute next" 
                  if Debug 
                    println("loop $count : x after derivative = ",x) 
                    println("loop $count : q after derivative = ",q) 
                  end
                  reComputeNextTime(Val(2),j, simt, nextStateTime, x, q, quantum)
                  if Debug println("loop $count : nextStateTime after reComputenexttime= ",nextStateTime) end
              
                end
            end
 
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
        if (x[i].coeffs[3]) != 0
            #= println("currentTime = ",currentTime)
            println("quantum[$i] = ",quantum[i])
            println("(x[$i].coeffs[3])*2 = ",(x[i].coeffs[3])*2) =#
          nextTime[i] = currentTime + sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
           #println("schedule state at= ",nextTime[i])
        else
          nextTime[i] = Inf
        end
    end
    function computeDerivative( ::Val{2} ,j::Int, f::Function ,x::Vector{Taylor1{Float64}} , q::Vector{Taylor1{Float64}}, t::Taylor1{Float64}  )
        #x[j].coeffs[2] =((differentiate(f(q,t),0)).coeffs[1] )/factorial(1) 
        x[j].coeffs[2] =(f(q,t)).coeffs[1]
        x[j].coeffs[3] =((differentiate(f(q,t),1)).coeffs[1] )/factorial(2)
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
end#end module

#-----------user space----------------
using TaylorSeries
using StaticArrays
using BenchmarkTools
using Plots;gr()

initCond=[1.0,2.0] #[x1zero,x2zero]
order=2
Debug=false
function f1(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    q[2]#*exp(t) #cos= 1.0 - 0.5 t² + 𝒪(t³) for order 2
end
function f2(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    -q[1]-q[2]+(t)^1.9 
end
jacobian=SVector{2,Function}(f1,f2) 
sol=modqss.qss_Integrate(initCond,jacobian,order)
#= println(sol[1])
println("x1")
println(sol[2][1])
println("x2")
println(sol[2][2]) =#


if !Debug
  temp1=[]
  temp2=[]
  for i=1:length(sol[2][1])
      push!(temp1,sol[2][1][i].coeffs[1])
      push!(temp2,sol[2][2][i].coeffs[1])
  end
  display(plot!(sol[1],temp1))
  display(plot!(sol[1],temp2))
end

using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2] #+t+1
      #  du[2] = 1/(1-t)#dt <= dtmin. Aborting. There is either an error in your model specification or the true solution is unstable
     # du[2]=log(2-t)  #DomainError with -0.01799528957476415:
                      # log will only return a complex result if called with a complex argument. Try log(Complex(x)).
     #  du[2]=abs(t-2)  # ok
    #  du[2]=(1-t)^3.2   #DomainError with -0.07677066154310119:
                                    #  Exponentiation yielding a complex result requires a complex argument.
                                    # Replace x^y with (x+0im)^y, Complex(x)^y, or similar.
        du[2]=(t)^1.9 
    end
    tspan = (0.0,2)
    u0 = [1.0,2.4]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end
odeDiffEquPackage()
println("done") 
readline()
#@btime modqss.qss_Integrate(initCond,jacobian,order)

