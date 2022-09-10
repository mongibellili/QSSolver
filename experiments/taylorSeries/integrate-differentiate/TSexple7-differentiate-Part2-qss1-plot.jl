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
            x[i].coeffs[k+1] =((differentiate(f[i](q,t),k-1)).coeffs[1] )/factorial(k) # /fact cuz i will store der/fac like the convention...to extract the derivatives multiply by fac  derderx=coef[3]*fac(2)
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
            computeNextTime(Val(1), i, initTime, nextStateTime, x, quantum)
        end
        if Debug println("intial nextstateTime= ",nextStateTime) end
        simt = initTime
        count=0
        while simt < ft #&& count < 10
           # count+=1
            simt, index = findmin(nextStateTime)
            if Debug println("loop$count : nextStateTime rightafter findmin = ",nextStateTime) end
            if Debug  println("simtime= ",simt," index=  ",index) end           
           # println(x[index](t))
           elapsed = simt - tx[index]
             #shift taylor
             for k=1:order #evaluate x at new time, and derivatives...remember coeffs store the over facto also              
                x[index].coeffs[k] =((differentiate(x[index](t+elapsed),k-1)).coeffs[1] )/factorial(k-1) # test it
           end
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
            computeNextTime(Val(1), index, simt, nextStateTime, x, quantum)
            if Debug println("loop $count : nextStateTime after computenexttime= ",nextStateTime) end
            for i = 1:length(dep[index])
                j = dep[index][i] # this line before the next line or vice versa gave the same bench results
                if j != 0
                  
                  elapsed = simt - tx[j]
                  
                  if elapsed > 0
                    #@timeit "integrate state" 
                    x[j].coeffs[1]=x[j](t+elapsed).coeffs[1] #evaluate x at new time
                    tx[j] = simt
                  end
                  #@timeit "compute deriv" 
                  computeDerivative(Val(1),j, f[j], x, q, t)
                  #@timeit "Recompute next" 
                  if Debug println("loop $count : x after derivative = ",x[j]) end
                  reComputeNextTime(Val(1),j, simt, nextStateTime, x, q, quantum)
                  if Debug println("loop $count : nextStateTime after reComputenexttime= ",nextStateTime) end
              
                end
            end
 
            for k = 1:states
                temp=Taylor1(order)
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

    function computeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}}, quantum::MVector{T,Float64})where{T}
        if (x[i].coeffs[2]) != 0
            #println(quantum[i])
          nextTime[i] = currentTime + abs(quantum[i] / (x[i].coeffs[2]))
          #println("schedule state at= ",nextTime[index])
        else
          nextTime[i] = Inf
        end
    end
    function computeDerivative( ::Val{1} ,j::Int, f::Function ,x::Vector{Taylor1{Float64}} , q::Vector{Taylor1{Float64}}, t::Taylor1{Float64}  )
        #x[j].coeffs[2] =((differentiate(f(q,t),0)).coeffs[1] )/factorial(1) 
        x[j].coeffs[2] =(f(q,t)).coeffs[1]
        #println("x[j].coeffs[2]= ",x[j].coeffs[2])
    end
    function reComputeNextTime(::Val{1}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}}, quantum::MVector{T,Float64})where{T}
        coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], -x[index].coeffs[2]]
        time1 = currentTime + minPosRoot(coef, Val(1))
        coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
        time2 = currentTime + minPosRoot(coef, Val(1))
        nextTime[index] = time1 < time2 ? time1 : time2
    end
    function minPosRoot(coeff::SVector{2,Float64}, ::Val{1}) # coming from val(1) means coef has x and derx only...size 2
        mpr=-1
            if coeff[2] == 0 
                mpr = Inf
            else 
                mpr =-coeff[1] / coeff[2];
            end
            if mpr < 0
                mpr = Inf
            end
           # println("mpr inside minPosRoot in utils= ",mpr)
        return mpr
    end
end#end module

#-----------user space----------------
using TaylorSeries
using StaticArrays
using BenchmarkTools
using Plots;gr()
using OrdinaryDiffEq
initCond=[1.0,2.0] #[x1zero,x2zero]
order=1
function f1(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    q[2]#*exp(t) #cos= 1.0 - 0.5 t² + 𝒪(t³) for order 2
end
function f2(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    -q[1]-q[2]
end
jacobian=SVector{2,Function}(f1,f2) 
sol=modqss.qss_Integrate(initCond,jacobian,order)
#= println(sol[1])
println("x1")
println(sol[2][1])
println("x2")
println(sol[2][2]) =#
temp1=[]
temp2=[]
for i=1:length(sol[2][1])
    push!(temp1,sol[2][1][i].coeffs[1])
    push!(temp2,sol[2][2][i].coeffs[1])
end
display(plot!(sol[1],temp1))
display(plot!(sol[1],temp2))









function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2] #+t+1
        du[2] = -u[1] - u[2] # +1#+10-t*t+t +cos(t)
    
    end
    tspan = (0.0,5)
    u0 = [1.0,2.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end


#odeDiffEquPackage()
println("done") 
readline()
#@btime modqss.qss_Integrate(initCond,jacobian,order)

