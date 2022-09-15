

#= function getcoeffs(expr::Expr)
  smbl=Symbol(expr,:[k])
  str=String(smbl)
  myexpr=Meta.parse(str)
  return myexpr
end
macro tayArth2(expr)
  lhsrhs=postwalk(x -> x isa Expr && x.head == :ref && x.args[1] != :d ? getcoeffs(x) : x, expr)
  computedExpr=quote
                  for k=0:2
                           $lhsrhs
                          # println($lhsrhs)
                   end
              end
  esc(computedExpr)
end =#


@inline function integrateState(::Val{2}, x::Taylor0{Float64},cacheT::Taylor0{Float64},elapsed::Float64)
  x.coeffs[1] = x(elapsed)
  differentiate!(cacheT,x)
  x.coeffs[2] = cacheT(elapsed)
  cacheT.coeffs.=0.0 #clear cache
end
# later improve cases where no compatible
#@inline #did not increase performance unlike integratestate because it contains a function (p=f(q,t)) 
#= function integrate!(p::Taylor0{T}, a::Taylor0{T}, x::T) where {T<:Number}
  #p.coeffs[1]=x
  p[0]=x
  @inbounds for i = 1:a.order# improves time: with 373ns without 545ns
      p[i] = a[i-1] / i
  end
  #return nothing
end =#   #this method is needed but when using taylor0 it exist and can be seen unlike taylo1

function computeDerivative1( ::Val{2} ,j::Int,tq::MVector{T,Float64} ,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},f::Function ,cacheT::Taylor0{Float64},simt::Float64 )where{T}
 #cacheT=f[j](q,d,t)
  cacheT=f(q,d,t)
  x[j].coeffs[2] =cacheT.coeffs[1]
  x[j].coeffs[3]=cacheT.coeffs[2]/2
   return nothing
 end

 function computeDerivative2( ::Val{2} ,j::Int ,x::Vector{Taylor0{Float64}},f::Taylor0{Float64}  )
  #cacheT=f[j](q,d,t)
   
   x[j].coeffs[2] =f.coeffs[1]
   x[j].coeffs[3]=f.coeffs[2]/2
    return nothing
  end
#=  function computeDerivative2( ::Val{2} ,j::Int ,x::Vector{Taylor0{Float64}} , tt::Taylor0{Float64} ,cacheT::Taylor0{Float64},simt::Float64 )
   x[j].coeffs[3] =((differentiate(tt,1)).coeffs[1] )/2#factorial(2)
   return nothing
 end =#

 #computeDerivative(Val(2), j,tq,x,q,d,t, f[j],cacheT,simt)


 function computeDerivative( ::Val{2} ,j::Int,tq::MVector{T,Float64} ,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},f::Function ,cacheT::Taylor0{Float64},simt::Float64 )where{T}
  # println("x$j before computeDer=",x)
    x[j].coeffs[2] =f(q,d,t).coeffs[1]
  #   println("x$j after computeDer=",x)
   x[j].coeffs[3] =((differentiate(f(q,d,t),1)).coeffs[1] )/2#factorial(2)
 
   #x[j].coeffs[2] =tt(simt)# does not work because this update tt which contains an already updated one q.
   #####differentiate!(cacheT,x)
  # x[j].coeffs[3]=differentiate(tt)(simt)/2
   #println("x$j after computeDerDER=",x)
   return nothing
 end


 function computeDerivative( ::Val{2} ,j::Int ,x::Vector{Taylor0{Float64}} , tt::Taylor0{Float64} ,cacheT::Taylor0{Float64},simt::Float64 )where{T}
 # println("x$j before computeDer=",x)
   x[j].coeffs[2] =(tt).coeffs[1]
 #   println("x$j after computeDer=",x)
  x[j].coeffs[3] =((differentiate(tt,1)).coeffs[1] )/2#factorial(2)

  #x[j].coeffs[2] =tt(simt)# does not work because this update tt which contains an already updated one q.
  #####differentiate!(cacheT,x)
 # x[j].coeffs[3]=differentiate(tt)(simt)/2
  #println("x$j after computeDerDER=",x)
  return nothing
end
function differentiate!(cache::Taylor0{Float64}, a::Taylor0{Float64})
  for k in eachindex(cache)
     # differentiate!(res, a, ord)
      if k < a.order
        @inbounds cache[k] = (k+1)*a[k+1]
    end
  end
  nothing
end

#= function integrate!( a::Taylor0{T}, x::S) where {T<:Number,S<:Number}
    
  @inbounds for i in a.order-1:-1:0
      a[i+1] = a[i] / (i+1)
  end
  a[0]=x
  return nothing
end =#   #seen by taylor0

function computeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
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
            #println("(x[$i].coeffs[3])*2 = ",(x[i].coeffs[3])*2)
            x[i].coeffs[3]=sign(x[i].coeffs[3])*(abs(quantum[i])/(absDeltaT*absDeltaT))/2# adjust second derivative if it is too high
            nextTime[i] = currentTime + tempTime
           # println("corrected (x[$i].coeffs[3])*2 = ",(x[i].coeffs[3])*2)
          end
         #println("schedule state at= ",nextTime[i])
      else
        nextTime[i] = Inf
      end
      return nothing
  end
  #function computeNextInputTime(::Val{2}, i::Int, currentTime::Float64,elapsed::Float64, t::Taylor0{Float64} ,f::Function,nextInputTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::MVector{T,Float64})where{T}
    function computeNextInputTime(::Val{2}, i::Int, currentTime::Float64,elapsed::Float64, tt::Taylor0{Float64} ,nextInputTime::Vector{Float64}, x::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T<:Int64}
    df=0.0
    oldDerDerX=((x[i].coeffs[3])*2.0)
    newDerDerX=(differentiate(tt).coeffs[1] )# do not put /factorial(2) cuz here we actually need derder not store the coeff
   # println("-------------------------------------------",f(q,t))
      if elapsed > 0.0
       # println("elapsed= ",elapsed) 
       # println("oldDerDerX= ",oldDerDerX) 
       # println("newDerDerX= ",newDerDerX) 
        df=(newDerDerX-oldDerDerX)/(elapsed/2)
        #df=newDerDerX
       # println("df= ",df) 
      else
        df= quantum[i]*1e12
      end
        
     if df!=0.0
      nextInputTime[i]=currentTime+cbrt((abs(quantum[i]/df)))
      
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
        return nothing
  end
#=   function computeDerivative( ::Val{2} ,j::Int, f::Function ,x::Vector{Taylor0{Float64}} , q::Vector{Taylor0{Float64}}, t::Taylor0{Float64}  )
      #x[j].coeffs[2] =((differentiate(f(q,t),0)).coeffs[1] )/factorial(1) 
      x[j].coeffs[2] =(f(q,t)).coeffs[1]
      x[j].coeffs[3] =((differentiate(f(q,t),1)).coeffs[1] )/2#factorial(2)
      #println("x[j].coeffs[2]= ",x[j].coeffs[2])
      return nothing
  end =#

  function reComputeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
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
      return nothing
  end
 

 function computeNextEventTime(j,ZCFun::Taylor0{Float64},oldsignValue,currentTime,  nextEventTime,  quantum) #later specify args

  #ZCFun=zcFunctions[j](q,d,tq,1) # 1 for order 1
  
  #ZCFun=ZCFun0(currentTime)
  if oldsignValue[j,1] != sign(ZCFun.coeffs[1])
    nextEventTime[j]=currentTime 
    
     #=  println("old sign of",j,"  =  ",oldsignValue[j,1])
      println("new sign of",j,"  =  ",sign(ZCFun.coeffs[1])) =#
 #=      println("from computeNextEventTime:")
      println("old value of",j,"  =  ",oldsignValue[j,2])
      println("new value of",j,"  =  ",(ZCFun.coeffs[1]))
      println("schedule event now at= ",nextEventTime[j]) =#
      printFew[]=5
    
  else
    nextEventTime[j] =currentTime + minPosRoot(ZCFun.coeffs, Val(2)) #Inf  # we can estimate the time. this is important when zc depends only on time
   # nextEventTime[j]=Inf
    
  #=   println("debug")
    println("oldsignValue= ",oldsignValue )
    println("ZCFun= ", ZCFun) =#
#=     if printfewLocal>0
      println("from computeNextEventTime:")
      println("nextEventTime= ",nextEventTime)
      println("ZCFun2= ",ZCFun)
     
      global printFew[]-=1
    end =#
  
  
  
    
      
  end
  oldsignValue[j,1]=sign(ZCFun.coeffs[1])#update the values
  oldsignValue[j,2]=ZCFun.coeffs[1]
  
end


function fbarier(x::Taylor0{Float64},cacheT::Taylor0{Float64})
  x.coeffs[2] = cacheT.coeffs[1]
  x.coeffs[3] = (cacheT).coeffs[2] / 2
end