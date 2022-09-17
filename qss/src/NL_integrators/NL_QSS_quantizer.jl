@inline function integrateState(::Val{2}, x::Taylor0{Float64},cacheT::Taylor0{Float64},elapsed::Float64) #@inline increased performance
  x.coeffs[1] = x(elapsed)
  differentiate!(cacheT,x)
  x.coeffs[2] = cacheT(elapsed)
  #cacheT.coeffs.=0.0 #clear cache
end


function computeDerivative1( ::Val{2} ,j::Int,tq::MVector{T,Float64} ,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},f::Function ,cacheT::Taylor0{Float64},simt::Float64 )where{T}
  cacheT=f(q,d,t)
  x[j].coeffs[2] =cacheT.coeffs[1]
  x[j].coeffs[3]=cacheT.coeffs[2]/2
   return nothing
 end

 function computeDerivative2( ::Val{2}  ,x::Taylor0{Float64},f::Taylor0{Float64}  )
  #cacheT=f[j](q,d,t)
   
   x.coeffs[2] =f.coeffs[1]
   x.coeffs[3]=f.coeffs[2]/2
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
    #newDerDerX=(differentiate(tt).coeffs[1] )# do not put /factorial(2) cuz here we actually need derder not store the coeff
    newDerDerX=(tt.coeffs[2])
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

#=   function reComputeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
      #coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],-(x[index].coeffs[3])*2]
     #=  println("coef[1]= ",coef[1])
      println("coef[2]= ",coef[2])
      println("coef[3]= ",coef[3]) =#
      if  1.4387457897838291<=currentTime < 1.43874578978382928
      @show x
      @show q
      @show quantum
      end
      coef=NTuple{3,Float64}((-(x[index].coeffs[3])*2,q[index].coeffs[2]-x[index].coeffs[2],q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index]))
      time1 = currentTime + minPosRoot(coef, Val(2))
      if  1.4387457897838291<=currentTime < 1.43874578978382928
        println("currenttime= ",currentTime)
        println("minPosRoot(coef, Val(2))= ",minPosRoot(coef, Val(2)))
       end
      coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],3)
     # println("coef[1]= ",coef[1])
      time2 = currentTime + minPosRoot(coef, Val(2))
    #  println("time2= ",time2)
    if  1.4387457897838291<=currentTime < 1.43874578978382928
      println("2nd inPosRoot= ",minPosRoot(coef, Val(2)))
     end
      nextTime[index] = time1 < time2 ? time1 : time2
      if  1.4387457897838291<=currentTime < 1.43874578978382928
        println("nextTime[$index]= ",nextTime[index])
       end
      return nothing
  end =#
 
  function reComputeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],-(x[index].coeffs[3])*2]
   #=  println("coef[1]= ",coef[1])
    println("coef[2]= ",coef[2])
    println("coef[3]= ",coef[3]) =#
    #= if  1.4387457897838291<=currentTime < 1.43874578978382928
      @show x
      @show q
      @show quantum
      end =#
    time1 = currentTime + minPosRoot(coef, Val(2))
#=     if  1.4387457897838291<=currentTime < 1.43874578978382928
      println("currenttime= ",currentTime)
      println("minPosRoot(coef, Val(2))= ",minPosRoot(coef, Val(2)))
     end =#
   # println("time1= ",time1)
    coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
   # println("coef[1]= ",coef[1])
    time2 = currentTime + minPosRoot(coef, Val(2))
  #  println("time2= ",time2)
#=   if  1.4387457897838291<=currentTime < 1.43874578978382928
    println("2nd inPosRoot= ",minPosRoot(coef, Val(2)))
   end =#
    nextTime[index] = time1 < time2 ? time1 : time2
 #=    if  1.4387457897838291<=currentTime < 1.43874578978382928
      println("nextTime[$index]= ",nextTime[index])
     end =#
    return nothing
end





 function computeNextEventTime(j::Int,ZCFunRes::Float64,oldsignValue,currentTime,  nextEventTime, quantum::Vector{Float64},printCounter::Vector{Int}) #later specify args

  #ZCFun=zcFunctions[j](q,d,tq,1) # 1 for order 1
  
  #ZCFun=ZCFun0(currentTime)
  if oldsignValue[j,1] != sign(ZCFunRes)
    nextEventTime[j]=currentTime 
    
     #=  println("old sign of",j,"  =  ",oldsignValue[j,1])
      println("new sign of",j,"  =  ",sign(ZCFun.coeffs[1])) =#
      #println("from computeNextEventTime:")
      if printCounter[1]>0 
      println("old value of",j,"  =  ",oldsignValue[j,2])
      println("new value of",j,"  =  ",(ZCFun.coeffs[1]))
      println("schedule event now at= ",nextEventTime[j])
      
      end
      
    
  else
   # nextEventTime[j] =currentTime + minPosRoot(ZCFun.coeffs, Val(2)) #Inf  # we can estimate the time. this is important when zc depends only on time
    nextEventTime[j]=Inf
    
  #=   println("debug")
    println("oldsignValue= ",oldsignValue )
    println("ZCFun= ", ZCFun) =#
   #=  if printCounter[1]>0 && printCounter[2] <20
      #println("from computeNextEventTime:")
      println("nextEventTime= ",nextEventTime)
      println("ZCFun2= ",ZCFun)
     
      printCounter[1]-=1
    end =#
  
  
  
    
      
  end
  oldsignValue[j,1]=sign(ZCFunRes)#update the values
  oldsignValue[j,2]=ZCFunRes
  
end
function computeNextEventTime(j::Int,ZCFun::Taylor0{Float64},oldsignValue,currentTime,  nextEventTime, quantum::Vector{Float64},printCounter::Vector{Int}) #later specify args
  if oldsignValue[j,1] != sign(ZCFun.coeffs[1])
    nextEventTime[j]=currentTime 
  else
    #nextEventTime[j] =currentTime + minPosRoot(ZCFun.coeffs, Val(2)) #Inf  # we can estimate the time. this is important when zc depends only on time   
    nextEventTime[j]=Inf
  end
  oldsignValue[j,1]=sign(ZCFun.coeffs[1])#update the values
  oldsignValue[j,2]=ZCFun.coeffs[1]
end


function fbarier(x::Taylor0{Float64},cacheT::Taylor0{Float64})
  x.coeffs[2] = cacheT.coeffs[1]
  x.coeffs[3] = (cacheT).coeffs[2] / 2
end