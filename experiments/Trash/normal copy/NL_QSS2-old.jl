function computeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}}, quantum::Vector{Float64})where{T}
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
  #function computeNextInputTime(::Val{2}, i::Int, currentTime::Float64,elapsed::Float64, t::Taylor1{Float64} ,f::Function,nextInputTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}}, quantum::MVector{T,Float64})where{T}
    function computeNextInputTime(::Val{2}, i::Int, currentTime::Float64,elapsed::Float64, t::Taylor1{Float64} ,f::Function,nextInputTime::Vector{Float64}, x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}}, quantum::Vector{Float64})where{T<:Int64}
    df=0.0
    oldDerDerX=((x[i].coeffs[3])*2.0)
    newDerDerX=((differentiate(f(q,t),1)).coeffs[1] )# do not put /factorial(2) cuz here we actually need derder not store the coeff
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
  function computeDerivative( ::Val{2} ,j::Int, f::Function ,x::Vector{Taylor1{Float64}} , q::Vector{Taylor1{Float64}}, t::Taylor1{Float64}  )
      #x[j].coeffs[2] =((differentiate(f(q,t),0)).coeffs[1] )/factorial(1) 
      x[j].coeffs[2] =(f(q,t)).coeffs[1]
      x[j].coeffs[3] =((differentiate(f(q,t),1)).coeffs[1] )/2#factorial(2)
      #println("x[j].coeffs[2]= ",x[j].coeffs[2])
      return nothing
  end
#=   function computeDerivative( ::Val{2} ,j::Int ,x::Vector{Taylor1{Float64}} , tt::Taylor1{Float64}  )
    #x[j].coeffs[2] =((differentiate(f(q,t),0)).coeffs[1] )/factorial(1) 
    x[j].coeffs[2] =(tt).coeffs[1]
    x[j].coeffs[3] =((differentiate(tt,1)).coeffs[1] )/2#factorial(2)
    #println("x[j].coeffs[2]= ",x[j].coeffs[2])
    return nothing
end =#
  function reComputeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}}, quantum::Vector{Float64})where{T}
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
  #= function minPosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
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
 =#