using StaticArrays
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
        println("simple= ",mpr)
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
        println("not simple= ",mpr)
    end
    return mpr
end



coef=@SVector[0.25,1.0 , 0.75]#[c,b,a] coef of quadratic eq
  time1 =  minPosRoot(coef, Val(2))
  println("time1",time1)





#= q=1.0091287092917527
derq=1.9452277442494834
x=1.0091287092917527
derx=1.9452277442494834
derderx=-3.0
quantum=0.0010091287092917527
currentTime=0.03659798118203518

coef=@SVector[q - x - quantum,derq-derx , -derderx]#[c,b,a] coef of quadratic eq
  time1 = currentTime + minPosRoot(coef, Val(2))
  println("time1",time1)
  coef=setindex(coef,q - x + quantum,1)
  time2 = currentTime + minPosRoot(coef, Val(2))
  println("time2",time2)
  nextTime = time1 < time2 ? time1 : time2 

display(nextTime) =#