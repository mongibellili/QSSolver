



function minPosRoot(coeff::SVector{2,Float64}, ::Val{1}) # coming from val(1) means coef has x and derx only...size 2
    mpr=-1
        if coeff[2] == 0 
            mpr = Inf
        else 
            mpr = -coeff[1] / coeff[2];
        end
        if mpr < 0
            mpr = Inf
        end
       # println("mpr inside minPosRoot in utils= ",mpr)
    return mpr
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
@inline function minPosRoot(coeffs::NTuple{3,Float64}, ::Val{2})
	_a = 1.0 / coeffs[1]
	b, c = -0.5coeffs[2] * _a, coeffs[3] * _a
	Δ = muladd(b, b, -c) # b * b - c
	if Δ < -4eps() # Complex roots
		Inf
	elseif Δ > 4eps() # Real roots
		if b > eps()
			c > eps() ? c / (b + sqrt(Δ)) : b + sqrt(Δ)
		elseif b < -eps()
			c < eps() ? c / (b - sqrt(Δ)) : Inf
		else
			sqrt(-c)
		end
	else # Double real root
		b > -eps() ? b : Inf
	end
end
function minPosRoot(coeff::Vector{Float64}, ::Val{1}) # coming from val(1) means coef has x and derx only...size 2
    mpr=-1
        if coeff[2] == 0 
            mpr = Inf
        else 
            mpr = -coeff[1] / coeff[2];
        end
        if mpr < 0
            mpr = Inf
        end
       # println("mpr inside minPosRoot in utils= ",mpr)
    return mpr
end
function minPosRoot(coeff::Vector{Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
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

#= function minPosRoot(coeff::SVector{N,Float64}, order ::Int) where{N}
    mpr=-1
    if order==1
        if coeff[2] == 0 
            mpr = Inf
        else 
            mpr = -coeff[1] / coeff[2];
            #println(mpr)
        end
        if mpr < 0
            mpr = Inf
        end

    else
        println("logicof order2")
    end
    return mpr
end =#



#=
function integrateState(j::Int,order::Int,elapsed::Float64,x::  Vector{Float64} )
   
    if order ==1
        x[2*j-1]=x[2*j-1]+ elapsed * (x[2*j])
       
       # push!(x[2*j-1],(x[2*j-1])+ elapsed * (x[2*j]))
    else
        println("keep coding next orders")
    end
    
end

function minPosRoot(coeff::Vector{Float64}, order ::Int)
    mpr=-1
    if order==1
        if coeff[2] == 0 
            mpr = Inf
        else 
            mpr = -coeff[1] / coeff[2];
            #println(mpr)
        end
        if mpr < 0
            mpr = Inf
        end

    else
        println("logicof order2")
    end
    return mpr
end

function plotX(simulator :: QSS_simulator)
    x=simulator.qssData.x
    t=simulator.qssTime.tx
    states=simulator.settings.states
    order=simulator.settings.order

   for i = 1:states
      # println("steps to solve eq$i= ",length(x[(order+1)*i-1]))
        #println(length(t[i]))
        display(plot!(t[i],x[(order+1)*i-order]))#,xlims = (0.12,0.15),ylims = (1.23,1.26)))
        #, xlims = (0,0.00005),ylims = (1,1.0001)))
        #display(plot!(t[i],x[(order+1)*i])) #can't plot derivative cuz length does not equal time Vector length
    end
    
    readline()
end
=#





#=function plotX(simulator :: QSS_simulator)
    x=simulator.qssData.x
    t=simulator.qssTime.tx
    states=simulator.settings.states
    order=simulator.settings.order

   for i = 1:states
      # println("steps to solve eq$i= ",length(x[(order+1)*i-1]))
        #println(length(t[i]))
        display(plot!(t[i],x[(order+1)*i-order]))#,xlims = (0.12,0.15),ylims = (1.23,1.26)))
        #, xlims = (0,0.00005),ylims = (1,1.0001)))
        #display(plot!(t[i],x[(order+1)*i])) #can't plot derivative cuz length does not equal time Vector length
    end
    
    readline()
end =#