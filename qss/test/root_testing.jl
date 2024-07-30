
using Plots
function quadRootv2(coeff::NTuple{3,Float64}) # call from mliqss1/simultupdate1
    mpr = (-1.0, -1.0) #size 2 to use mpr[2] in quantizer
    a = coeff[1]
    b = coeff[2]
    c = coeff[3]
    if a == 0.0 || ((1e7 * abs(a)) < abs(b) && a*b>0)# coef3 is the coef of t^2
      if b != 0.0
        if 0 < -c / b < 1e5  # neglecting a small 'a' then having a large h (in this case h=1e7 is large)would cause an error  'a*h*h' becomes large
          mpr = (-c / b, -1.0)
        end
      end
    elseif b == 0.0
      if -c / a > 0
        mpr = (sqrt(-c / a), -1.0)
      end
    elseif c == 0.0
      mpr = (-1.0, -b / a)
    else
      #double disc;
      Δ = 1.0 - 4.0 * c * a / (b * b)
      if Δ > 0.0
        sq = sqrt(Δ)
        r1 = -0.5 * (1.0 + sq) * b / a
        r2 = -0.5 * (1.0 - sq) * b / a
        mpr = (r1, r2)
      elseif Δ == 0.0
        r1 = -0.5 * b / a
        mpr = (r1, r1 - 1e-12)
      end
    end
    return mpr
  end

  function plotzc(a,b,c)
       zc(x)=a*x^2+b*x+c
       p1=plot(zc)
       savefig(p1, "plotzc.png")
  end
  function plotf(βi, ci, αi, bi )
    f(x)=(ci*x^2+bi*x)/(βi*x^2+αi*x+1)
    p1=plot(f, xlims=(-1e12, 1e17),ylims=(-1e1, 1e1))
    savefig(p1, "plotf7.png")
end
  function evzc(x,a,b,c)
    return a*x^2+b*x+c
    
  end
  function evf(x,βi, ci, αi, bi )
    return (ci*x^2+bi*x)/(βi*x^2+αi*x+1)
   
end
  coefi = (0.4788446667371318, -4.734219817347005e11, -1.0e-5)
 βi, ci, αi, bi =  0.0, 0.4788446667371318, 4.734219817110833e16, -23.6171875
 #plotf(βi, ci, αi, bi )

  #@show evzc(9.886754820944719e11,coefi...)
  @show evf(9.886754820944719e11,βi, ci, αi, bi)
  #@show coefi[3]/coefi[2]
  #  @show quadRootv2(coefi)
  #plotzc(coefi...)