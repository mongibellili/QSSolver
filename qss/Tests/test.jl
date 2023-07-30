ex=:((2,9),i)
dump(ex)
  #= if i==1
        exacteA(qv,cacheA,1,2);a12=cacheA[1]
        utemp=x1-a*qaux[i][1]-a12*qv[2][0] #dx1=a11 q1 + a12 q2 +utemp
        utem2=x2-a*q1-a12*qv[2][1]
        if abs(utemp-1600.0)>1e-2
            @show utemp,i
            @show a,a12
        end
        if abs(utem2-0.0)>1e-2
            @show utem2,i ############this might appear because q1 updated and x2 did not
            @show a,a12
        end
   else
        exacteA(qv,cacheA,2,1);a21=cacheA[1]
        utemp=x1-a*q-a21*qv[1][0]#dx2=a22 q2 + a21 q1 +utemp
        utem2=x2-a*q1-a21*qv[1][1]
        if abs(utemp-0.2)>1e-2
            @show utemp,i
            @show a,a21
        end
        if abs(utem2-0.0)>1e-2
            @show utem2,i
            @show a,a21
        end
    end =#

  # uv[i][i][2]= u2
   # tu[i]=simt  
    # olddx[i][2]=2*x2# 