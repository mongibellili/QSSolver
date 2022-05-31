 using BenchmarkTools
using StaticArrays
using SymEngine
using qss

#= function normalQSS() 

   psettings = ProblemSettings(5.0,saveat(0.1),qss1())
   myProb=@odeProblem(begin
        u=[1.0,2.0,0.5] 
        d=[1.0,0.5]
        du1=u2+2.0     # du1...2...N are expected to be in order....later can fix this
        du2=-u1-u2
        du3=u3+u2+d1
        if u1+0.7 >0 
            d1=0.1
        else
            d1=1.0
        end
        if u2+d1 >0 
            d2=0.33
            u3=2.2
        else
            d2=1.0
        end
    end)  
   sol=QSS_Solve(psettings,myProb)   
end


#@btime normalQSS() 
normalQSS() 
 =#
 function normalQSS() 

    psettings = ProblemSettings(5.0,saveat(0.1),qss1())
    myProb=@odeProblem(begin
         u=[1.0,2.0,0.5] 
         d=[1.0,0.5]
         du1=u3-u2+2.0     # du1...2...N are expected to be in order....later can fix this
         du2=-u1-3u2-d1+t
         du3=u3+u2+d1+d2
         if u1+0.7+d2 >0   #still have not added 'usercode checking'....throw error msg
             d1=0.1
         else
             d1=1.0
             u2=4.5
         end
         if u2+d1 >0 
             d2=0.33
             u3=2.2
         else
             d2=1.0
         end
     end)  
    sol=QSS_Solve(psettings,myProb)   
 end
 
 
 #@btime normalQSS() 
 normalQSS() 
 