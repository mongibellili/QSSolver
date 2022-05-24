  
macro creatStruct(name,arg) 
  v=[]
  code11=:($name)
  push!(v,code11)
  code2=(quote
         # struct $name
          #$(:($arg.args[1]::$arg.args[2]))
          $(arg.args[1])::$(arg.args[2])
          #end # this should go at the end
        end)
   #v=[code11,code2]
   push!(v,code2)
   esc(Expr(:block,v...))    
end 
println(@macroexpand @creatStruct test12 (arg1,Int))