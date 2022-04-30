function check(str,arg;type=DataType,max=nothing,min=nothing,description="")
    @argcheck typeof(arg)==type
    @argcheck arg>min
    @argcheck arg<max
    @argcheck typeof(description)==String
   return arg
end


 function constr(name,arg,field)
        return :(function $name($arg,$field)
            new(check($name,$arg,$field))
        end)
    end
    
macro creatStruct(name,arg)
   code = Base.remove_linenums!(quote
     struct $name
        end
     end)
     print(arg)
   append!(code.args[1].args[3].args,[constr(name,arg.args[1].args[1],arg.args[1].args[2])])
   code 
end 

macro myStruct(name,arg)
    :(@creatStruct $name $arg)
end 





@myStruct test12 (
  (arg1,(max=10))
)