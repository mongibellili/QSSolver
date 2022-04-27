macro make_struct(struct_name, schema...)
    fields=[:($(entry.args[1])::$(entry.args[2])) for entry in schema]
    esc(quote struct $struct_name
        $(fields...)
        end
    end)
end
#println(@macroexpand @make_struct STRUCT_NAME (x,Int) (y, Float64))
@make_struct sim (x, Int) (y, Float64) (wRegister::MVector{2,Float64}(undef)
p=sim(1,2.0)
display(p.x)