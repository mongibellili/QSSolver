rhs=((schema.args[i].args[2]))
#println(typeof(rhs))
basi = convert(Basic, rhs)#:(2u1+3u2))
#extract jaco components
for j=1:T            # the number of diffequ always coincides with the number of continuous vars?
    coef=diff(basi,u[j])
    push!(jacArr,coef)
end