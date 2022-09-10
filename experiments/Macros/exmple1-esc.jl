    #--------------------------------------------------#
    #                   esc                            #
    #--------------------------------------------------#
#################example1#####################"
#= macro containervar(container,elemt)
    return esc(:($(Symbol(container,elemt))=$container[$elemt]))
end
word="the-second-letter-is-h"
elemt=2
display(@macroexpand @containervar word elemt)#:(wordelemt = word[elemt])
myexpr= @containervar word elemt
display(myexpr)#'h': ASCII/Unicode U+0068
 =#
#################example2#####################"
#since the  macro is run before runtime, when a var elemt inside, it will be used
#= macro containervar(container,elemt)
    elemt=3
    return esc(:($(Symbol(container,elemt))=$container[$elemt]))
end
letter=2
word="the-second-letter-is-h"
display(@macroexpand @containervar word letter)#:(word3 = word[3])
myexpr= @containervar word letter
display(myexpr)#'e': ASCII/Unicode U+0065 =#


#################example2#####################"
#= #removed interpolation $
macro containervar(container,elemt)
    return esc(:($(Symbol(container,elemt))=$container[elemt]))#when  elemt not interpolated it will look for it outside even if it exists inside, if not found then error
end
elemt=2  
word="the-second-letter-is-h"
display(@macroexpand @containervar word 2)#:(word2 = word[2])
myexpr= @containervar word 2
display(myexpr)#'h': ASCII/Unicode U+0068 =#


#################example2#####################"
#removed interpolation $
#= macro containervar(container,elemt)
    elemt=3
    return esc(:($(Symbol(container,elemt))=$container[$elemt]))#when  elemt not interpolated it will look for it outside, if not found then error
end
elemt=2
word="the-second-letter-is-h"
display(@macroexpand @containervar word 2)#:(word2 = word[2])
myexpr= @containervar word 2
display(myexpr)#'h': ASCII/Unicode U+0068 =#

    #--------------------------------------------------#
    #                  no esc                          #
    #--------------------------------------------------#

    #################example3#####################"
# removed esc: this is the most interesting example. interpolation has to be there. notice without esc it uses the inside var
macro containervar(container,elemt)
    elemt=3
    return (:($(Symbol(container,elemt))=$container[$elemt]))#when  elemt not interpolated it will look for it outside
end
elemt=2  
word="the-second-letter-is-h"
display(@macroexpand @containervar word 2)#:(word2 = word[2])
myexpr= @containervar word 2
display(myexpr)#'h': ASCII/Unicode U+0068

#################example3#####################"
#removed interpolation $ and removed esc
#= macro containervar(container,elemt)
    elemt=3
    return (:($(Symbol(container,elemt))=$container[elemt]))#when  elemt not interpolated it will look for it outside
end
elemt=2  
word="the-second-letter-is-h"
display(@macroexpand @containervar word 2)#:(word2 = word[2])
myexpr= @containervar word 2
display(myexpr)#'h': ASCII/Unicode U+0068 =#

#################example3#####################"
#no esc
#= macro containervar(container,elemt)
    container="localVar-secondletteris-o"
    return (:($(Symbol(container,elemt))=$container[$elemt]))
end
word="the-second-letter-is-h"
display(@macroexpand @containervar word 2)#:(var"#2#localVar-secondletteris-o2" = ("localVar-secondletteris-o")[2])
myexpr= @containervar word 2
display(myexpr)#'o': ASCII/Unicode U+006F =#
