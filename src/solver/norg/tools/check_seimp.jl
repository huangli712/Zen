
using DelimitedFiles

# println("The out put is in: "*pwd()*"/out")# for customers.

# mkdir("out")
b = readdir()
for a in b
    m = readdlm(a)
    for i=28:52
        if(m[2, i]>0) println(a,"  ",i,"  ",m[1, i],"  ",m[2, i])
        end
    end
end