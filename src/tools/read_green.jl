# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Use the ZEN Framework
using ZenCore
using Printf

cfg = inp_toml("../SrVO3.toml", true)
rev_dict(cfg)
ai = GetImpurity()

fmesh, Delta = read_delta(ai, "dmft.green")
_, qdim, nmesh, nspin, nsite = size(Delta)

file = "dmft.green.diag"
open(file, "w") do fout
    for m = 1:nmesh
        @printf(fout, "%6i%16.8f", m, fmesh[m])
        for q = 1:qdim
            z = Delta[q,q,m,1,1]
            @printf(fout, "%16.8f%16.8f", real(z), imag(z))
        end
        println(fout)
    end
end
