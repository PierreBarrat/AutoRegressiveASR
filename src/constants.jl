const AA_IQTREE_ALPHABET = "ARNDCQEGHILKMFPSTWYV"


using BioSequenceMappings
const alphabet_iqtree = Alphabet(AA_IQTREE_ALPHABET)

function parse_iqtree_model_file(model)
    lines = readlines(srcdir("iqtree_models.txt"))
    istart = findfirst(s -> occursin("model $(model)=", s), lines) + 1
    iend = findfirst(s -> occursin(";", s), lines[istart:end]) + istart - 1
    # [istart, iend-1] --> lower diag of rate matrix (not Q but Q/π)
    # [iend] --> equilibrium frequencies
    q = 20
    R = zeros(q, q)
    for (a, l) in enumerate(lines[istart:iend-1])
        R[a+1, 1:a] .= map(x -> parse(Float64, x), split(l, " "))
    end
    R .= R + R'

    p = map(x -> parse(Float64, x), split(lines[iend][1:end-1], " "))
    return R, p
end

