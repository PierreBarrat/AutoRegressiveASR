function project_path(pth::AbstractString, project::AbstractString = projectdir())
    return relpath(abspath(pth), project)
end

function pluto_ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
    m
end

function load_generative_model(file::AbstractString)
    ext = splitext(file)[2]
    return if ext == ".jld2"
        @info "Assume $(basename(file)) is an `ArNet`"
        JLD2.load(file)["arnet"]
    else
        @info "Assume $(basename(file)) is a `DCAGraph`"
        DCAGraph(file)
    end
end

# loglikelihood(sequence, model::DCAGraph) = -DCAToosl.energy(model, sequence)
# loglikelihood(sequence, model::ArNet) = ArDCA.loglikelihood(sequence, model)
# local field
function ar_local_field(s, i, arnet)
    j = findfirst(==(i), arnet.idxperm) - 1
    if j == 0
        return arnet.p0
    end

    H = copy(arnet.H[j])
    for k in 1:(j-1)
        H += arnet.J[j][:, s[arnet.idxperm[k]], k]
    end
    return ArDCA.softmax(H)
end
