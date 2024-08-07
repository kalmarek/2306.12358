include(joinpath(@__DIR__, "preamble.jl"))
include(joinpath(@__DIR__, "..", "src", "G₂_gens.jl"))
G, roots, Weyl = G₂_roots_weyl()
@info "Running Adj - λ·Δ sum of squares decomposition for G₂"

@info "computing group algebra structure"
RG, S, sizes = @time PropertyT.group_algebra(G, halfradius = HALFRADIUS)

@info "computing WedderburnDecomposition"
wd = let Σ = Weyl, RG = RG
    act = PropertyT.AlphabetPermutation{eltype(Σ),Int64}(
        Dict(g => g for g in Σ),
    )

    @time wedderburn_decomposition(RG, Σ, act, 1:sizes[HALFRADIUS])
end
@info wd

function desubscriptify(symbol::Symbol)
    digits = [
        Int(l) - 0x2080 for
        l in reverse(string(symbol)) if 0 ≤ Int(l) - 0x2080 ≤ 9
    ]
    res = 0
    for (i, d) in enumerate(digits)
        res += 10^(i - 1) * d
    end
    return res
end

function PropertyT.grading(g::MatrixGroups.MatrixElt, roots = roots)
    id = desubscriptify(g.id)
    return roots[id]
end

Δ = RG(length(S)) - sum(RG(s) for s in S)
Δs = PropertyT.laplacians(
    RG,
    S,
    x -> (gx = PropertyT.grading(x); Set([gx, -gx])),
)

elt = PropertyT.Adj(Δs)
@assert elt == Δ^2 - PropertyT.Sq(Δs)
unit = Δ

@time model, varP = PropertyT.sos_problem_primal(
    elt,
    unit,
    wd;
    upper_bound = UPPER_BOUND,
    augmented = true,
    show_progress = isinteractive(),
)

certified, λ = solve_in_loop(
    model,
    wd,
    varP;
    logdir = "./log/G2/r=$HALFRADIUS/Adj-$(UPPER_BOUND)Δ",
    optimizer = scs_optimizer(;
        eps = 1e-9,
        max_iters = 100_000,
        accel = 50,
        alpha = 1.95,
    ),
    data = (elt = elt, unit = unit, halfradius = HALFRADIUS),
)

if certified && λ > 0
    @info "Certified result: in group G₂ we have Adj -λΔ ≥ 0" λ R = HALFRADIUS
else
    @info "Could NOT certify the result:" certified λ
end
