include(joinpath(@__DIR__, "preamble.jl"))
include(joinpath(@__DIR__, "..", "src", "G₂_gens.jl"))
G, roots, Weyl = G₂_roots_weyl()
@info "Running Δ² - λ·Δ sum of squares decomposition for G₂"

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

Δ = RG(length(S)) - sum(RG(s) for s in S)
elt = Δ^2
unit = Δ

@info "defining the optimization problem"
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
    logdir = "./log/G2/r=$HALFRADIUS/Δ²-$(UPPER_BOUND)Δ",
    optimizer = scs_optimizer(;
        eps = 1e-10,
        max_iters = 50_000,
        accel = 50,
        alpha = 1.95,
    ),
    data = (elt = elt, unit = unit, halfradius = HALFRADIUS),
)

if certified && λ > 0
    Κ(λ, S) = round(sqrt(2λ / length(S)), Base.RoundDown; digits = 5)
    @info "Certified result: G₂ has property (T):" λ R = HALFRADIUS Κ(λ, S)
else
    @info "Could NOT certify the result:" certified λ
end
