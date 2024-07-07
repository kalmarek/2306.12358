include(joinpath(@__DIR__, "preamble.jl"))

G = MatrixGroups.SpecialLinearGroup{N}(Int8)
@info "Running Adj - λ·Δ sum of squares decomposition for " G

@info "computing group algebra structure"
RG, S, sizes = @time PropertyT.group_algebra(G, halfradius = HALFRADIUS)

@info "computing WedderburnDecomposition"
wd = let RG = RG, N = N
    G = StarAlgebras.object(RG)
    P = PermGroup(perm"(1,2)", Perm(circshift(1:N, -1)))
    Σ = Groups.Constructions.WreathProduct(PermGroup(perm"(1,2)"), P)
    act = PropertyT.action_by_conjugation(G, Σ)

    @time wedderburn_decomposition(RG, Σ, act, 1:sizes[HALFRADIUS])
end
@info wd

Δ = RG(length(S)) - sum(RG(s) for s in S)
Δs = PropertyT.laplacians(
    RG,
    S,
    x -> (gx = PropertyT.grading(x); Set([gx, -gx])),
)

elt = PropertyT.Adj(Δs, :A₂)
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
    logdir = "./log/SL($N,Z)/r=$HALFRADIUS/Δ²-$(UPPER_BOUND)Δ",
    optimizer = scs_optimizer(;
        eps = 1e-10,
        max_iters = 50_000,
        accel = 50,
        alpha = 1.95,
    ),
    data = (elt = elt, unit = unit, halfradius = HALFRADIUS),
)

if certified && λ > 0
    @info "Certified result: in group $G we have Adj_A₂ -λΔ ≥ 0" λ R =
        HALFRADIUS
else
    @info "Could NOT certify the result:" certified λ
end
