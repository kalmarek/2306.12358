using Dates
using Serialization
using Logging

import JuMP

function get_solution(model, wd, varP, eps = 1e-10)
    λ = JuMP.value(model[:λ])

    @info "reconstructing the solution"
    Q = @time let wd = wd, Ps = [JuMP.value.(P) for P in varP], eps = eps
        PropertyT.__droptol!.(Ps, 100eps)
        Qs = real.(sqrt.(Ps))
        PropertyT.__droptol!.(Qs, eps)
        PropertyT.reconstruct(Qs, wd)
    end

    solution = Dict(:λ => λ, :Q => Q)

    return solution
end

function run_optimization(log_file, model, optimizer, warm)
    @info "Current logfile is $log_file."
    isdir(dirname(log_file)) || mkpath(dirname(log_file))

    status, warm = @time PropertyT.solve(log_file, model, optimizer, warm)

    return status, warm
end

function solve_in_loop(model::JuMP.Model, args...; logdir, optimizer, data)
    @info "logging to $logdir"
    status = JuMP.UNKNOWN_RESULT_STATUS
    warm = try
        solution = deserialize(joinpath(logdir, "solution.sjl"))
        warm = solution[:warm]
        @info "trying to warm-start model with λ=$(solution[:λ])..."
        warm
    catch
        nothing
    end
    old_λ = 0.0
    certified = false
    while status != JuMP.OPTIMAL
        date_str = string(now())
        if Sys.iswindows()
            date_str = replace(date_str, ':' => '_')
        end
        log_file = joinpath(logdir, "solver_$date_str.log")
        status, warm = run_optimization(log_file, model, optimizer, warm)
        solution = get_solution(model, args...)
        solution[:warm] = warm
        serialize(joinpath(logdir, "solution_$date_str.sjl"), solution)
        serialize(joinpath(logdir, "solution.sjl"), solution)

        certified, λ = open(log_file; append = true) do io
            with_logger(SimpleLogger(io)) do
                return PropertyT.certify_solution(
                    data.elt,
                    data.unit,
                    solution[:λ],
                    solution[:Q];
                    halfradius = data.halfradius,
                )
            end
        end

        if certified == true
            @info "Certification done with λ = $λ" status
        end

        if status == JuMP.OPTIMAL
            return certified, λ
        else
            relative_change = abs(λ - old_λ) / (abs(λ) + abs(old_λ))
            @info "Relative improvement for λ" relative_change
            if relative_change < 1e-9
                @info "No progress detected, breaking" λ relative_change status
                return certified, λ
            end
        end
        old_λ = λ
    end

    return certified, old_λ
end
