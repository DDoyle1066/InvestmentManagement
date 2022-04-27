module ConstrOpt
using Random
using Statistics
using StatsBase
using Convex, SCS
using ProgressBars
using DataFrames
using Plots
include("../data/clean_data.jl")
function generate_eff_frontier(scen_id, all_rets; sample_returns=true)
    if sample_returns
        sampled_indices = sample(MersenneTwister(scen_id), 1:nrow(all_rets), nrow(all_rets))
        sampled_returns = Matrix(all_rets[sampled_indices, Not(:date)])
    else
        sampled_returns = Matrix(all_rets[:, Not(:date)])
    end
    num_candidates = 20
    tickers = names(all_rets[:, Not(:date)])
    all_μ = Vector{Float64}(undef, num_candidates + 2)
    all_σ = Vector{Float64}(undef, num_candidates + 2)
    all_allocations = Matrix{Float64}(undef, num_candidates + 2, length(tickers))
    μ = mean(sampled_returns, dims=1)
    Σ = cov(sampled_returns)
    allocation = Variable(length(tickers))
    add_constraint!(allocation, allocation >= 0)
    add_constraint!(allocation, sum(allocation) <= 1)
    μ_hat = dot(allocation, μ)
    Σ_hat = quadform(allocation, Σ)
    # port_returns = sampled_returns * allocation

    problem_mean = maximize(μ_hat)
    solve!(problem_mean, SCS.Optimizer, silent_solver=true)
    all_μ[end] = evaluate(μ_hat)
    all_σ[end] = sqrt(evaluate(Σ_hat))
    all_allocations[end, :] = evaluate(allocation) |> vec
    problem_sd = minimize(Σ_hat)
    solve!(problem_sd, SCS.Optimizer, silent_solver=true)
    all_μ[1] = evaluate(μ_hat)
    all_σ[1] = sqrt(evaluate(Σ_hat))
    all_allocations[1, :] = evaluate(allocation) |> vec
    candidate_weights = float.(1:num_candidates) ./ num_candidates .- 0.5 / num_candidates
    candidate_Σs = @. (candidate_weights * all_σ[end] + (1 - candidate_weights) * all_σ[1])^2
    for i in 1:num_candidates
        problem = maximize(μ_hat, Σ_hat <= candidate_Σs[i])
        # problem_mean.constraints = [Σ_hat <= candidate_Σs[i]]
        solve!(problem, SCS.Optimizer, silent_solver=true, warmstart=true)
        all_μ[i+1] = evaluate(μ_hat)
        all_σ[i+1] = sqrt(evaluate(Σ_hat))
        all_allocations[i+1, :] = evaluate(allocation) |> vec
    end

    return_df = hcat(DataFrame(mean=all_μ .* 12, sd=all_σ .* sqrt(12), scen_id=scen_id, risk_level=1:(num_candidates+2)),
        DataFrame(all_allocations, tickers))
    return return_df
end

function generate_all_eff_frontiers(num_sims)
    all_eff_frontiers = Vector{DataFrame}(undef, num_sims)
    iter = ProgressBar(1:num_sims)
    Threads.@threads for i in iter
        all_rets = CleanData.clean_data(i)
        all_eff_frontiers[i] = generate_eff_frontier(i, all_rets)
    end
    return vcat(all_eff_frontiers...)
end
function consolidate_frontier_scens(constr_eff_front)
    return combine(groupby(constr_eff_front, :risk_level),
        names(constr_eff_front) .=> mean .=> names(constr_eff_front)
    )[:, Not(:scen_id)]
end
# consol_eff_frontier = combine(groupby(consol_eff_frontier, :risk_level),
#     names(consol_eff_frontier) .=> mean .=> names(consol_eff_frontier)
# )[:, Not(:scen_id)]
function plot_eff_front(constr_eff_front_agg)
    return plot(constr_eff_front_agg.sd, constr_eff_front_agg.mean)
end
function plot_allocations(constr_eff_front_agg)
    return heatmap(Matrix(constr_eff_front_agg[:, Not([:mean, :sd, :risk_level])]))
end
function see_allocations_at_risk_level(constr_eff_front_agg, risk_ind)
    risk_ret = filter(:risk_level => ==(risk_ind), constr_eff_front_agg)[:, [:mean, :sd]]
    risk_df = filter(:risk_level => ==(risk_ind), constr_eff_front_agg)[:, Not([:risk_level, :mean, :sd])]
    risk_df = stack(risk_df, variable_name=:ticker, value_name=:allocation)
    sort!(risk_df, order(:allocation, rev=true))
    risk_df.cum_alloc = cumsum(risk_df.allocation)
    risk_df = risk_df[risk_df.allocation.>=0.01, :]
    transform!(risk_df, :allocation => (x -> round.(x / sum(x), digits=4)) => :allocation_adjusted)
    return risk_df, risk_ret
end
# heatmap(Matrix(consol_eff_frontier[:, Not([:mean, :sd, :risk_level])]))
# # function view_at_risk_level 
# risk_ind = 22
# risk_df = filter(:risk_level => ==(risk_ind), consol_eff_frontier)[:, Not([:risk_level, :mean, :sd])]
# risk_df = stack(risk_df, variable_name=:ticker, value_name=:allocation)
# sort!(risk_df, order(:allocation, rev=true))
# risk_df.cum_alloc = cumsum(risk_df.allocation)
# risk_df = risk_df[risk_df.allocation.>=0.01, :]
# transform!(risk_df, :allocation => (x -> round.(x / sum(x), digits=4)) => :allocation_adjusted)
end