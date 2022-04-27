import CSV
using DataFrames
include("src/data/load_data.jl")
include("src/model/constrained_optimization.jl")
if !isfile("data/raw/full_etf_prices.csv") | !isfile("data/raw/full_yields.csv")
    LoadData.load_data()
end

#### Constrained optimization approach
constr_eff_front_file_path = "data/model_results/constr_opt_eff_front.csv"
if !isfile(constr_eff_front_file_path)
    @time constr_eff_front = ConstrOpt.generate_all_eff_frontiers(1000)
    CSV.write(constr_eff_front_file_path, constr_eff_front)
else
    constr_eff_front = CSV.read(constr_eff_front_file_path, DataFrame)
end
constr_eff_front_agg = ConstrOpt.consolidate_frontier_scens(constr_eff_front)
ConstrOpt.plot_eff_front(constr_eff_front_agg)
ConstrOpt.plot_allocations(constr_eff_front_agg)
ConstrOpt.see_allocations_at_risk_level(constr_eff_front_agg, 5)
using Plots
using GLM
using Statistics
using Dates
using MultivariateStats
using StatsPlots
using GaussianMixtures
using StatsBase
using Random
using Distributions
using ProgressBars
include("src/data/clean_data.jl")
# fix this to apply all preprocessing inside function
file_path_treas = "data/raw/full_yields_treasury.csv"
file_path_corp = "data/raw/full_yields_corp.csv"
yields_treas = CSV.read(file_path_treas, DataFrame)
yields_corp = CSV.read(file_path_corp, DataFrame)
if typeof(yields_treas.Date[1]) != Date
    yields_treas.Date = parse.(Date, yields_treas.Date)
end
if typeof(yields_corp.Date[1]) != Date
    yields_corp.Date = parse.(Date, yields_corp.Date)
end
sort!(yields_corp, order(:Date, rev=true))
sort!(yields_treas, order(:Date, rev=true))
(nrow(yields_corp) != nrow(yields_treas)) | (sum(yields_corp.Date .!= yields_treas.Date) > 0)
yields = hcat(yields_corp, yields_treas[:, 2:end])
cleaned_rets = CleanData.clean_data(1)
durations = Int.(1:30)
half_durations = durations .- 0.5
column_names = ["Date"]
append!(column_names, "HQMDur" .* string.(durations))
append!(column_names, "HQMDur" .* string.(half_durations))
append!(column_names, "USTDur" .* string.(durations))
append!(column_names, "USTDur" .* string.(half_durations))
yields = yields[1:nrow(cleaned_rets)+1, column_names]
yields_mat = Matrix(yields[:, Not(:Date)])
yields_cont = @. log((1 + yields_mat / 200)^2)
δ_yields = yields_cont[2:end, :] .- yields_cont[1:(end-1), :]
yields_pca = fit(PCA, δ_yields |> transpose, maxoutdim=7)
plot(yields_pca.proj[:, 1])
plot(yields_pca.proj[:, 2])
plot(yields_pca.proj[:, 3])
plot(yields_pca.proj[:, 4])
yields_components = MultivariateStats.transform(yields_pca, δ_yields |> transpose) |> transpose |> Array
density(yields_components)

start_durations = [parse(Float64, x[7:end]) for x in names(yields)[2:end]]
end_durations = start_durations .- 1 / 12
function generate_bond_returns(yields_cont, start_durations, end_durations)
    start_prices = hcat([exp.(.-yields_cont[:, i] .* start_durations[i]) for i in 1:length(start_durations)]...)
    end_prices = hcat([exp.(.-yields_cont[:, i] .* end_durations[i]) for i in 1:length(start_durations)]...)
    returns = @. log(end_prices[1:(end-1), :] / start_prices[2:end, :])
    return returns
end
returns = generate_bond_returns(yields_cont, start_durations, end_durations)
yields_returns = hcat(yields[2:end, [:Date]],
    DataFrame(returns, Symbol.(names(yields[:, Not(:Date)]))))
all_krd_indices = []
all_tickers = names(cleaned_rets[:, Not(:date)])
all_r2s = []
all_final_regressions = []
all_residuals = []
r2_threshold = 0.5
for i in 1:length(all_tickers)
    cand_returns = cleaned_rets[:, all_tickers[i]]
    # cand_returns = cleaned_rets[:, :VWEHX]
    krd_indices = []
    cand_cor = cor(cand_returns, returns) |> vec
    candidate_krd_index = argmax(abs.(cand_cor))
    regression = lm(returns[:, vcat(krd_indices, candidate_krd_index)], cand_returns)
    reg_confint = confint(regression, 0.999)
    insignificant = sign.(reg_confint[:, 1]) .!= sign.(reg_confint[:, 2])
    while (sum(insignificant) == 0) & (length(cand_cor) != length(krd_indices))
        push!(krd_indices, candidate_krd_index)
        cand_cor = cor(residuals(regression), returns) |> vec
        candidate_krd_index = argmax(abs.(cand_cor))
        regression = lm(returns[:, vcat(krd_indices, candidate_krd_index)], cand_returns)
        reg_confint = confint(regression, 0.999)
        insignificant = sign.(reg_confint[:, 1]) .!= sign.(reg_confint[:, 2])
    end
    if length(krd_indices) > 0
        regression = lm(returns[:, krd_indices], cand_returns)
        reg_r2 = r2(regression)
        if reg_r2 > r2_threshold
            push!(all_r2s, reg_r2)
            push!(all_final_regressions, lm(returns[:, krd_indices], cand_returns))
            push!(all_residuals, residuals(regression))
        else
            push!(all_r2s, 0.0)
            push!(all_final_regressions, nothing)
            push!(all_residuals, cand_returns)
        end
    else
        push!(all_r2s, 0.0)
        push!(all_final_regressions, nothing)
        push!(all_residuals, cand_returns)
    end
    push!(all_krd_indices, krd_indices)
end
lower = 0.5
upper = 1
all_tickers[lower.<all_r2s.<upper]
all_final_regressions[lower.<all_r2s.<upper]
all_krd_indices[lower.<all_r2s.<upper]
return_residuals = hcat(all_residuals...)
# yields_gmm = GMM(4, yields_components, kind = :full)
# em!(yields_gmm, yields_components, nIter = 10)
# yields_gmm.μ
# yields_gmm.Σ
# yields_gmm.w
# avll(yields_gmm, yields_components)
return_gmm = GMM(2, return_residuals, kind=:full, method=:split)
em!(return_gmm, return_residuals, nIter=10, varfloor=1e-12)
return_gmm.μ
return_gmm.w
covars(return_gmm)[1]
return_gmm.Σ[2] |> Matrix |> heatmap
state_weights = gmmposterior(return_gmm, return_residuals)[1]
states = [argmax(state_weights[i, :]) for i in 1:size(state_weights)[1]]
p₁₂ = sum((states[2:end] .== 1) .& (states[1:(end-1)] .== 2)) / sum((states[2:end] .== 1))
p₂₁ = sum((states[2:end] .== 2) .& (states[1:(end-1)] .== 1)) / sum((states[2:end] .== 2))
trans_mat = [[1 - p₁₂, p₁₂] [p₂₁, 1 - p₂₁]]
starting_state = argmax(state_weights[1, :])
yield_comp_cor = cor(yields_components, return_residuals)
heatmap(yield_comp_cor)
quantile(yield_comp_cor |> vec, float.((1:100) .- 0.5) / 100) |> plot
num_time_steps = 15 * 12
seed = 1
starting_yield = yields_cont[1, :]
start_date = yields.Date[1]
function forecast_fund_returns(num_time_steps, seed, starting_yield,
    yields_components, yields_pca, start_durations, end_durations,
    all_final_regressions, all_krd_indices,
    return_gmm, starting_date, all_tickers)

    rng = MersenneTwister(seed)
    forecasted_yield_pcs = randn(rng, num_time_steps + 1, size(yields_components)[2]) .* std(yields_components, dims=1)
    forecasted_yield_changes = reconstruct(yields_pca, forecasted_yield_pcs |> transpose) |> transpose |> Matrix
    forecasted_yields = cumsum(forecasted_yield_changes, dims=1) .+ (starting_yield |> transpose)
    forecasted_bond_returns = generate_bond_returns(forecasted_yields[1:end|>reverse, :], start_durations, end_durations)[1:end|>reverse, :]

    mean_fund_returns = []
    zero_vec = zeros(num_time_steps)
    for i in 1:length(all_final_regressions)
        if all_final_regressions[i] == nothing
            push!(mean_fund_returns, zero_vec)
        else
            push!(mean_fund_returns, predict(all_final_regressions[i], forecasted_bond_returns[:, all_krd_indices[i]]))
        end
    end
    mean_fund_returns = hcat(mean_fund_returns...)
    state_1_dist = MvNormal(means(return_gmm)[1, :], covars(return_gmm)[1])
    state_2_dist = MvNormal(means(return_gmm)[2, :], covars(return_gmm)[2])
    state_1_noise = rand(rng, state_1_dist, num_time_steps) |> transpose |> Matrix
    state_2_noise = rand(rng, state_2_dist, num_time_steps) |> transpose |> Matrix

    forecasted_states = Int[]
    trans_pulls = rand(rng, num_time_steps)
    current_state = starting_state
    for i in 1:num_time_steps
        if (current_state == 1) & (trans_pulls[i] < p₁₂)
            current_state = 2
        elseif (current_state == 2) & (trans_pulls[i] < p₂₁)
            current_state = 1
        end
        push!(forecasted_states, current_state)
    end
    final_noise = hcat([forecasted_states[i] == 1 ? state_1_noise[i, :] : state_2_noise[i, :]
                        for i in 1:num_time_steps]...) |> transpose |> Matrix
    final_returns = mean_fund_returns .+ final_noise
    final_returns = hcat(DataFrame(date=Dates.lastdayofmonth.(starting_date .+ Dates.Month.(1:num_time_steps))),
        DataFrame(final_returns, all_tickers))
    return final_returns
end

num_sims = 1000
all_eff_frontiers = Vector{DataFrame}(undef, num_sims)
iter = ProgressBar(1:num_sims)
Threads.@threads for i in iter
    all_rets = forecast_fund_returns(num_time_steps, i, starting_yield,
        yields_components, yields_pca, start_durations, end_durations,
        all_final_regressions, all_krd_indices,
        return_gmm, start_date, all_tickers)
    all_eff_frontiers[i] = ConstrOpt.generate_eff_frontier(i, all_rets, sample_returns=false)
end
all_eff_frontiers = vcat(all_eff_frontiers...)
constr_eff_front_agg = ConstrOpt.consolidate_frontier_scens(all_eff_frontiers)
ConstrOpt.plot_eff_front(constr_eff_front_agg)
ConstrOpt.plot_allocations(constr_eff_front_agg)
risk_level = 9
target_portfolio = ConstrOpt.see_allocations_at_risk_level(constr_eff_front_agg, risk_level)
allocation = Matrix(constr_eff_front_agg[risk_level, Not([:risk_level, :mean, :sd])] |> DataFrame) |> vec
naive_allocation = constr_eff_front_agg[risk_level, Not([:risk_level, :mean, :sd])] |> DataFrame
naive_allocation .= 0
naive_allocation[1, :VCIT] = 0.103
naive_allocation[1, :VGT] = 0.048
naive_allocation[1, :VHT] = 0.057
naive_allocation[1, :VWEHX] = 0.14
naive_allocation[1, :VWOB] = 0.272
naive_allocation[1, :VYMI] = 0.143
naive_allocation[1, :VTI] = 0.237
naive_allocation = Matrix(naive_allocation) |> vec
historical_returns = CleanData.clean_data(1)
port_return_hist = Matrix(historical_returns[:, Not(:date)]) * allocation
port_return_naive = Matrix(historical_returns[:, Not(:date)]) * naive_allocation
plot(exp.(cumsum(reverse(port_return_hist[(end-24):end]))), label="Historical Performance", legend=:topleft)
plot!(exp.(cumsum(reverse(port_return_naive[(end-24):end]))), label="Naive Performance")
all_simulated_returns = []
for i in 1:num_sims
    all_rets = forecast_fund_returns(num_time_steps, i, starting_yield,
        yields_components, yields_pca, start_durations, end_durations,
        all_final_regressions, all_krd_indices,
        return_gmm, start_date, all_tickers)
    portfolio_return = Matrix(all_rets[:, Not(:date)]) * allocation
    cum_port_ret = exp.(cumsum(portfolio_return))
    push!(all_simulated_returns, cum_port_ret)
end
sim_ret_mat = hcat(all_simulated_returns...)
avg_sim_ret = mean(sim_ret_mat, dims=2)
sd_sim_ret = std(sim_ret_mat, dims=2)
mean(sim_ret_mat[num_years_till_house*12, :] .< 1)
plot(avg_sim_ret, color=:black, label="Average Return")
plot!(sim_ret_mat, alpha=0.05, label=:none, color=:black)
hline!([1], label=:none, color=:red)
plot(avg_sim_ret, color=:black, label="Average Return")
plot!(avg_sim_ret .+ sd_sim_ret .* 1.65, color=:black, label="Upper 95% CI", alpha=0.5, linestyle=:dash)
plot!(avg_sim_ret .- sd_sim_ret .* 1.65, color=:black, label="Lower 95% CI", alpha=0.5, linestyle=:dash)
x
all_nonhit_chance = Vector{Float64}(undef, nrow(constr_eff_front_agg))
num_years_till_house = 7
for j in ProgressBar(1:nrow(constr_eff_front_agg))
    allocation = Matrix(constr_eff_front_agg[j, Not([:risk_level, :mean, :sd])] |> DataFrame) |> vec
    all_simulated_returns = []
    for i in 1:num_sims
        all_rets = forecast_fund_returns(num_time_steps, i, starting_yield,
            yields_components, yields_pca, start_durations, end_durations,
            all_final_regressions, all_krd_indices,
            return_gmm, start_date, all_tickers)
        portfolio_return = Matrix(all_rets[:, Not(:date)]) * allocation
        cum_port_ret = exp.(cumsum(portfolio_return))
        push!(all_simulated_returns, cum_port_ret)
    end
    sim_ret_mat = hcat(all_simulated_returns...)
    avg_sim_ret = mean(sim_ret_mat, dims=2)
    all_nonhit_chance[j] = mean(sim_ret_mat[num_years_till_house*12, :] .< 1)
end
plot(all_nonhit_chance)

