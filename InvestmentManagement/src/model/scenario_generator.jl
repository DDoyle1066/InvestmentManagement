module ScenGen
using GLM
using Statistics
using Dates
using MultivariateStats
using StatsPlots
using GaussianMixtures
using StatsBase
using Random
using Distributions
using DataFrames
using Suppressor
import CSV
function get_relevant_yields(file_path_treas, file_path_corp, durations, cleaned_rets)
    # fix this to apply all preprocessing inside function
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
    min_date = max(yields_corp.Date[end], yields_treas.Date[end])
    yields_corp = yields_corp[yields_corp.Date.>=min_date, :]
    yields_treas = yields_treas[yields_treas.Date.>=min_date, :]
    @assert nrow(yields_corp) == nrow(yields_treas)
    @assert sum(yields_corp.Date .!= yields_treas.Date) == 0
    yields = hcat(yields_corp, yields_treas[:, 2:end])
    column_names = ["Date"]
    append!(column_names, "HQMDur" .* string.(durations))
    append!(column_names, "USTDur" .* string.(durations))
    yields = yields[1:nrow(cleaned_rets)+1, column_names]
    return yields
end
function generate_bond_returns(yields_cont, start_durations, end_durations)
    start_prices = hcat([exp.(.-yields_cont[:, i] .* start_durations[i]) for i in 1:length(start_durations)]...)
    end_prices = hcat([exp.(.-yields_cont[:, i] .* end_durations[i]) for i in 1:length(start_durations)]...)
    returns = @. log(end_prices[1:(end-1), :] / start_prices[2:end, :])
    return returns
end
function get_yields_info(yields)
    yields_mat = Matrix(yields[:, Not(:Date)])
    yields_cont = @. log((1 + yields_mat / 200)^2)
    δ_yields = yields_cont[2:end, :] .- yields_cont[1:(end-1), :]
    yields_pca = fit(PCA, δ_yields |> transpose, maxoutdim=7)
    yields_components = MultivariateStats.transform(yields_pca, δ_yields |> transpose) |> transpose |> Array
    start_durations = [parse(Float64, x[7:end]) for x in names(yields)[2:end]]
    end_durations = start_durations .- 1 / 12
    bond_returns = generate_bond_returns(yields_cont, start_durations, end_durations)
    return yields_pca, yields_cont, yields_components, bond_returns, start_durations, end_durations
end

function get_etf_generator(cleaned_rets, bond_returns)
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
        cand_cor = cor(cand_returns, bond_returns) |> vec
        candidate_krd_index = argmax(abs.(cand_cor))
        regression = lm(bond_returns[:, vcat(krd_indices, candidate_krd_index)], cand_returns)
        reg_confint = confint(regression, 0.999)
        insignificant = sign.(reg_confint[:, 1]) .!= sign.(reg_confint[:, 2])
        while (sum(insignificant) == 0) & (length(cand_cor) != length(krd_indices))
            push!(krd_indices, candidate_krd_index)
            cand_cor = cor(residuals(regression), bond_returns) |> vec
            candidate_krd_index = argmax(abs.(cand_cor))
            regression = lm(bond_returns[:, vcat(krd_indices, candidate_krd_index)], cand_returns)
            reg_confint = confint(regression, 0.999)
            insignificant = sign.(reg_confint[:, 1]) .!= sign.(reg_confint[:, 2])
        end
        if length(krd_indices) > 0
            regression = lm(bond_returns[:, krd_indices], cand_returns)
            reg_r2 = r2(regression)
            if reg_r2 > r2_threshold
                push!(all_r2s, reg_r2)
                push!(all_final_regressions, lm(bond_returns[:, krd_indices], cand_returns))
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
    return_residuals = hcat(all_residuals...)
    return_gmm = GMM(2, return_residuals, kind=:full, method=:split)
    em!(return_gmm, return_residuals, nIter=10, varfloor=1e-12)
    state_weights = gmmposterior(return_gmm, return_residuals)[1]
    states = [argmax(state_weights[i, :]) for i in 1:size(state_weights)[1]]
    p₁₂ = sum((states[2:end] .== 1) .& (states[1:(end-1)] .== 2)) / sum((states[2:end] .== 1))
    p₂₁ = sum((states[2:end] .== 2) .& (states[1:(end-1)] .== 1)) / sum((states[2:end] .== 2))
    starting_state = argmax(state_weights[1, :])
    return all_final_regressions, all_krd_indices, p₂₁, p₁₂, return_gmm, all_tickers, starting_state
end

function forecast_fund_returns(num_time_steps, seed, starting_yield,
    yields_components, yields_pca, start_durations, end_durations,
    all_final_regressions, all_krd_indices, p₂₁, p₁₂, starting_state,
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
end

