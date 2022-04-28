import CSV
using DataFrames
include("src/data/load_data.jl")
include("src/model/constrained_optimization.jl")
treas_urls = ["https://home.treasury.gov/system/files/226/tnceom_03_07.xls",
    "https://home.treasury.gov/system/files/226/tnceom_08_12.xls",
    "https://home.treasury.gov/system/files/226/tnceom_13_17.xls",
    "https://home.treasury.gov/system/files/226/tnceom_18_22.xls"]
corp_urls = ["https://home.treasury.gov/system/files/226/hqmeom_04_08.xls",
    "https://home.treasury.gov/system/files/226/hqmeom_09_13.xls",
    "https://home.treasury.gov/system/files/226/hqmeom_14_18.xls",
    "https://home.treasury.gov/system/files/226/hqmeom_19_23.xls"]

file_path_treas = LoadData.load_yield_curves(treas_urls, "UST", "data/raw/full_yields_treasury.csv")
file_path_corp = LoadData.load_yield_curves(corp_urls, "HQM", "data/raw/full_yields_corp.csv")
if !isfile("data/raw/full_etf_prices.csv") | !isfile("data/raw/full_yields.csv")
    LoadData.load_etfs()
end
#### Constrained optimization approach
constr_eff_front_file_path = "data/model_results/constr_opt_eff_front.csv"
constr_eff_front_forecast_file_path = "data/model_results/constr_opt_eff_front_forecast.csv"
if !isfile(constr_eff_front_file_path)
    constr_eff_front = ConstrOpt.generate_all_eff_frontiers(1000)
    CSV.write(constr_eff_front_file_path, constr_eff_front)
else
    constr_eff_front = CSV.read(constr_eff_front_file_path, DataFrame)
end
if !isfile(constr_eff_front_forecast_file_path)
    constr_eff_front_forecast = ConstrOpt.gen_forecast_eff_front(1000, file_path_treas, file_path_corp)
    CSV.write(constr_eff_front_forecast_file_path, constr_eff_front_forecast)
else
    constr_eff_front_forecast = CSV.read(constr_eff_front_forecast_file_path, DataFrame)
end
constr_eff_front_forecast = ConstrOpt.gen_forecast_eff_front(10, file_path_treas, file_path_corp)
constr_eff_front_agg = ConstrOpt.consolidate_frontier_scens(constr_eff_front_forecast)
CSV.write(constr_eff_front_file_path, constr_eff_front_forecast)
ConstrOpt.plot_eff_front(constr_eff_front_agg)
ConstrOpt.plot_allocations(constr_eff_front_agg)
ConstrOpt.see_allocations_at_risk_level(constr_eff_front_agg, 5)
str_eff_front_agg = ConstrOpt.consolidate_frontier_scens(all_eff_frontiers)
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

