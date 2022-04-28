@info "Loading packages and modules"
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
@info "Loading treasury data"
file_path_treas = LoadData.load_yield_curves(treas_urls, "UST", "data/raw/full_yields_treasury.csv")
file_path_corp = LoadData.load_yield_curves(corp_urls, "HQM", "data/raw/full_yields_corp.csv")
etf_path = "data/raw/full_etf_prices.csv"
if !isfile(etf_path)
    @info "Loading ETF data from AlphaVantage. Note: this is slow for free keys since downloads are limited to 5 requests per minute"
    LoadData.load_etfs(etf_path)
end
#### Constrained optimization approach
constr_eff_front_file_path = "data/model_results/constr_opt_eff_front.csv"
constr_eff_front_forecast_file_path = "data/model_results/constr_opt_eff_front_forecast.csv"
if !isfile(constr_eff_front_file_path)
    @info "Running optimization on historical data"
    constr_eff_front = ConstrOpt.generate_all_eff_frontiers(1000)
    CSV.write(constr_eff_front_file_path, constr_eff_front)
else
    constr_eff_front = CSV.read(constr_eff_front_file_path, DataFrame)
end
if !isfile(constr_eff_front_forecast_file_path)
    @info "Running optimization on forecast data"
    constr_eff_front_forecast = ConstrOpt.gen_forecast_eff_front(1000, file_path_treas, file_path_corp)
    CSV.write(constr_eff_front_forecast_file_path, constr_eff_front_forecast)
else
    constr_eff_front_forecast = CSV.read(constr_eff_front_forecast_file_path, DataFrame)
end
@info "Consolidating analysis"
CSV.write("data/model_results/constr_opt_eff_front_agg.csv", ConstrOpt.consolidate_frontier_scens(constr_eff_front))
CSV.write("data/model_results/constr_opt_eff_front_forecast_agg.csv", ConstrOpt.consolidate_frontier_scens(constr_eff_front_forecast))
