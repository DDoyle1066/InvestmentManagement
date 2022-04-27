module CleanData
using DataFrames
using Statistics
using Random
using Dates
using GLM
import CSV

function generate_returns(ticker_df)
    ticker_df = ticker_df[2:end, :]
    prices = ticker_df[:, "adjusted close"]
    returns = @. log(prices[1:(end-1)] / prices[2:end])
    return DataFrame(ticker=ticker_df.ticker[1], returns=returns, date=ticker_df.timestamp[2:end])
end
function fill_in_least_missing!(ret_etf, column_offset, scenario_offset)
    num_missing_by_col = sum(Matrix(ismissing.(ret_etf[:, Not(:date)])), dims=1) |> vec
    lowest_value = minimum(num_missing_by_col[num_missing_by_col.>0])
    column_index = (1:length(num_missing_by_col))[num_missing_by_col.==lowest_value][1] + 1
    fill_column_indices = (1:length(num_missing_by_col))[num_missing_by_col.==0] .+ 1
    returns_to_fill = ret_etf[:, column_index]
    returns_to_fill_with = Float64.(Matrix(ret_etf[:, fill_column_indices]))
    cors = cor(returns_to_fill[1:(end-lowest_value)], returns_to_fill_with[1:(end-lowest_value), :]) |> vec
    cor_thresh = quantile(abs.(cors), 1 - 5 / length(cors))
    regress_indices = (1:length(cors))[abs.(cors).>cor_thresh]
    regression = lm(returns_to_fill_with[1:(end-lowest_value), regress_indices], Float64.(returns_to_fill[1:(end-lowest_value)]))
    reg_confint = confint(regression)
    insignificant = sign.(reg_confint[:, 1]) .!= sign.(reg_confint[:, 2])
    while sum(insignificant) > 0 & length(regress_indices) > 1
        ranges = reg_confint[:, 2] .- reg_confint[:, 1]
        worst_range = insignificant .& (ranges .== maximum(ranges[insignificant]))
        regress_indices = regress_indices[.!worst_range]
        regression = lm(returns_to_fill_with[1:(end-lowest_value), regress_indices], Float64.(returns_to_fill[1:(end-lowest_value)]))
        reg_confint = confint(regression)
        insignificant = sign.(reg_confint[:, 1]) .!= sign.(reg_confint[:, 2])
    end
    random_levels = rand(MersenneTwister(20220324 + column_offset + scenario_offset * 10000), lowest_value) .* 2 .- 1
    return_fills = [predict(regression, returns_to_fill_with[[end - lowest_value + i], regress_indices], interval=:prediction, level=random_levels[i])[:upper][1] for i in 1:lowest_value]
    ret_etf[(end-lowest_value+1):end, column_index] .= return_fills
    return ret_etf
end
function clean_data(scenario_offset)
    etfs = CSV.read("data/raw/full_etf_prices.csv", DataFrame)
    if !isa(etfs.timestamp[1], Date)
        etfs.timestamp = Date.(etfs.timestamp, DateFormat("m/d/y"))
    end
    transform!(groupby(etfs, :ticker), :volume => length => :num_obs)
    num_obs_by_ticker = combine(groupby(etfs, :ticker), :volume => length => :num_obs)
    months_to_include = median(num_obs_by_ticker.num_obs)
    # remove any tickers with less than 5 years of observations
    filter!(:num_obs => >=(60), etfs)
    sort!(etfs, [:ticker, order(:timestamp, rev=true)])
    ret_etf = combine(groupby(etfs, :ticker), generate_returns)
    ret_etf = unstack(ret_etf, :date, :ticker, :returns)
    ret_etf = ret_etf[1:Int(months_to_include), :]
    column_offset = 0
    while sum(ismissing.(Matrix(ret_etf[:, Not(:date)]))) > 0
        fill_in_least_missing!(ret_etf, column_offset, scenario_offset)
        column_offset += 1
    end
    # clean up 
    return hcat(ret_etf[:, [:date]], DataFrame(Float64.(Matrix(ret_etf[:, Not(:date)])), names(ret_etf)[2:end]))
end
end