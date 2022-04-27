module LoadData
using AlphaVantage
using DataFrames, Dates
using ProgressBars
import CSV

function load_data()
    # Get daily S&P 500 data
    tickers = String.(CSV.read("Tickers.csv", DataFrame)[:, :ticker])
    raw_month_etf_resp = []
    for i in ProgressBar(1:length(tickers))
        resp = time_series_monthly_adjusted(tickers[i])
        push!(raw_month_etf_resp, resp)
    end
    df_month_etf = [DataFrame(raw_month_etf_resp[i]) for i in 1:length(raw_month_etf_resp)]
    [df_month_etf[i][!, :ticker] .= tickers[i] for i in 1:length(tickers)]
    raw_month_yield_resp = []
    maturities = ["3month", "5year", "10year", "30year"]
    for i in ProgressBar(1:length(maturities))
        resp = treasury_yield("monthly", maturities[i])
        push!(raw_month_yield_resp, resp)
    end
    df_month_yield = [DataFrame(raw_month_yield_resp[i]) for i in 1:length(raw_month_yield_resp)]
    [df_month_yield[i][!, :maturity] .= maturities[i] for i in 1:length(maturities)]
    full_yield = vcat(df_month_yield...)
    full_etf = vcat(df_month_etf...)
    CSV.write("data/raw/full_etf_prices.csv", full_etf)
    CSV.write("data/raw/full_yields.csv", full_yield)
end
end