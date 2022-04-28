module LoadData
using AlphaVantage
using DataFrames, Dates
using ProgressBars
import ExcelFiles
import CSV

function load_etfs(etf_path)
    # Get daily S&P 500 data
    tickers = String.(CSV.read("Tickers.csv", DataFrame)[:, :ticker])
    raw_month_etf_resp = []
    for i in ProgressBar(1:length(tickers))
        resp = time_series_monthly_adjusted(tickers[i])
        push!(raw_month_etf_resp, resp)
    end
    df_month_etf = [DataFrame(raw_month_etf_resp[i]) for i in 1:length(raw_month_etf_resp)]
    [df_month_etf[i][!, :ticker] .= tickers[i] for i in 1:length(tickers)]
    full_etf = vcat(df_month_etf...)
    CSV.write(etf_path, full_etf)
end

function load_yield_curves(input_urls, yield_label, output_file)
    download_dir = tempname()
    mkdir(download_dir)
    all_dfs = []
    for i in 1:length(input_urls)
        file_path = download(input_urls[i], "$download_dir/treas_$i.xls")
        rates_df = DataFrame(ExcelFiles.load(file_path, "Sheet1"))
        starting_row = 6
        starting_col = 3
        rates_mat = float.(rates_df[starting_row:end, starting_col:end]) |> Matrix |> transpose
        dates = Vector{Date}(undef, size(rates_mat)[1])
        month_lookup = Dict((monthabbr(x) => x < 10 ? "0" * string(x) : string(x) for x in 1:12))
        year = rates_df[3, starting_col]
        dates[1] = year * "-" * month_lookup[rates_df[4, starting_col]] |> Date |> lastdayofmonth
        for i in 2:length(dates)
            if !ismissing(rates_df[3, starting_col+i-1])
                year = rates_df[3, starting_col+i-1]
            end
            dates[i] = year * "-" * month_lookup[rates_df[4, starting_col+i-1]] |> Date |> lastdayofmonth
        end
        column_names = yield_label .* "Dur" .* string.(rates_df[starting_row:end, 1])
        return_df = hcat(DataFrame(Date=dates), DataFrame(rates_mat, column_names))
        push!(all_dfs, return_df)
    end
    full_df = vcat(all_dfs...)
    full_df = full_df[.!ismissing.(full_df[:, 2]), :]
    CSV.write(output_file, full_df)
end
end