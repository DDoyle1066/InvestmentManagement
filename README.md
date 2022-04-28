# Investment Management
This project pulls monthly market data for a set of specified indices, calculates the efficient frontiers based on historical data and forecast data and reports that information back
## Disclaimer
This project is not intended to give financial advice.
Rather it calculates optimal diversified portfolios based on historical data and forecasts based on that historical data 
and is not intended as a guarantee of future performance for any individual fund or portfolio of funds.
## Installation and runtime instructions
1. Install docker and git
2. Clone repo
```
git clone https://github.com/DDoyle1066/InvestmentManagement/
cd ./InvestmentManagement/
```
3. Get an AlphaVantage key from https://www.alphavantage.co/support/#api-key and replace the environment variable on line 37 of the Dockerfile with your key
4. Make docker containter, spin it up and run process
```
make docker-build
make docker-up
make docker-run-main
# output can be found at model_results/* It is in the form of portfolio allocation by ticker along with summary statistics for that portfolio on the left hand side
# The *_agg versions aggregate information acorss the scenarios found in the other files. the *_forceast file contains information from the forecast data methodology
```
  - Downloading ticker data is slow unless you spring for a premium alhpavantage key. The free version is limited to 5 requests per minute
  - You may run into errors if market data has been updated for the most recent month end but treasury data has not.
  This process should be stable the second or third market day of each month but might not on the very first because of data disconnects
5. To run the process with different tickers you can modify the Tickers.csv file to reflect what funds you want to examine
6. To rerun the process (can be done at month end to incorporate most recent market data), delete the following files:
```
# downloaded from alphavantage. Deleting these triggers a rerun of the download of ETFs
data/raw/full_etf_prices.csv
# results from prior model runs
data/model_results/constr_opt_eff_front.csv
data/model_results/constr_opt_eff_front_forecast.csv
```
Run `make docker-run-main` when you are ready to rerun
## Methodology
This project's steps are:
1. Download data from AlphaVantage and Treasury website
2. Calculate historical returns, remove indices with less than 5 years of data filter the remaining indices to only have about 15 years of data and for any indices with missing values (due to not existing for 15 years) backfill old returns by randomly sampling based on relationship between other indices where data exists
   1. For instance, if ETF A and ETF B are highly and statistically significantly correlated, ETF A has been around for 15 years but ETF B has only been around for 10, then returns for years 10-15 for ETF B will be filled based on ETF A's performance and its relationship with B along with random noise to not mess up correlations
3. Follow steps [here](https://www.newfrontieradvisors.com/media/1172/introduction-to-resampled-effiency.pdf) to calculate resampled effiicent frontier based on historical data
4. Create a process for forecasting future returns based on historical data, that involves:
    a. Forecasting changes in yields. This is done by running PCA on changes in the yield curve to avoid messcy correlations between different segments of the yield curve. 
    b. Estimating changes in ETFs based on changes in yield curve. This is done by regressing against a bond that solely invests in a single point in the yield curve for multiple points until no further significant relationship is found. Further, any regressions with an R<sup>2</sup> < 0.5 are discarded on the grounds that the predictions are likely spurious correlations and not reflective of systemic characteristics of that fund
    c. Using a Gaussian Mixtures model with 2 mixtures to estimate noise of residuals
        1.  This frameork has an appealing intuition where the market tends to rise slowly and steadily in one regime and crash quickly and violently in another
        2.  This can be improved by limiting correlations to avoid spurious correlaitions from limited data
    d. The main purpose of this step is remove bias from rate movements. If rates moved up or down bond funds with long durations will over or under perform relative to a long term expectation
5. Rerun 3, this time based on forecasted data from projected yields and projected returns given those yields  






