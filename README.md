# Investment Management
## Disclaimer
This project is not intended to give financial advice.
Rather it calculates optimal diversified portfolios based on historical data and forecasts based on that historical dat 
and is not intended as a guarantee of future perofrmance for any individual fund or portfolio of funds.
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
```
  a. Downloading ticker data is slow unless you spring for a premium alhpavantage key. The free version is limited to 5 requests per minute
  b. You may run into errors if market data has been updated for the most recent month end but treasury data has not.
  This process should be stable the second or third markt day of each month but might not on the very first because of data disconnects
5. To run the process with different tickers you can modify the Tickers.csv file to reflect what funds you want to examine
6. To rerun the process (can be done at month end to incorporate most recent market data), delete the following files:
```
# downloaded from alphavantage. Deleting these triggers a rerun of the download of ETFs
data/raw/full_etf_prices.csv
data/raw/full_yields.csv
# results from prior model runs
data/model_results/constr_opt_eff_front.csv
data/model_results/constr_opt_eff_front_forecast.csv
```
Run `make docker-run-main` when you are ready to rerun
