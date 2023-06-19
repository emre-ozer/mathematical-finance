import datetime as dt
import numpy as np
from math import *
from scipy.optimize import fsolve
from time import sleep

def phi(x):
    #'Cumulative distribution function for the standard normal distribution'
    return (1.0 + erf(x / sqrt(2.0))) / 2.0

def IV_call(vol,S,K,T,r,d,market_price):
    a1 = (log(S/K) + (r - d + 0.5*vol*vol)*T) / (vol*sqrt(T))
    a2 = (log(S/K) + (r - d - 0.5*vol*vol)*T) / (vol*sqrt(T))
    return exp(-d*T)*S*phi(a1) - K*exp(-r*T)*phi(a2) - market_price

def implied_vol(vol_guess, params):
    keys, vals = zip(*params.items())
    return fsolve(IV_call, vol_guess, args=vals)[0]

"""
Assume we fetch the following information:
- Strike
- Spot
- Interest rate
- Expiry date (as datetime.date object)
- Query date (as datetime.date object)
- Market price of option
- Dividend yield

Then, we compute the following:
- Time to expiry
- Implied volatility
- Greeks

One complication: the dividend payments are made on specific dates, and one must take into account how this works with respect to different expiries.
"""


"""
POLYGON STUFF FROM HERE

Here, we fetch the information used above.
"""
import polygon

reference_client = polygon.ReferenceClient('2ufZ2jOxhppcuNXiPM3XOuG3WiVizkVs')
stocks_client = polygon.StocksClient('2ufZ2jOxhppcuNXiPM3XOuG3WiVizkVs')
options_client = polygon.OptionsClient('2ufZ2jOxhppcuNXiPM3XOuG3WiVizkVs')

def approx_dividend_rate(ticker, current_price):
    """
    The reason we are approximating is to minimize the number of API queries.
    Query number = 1
    """
    report = reference_client.get_stock_dividends(ticker, order='desc', limit=1).get('results')
    if not report:
        return 0
    else:
        report = report[0]
    freq = report.get('frequency')
    div_amount = report.get('cash_amount')
    return freq * div_amount / current_price

"""
strike = float(input('Strike: '))
opt_type = 'call'
ticker = 'AAPL'
spot = stocks_client.get_previous_close(ticker).get('results')[0].get('vw')
dividend = approx_dividend_rate(ticker, spot)
expiry_date = dt.date(2023,8,18)
previous_date = dt.date.today() - dt.timedelta(days=1)
t = (expiry_date - previous_date).days / 365
opt_ticker = polygon.build_option_symbol(ticker, expiry_date, opt_type, strike)
opt_price = options_client.get_previous_close(opt_ticker).get('results')[0].get('vw')


params = {'S':spot, 'K':strike, 'T':t, 'r':0.0518 ,'d':dividend, 'market_price':opt_price}
"""

strikes = np.arange(100,171,5)

def construct_smile(ticker):
    query_count = 0
    opt_type = 'call'
    spot = stocks_client.get_previous_close(ticker).get('results')[0].get('vw')
    dividend = approx_dividend_rate(ticker, spot)
    expiry_date = dt.date(2023,8,18)
    previous_date = dt.date.today() - dt.timedelta(days=1)
    t = (expiry_date - previous_date).days / 365
    f = open("SM_{0}.csv".format(ticker), "w")
    f.write('K,IV')
    for K in strikes:
        if (query_count == 5):
            query_count = 0
            print('Query limit reached, waiting for 60 seconds.')
            sleep(60)
        opt_ticker = polygon.build_option_symbol(ticker, expiry_date, opt_type, K)
        opt_price = options_client.get_previous_close(opt_ticker)
        if (opt_price.get('results') == None):
            print('Could not fetch option previous close for {0}'.format(opt_ticker))
        else :
            opt_price = opt_price.get('results')[0].get('vw')
            params = {'S':spot, 'K':K, 'T':t, 'r':0.0518 ,'d':dividend, 'market_price':opt_price}
            output = str(K) + ',' + str(implied_vol(0.3,params))
            f.write('\n' + output)
        query_count += 1
    f.close()

construct_smile('GOOG')

"""
print(params)
output = str(strike) + ',' + str(implied_vol(0.3,params))
print(output)

with open("smile.csv", "a") as myfile:
    myfile.write('\n' + output)
"""


"""
TO DO: 
1 - Implement Greeks, maybe do OOP (for snapshot)
2 - Implemenet systematic way to compute smile
"""