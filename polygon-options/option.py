import datetime as dt
import numpy as np
from math import *
from scipy.optimize import fsolve
from time import sleep
import polygon

client_ID = '2ufZ2jOxhppcuNXiPM3XOuG3WiVizkVs'

reference_client = polygon.ReferenceClient(client_ID)
stocks_client = polygon.StocksClient(client_ID)
options_client = polygon.OptionsClient(client_ID)

class Option:

    def __init__(self, underlying, strike, type, expiry_date):
        """
        underlying: ticker of underlying stock (str)
        strike: strike price (float)
        type: call or put ('C' or 'P')
        expiry_date: datetime.date object, expiry date of contract
        """
        self.underlying = underlying
        self.strike = strike
        self.expiry_date = expiry_date
        self.type = type
        previous_date = dt.date.today() - dt.timedelta(days=1)
        self.T = t = (expiry_date - previous_date).days / 365
        self.spot = stocks_client.get_previous_close(underlying).get('results')[0].get('vw')
        self.d = self.__approx_dividend_rate(underlying, self.spot)
        self.r = 0.0518
        self.ticker = polygon.build_option_symbol(underlying, expiry_date, type, strike)
        self.exercise_style = reference_client.get_option_contract(self.ticker).get('results').get('exercise_style')
        self.market_price = options_client.get_previous_close(self.ticker).get('results')[0].get('vw')
        self.implied_vol = self.__compute_IV(0.2)
    
    def __phi(self, x):
        """
        Cumulative distribution function for the standard normal distribution
        """
        return (1.0 + erf(x / sqrt(2.0))) / 2.0
    
    def __phiprime(self, x):
        return exp(-x*x/2)/sqrt(2*pi)

    def __BS_call(self, vol):
        S = self.spot
        K = self.strike
        r = self.r
        d = self.d
        T = self.T
        a1 = (log(S/K) + (r - d + 0.5*vol*vol)*T) / (vol*sqrt(T))
        a2 = (log(S/K) + (r - d - 0.5*vol*vol)*T) / (vol*sqrt(T))
        return exp(-d*T)*S*self.__phi(a1) - K*exp(-r*T)*self.__phi(a2) - self.market_price

    def __BS_put(self, vol):
        S = self.spot
        K = self.strike
        r = self.r
        d = self.d
        T = self.T
        a1 = (log(S/K) + (r - d + 0.5*vol*vol)*T) / (vol*sqrt(T))
        a2 = (log(S/K) + (r - d - 0.5*vol*vol)*T) / (vol*sqrt(T))
        return - exp(-d*T)*S*self.__phi(-a1) + K*exp(-r*T)*self.__phi(-a2) - self.market_price

    def __compute_IV(self, vol_guess):
        if (self.type == 'C'):
            return fsolve(self.__BS_call, vol_guess)[0]
        else:
            return fsolve(self.__BS_put, vol_guess)[0]
    
    def __approx_dividend_rate(self, ticker, current_price):
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
    
    def __Delta(self):
        S = self.spot
        K = self.strike
        r = self.r
        d = self.d
        T = self.T
        vol = self.implied_vol
        a1 = (log(S/K) + (r - d + 0.5*vol*vol)*T) / (vol*sqrt(T))
        a2 = (log(S/K) + (r - d - 0.5*vol*vol)*T) / (vol*sqrt(T))
        if self.type == 'C':
            return exp(-d*T)*self.__phi(a1) + (exp(-d*T)*S*exp(-a1*a1/2) - exp(-r*T)*K*exp(-a2*a2/2)) / (S*vol*sqrt(2*pi*T))
        else:
            return exp(-d*T)*self.__phi(a1) + (exp(-d*T)*S*exp(-a1*a1/2) - exp(-r*T)*K*exp(-a2*a2/2)) / (S*vol*sqrt(2*pi*T)) - exp(-d*T)

    def __Gamma(self):
        S = self.spot
        K = self.strike
        r = self.r
        d = self.d
        T = self.T
        vol = self.implied_vol
        a1 = (log(S/K) + (r - d + 0.5*vol*vol)*T) / (vol*sqrt(T))
        a2 = (log(S/K) + (r - d - 0.5*vol*vol)*T) / (vol*sqrt(T))
        return 1/(S*vol*sqrt(2*pi*T))*( exp(-d*T)*exp(-a1*a1/2) - a1*exp(-d*T)*exp(-a1*a1/2)/(vol*sqrt(T)) + K/S*exp(-r*T)*exp(-a2*a2/2) + K/S*a2*exp(-r*T)*exp(-a2*a2/2)/(vol*sqrt(T)) )
    
    def __Vega(self):
        S = self.spot
        K = self.strike
        r = self.r
        d = self.d
        T = self.T
        vol = self.implied_vol
        a1 = (log(S/K) + (r - d + 0.5*vol*vol)*T) / (vol*sqrt(T))
        a2 = (log(S/K) + (r - d - 0.5*vol*vol)*T) / (vol*sqrt(T))
        b1 = -a1/vol + sqrt(T)
        b2 = -a2/vol - sqrt(T)
        return exp(-d*T)*S*self.__phiprime(a1)*b1 - exp(-r*T)*K*self.__phiprime(a2)*b2

    def __Rho(self):
        S = self.spot
        K = self.strike
        r = self.r
        d = self.d
        T = self.T
        vol = self.implied_vol
        a1 = (log(S/K) + (r - d + 0.5*vol*vol)*T) / (vol*sqrt(T))
        a2 = (log(S/K) + (r - d - 0.5*vol*vol)*T) / (vol*sqrt(T))
        rho_C = sqrt(T)/vol*(exp(-d*T)*S*self.__phiprime(a1) - exp(-r*T)*K*self.__phiprime(a2)) + T*exp(-r*T)*K*self.__phi(a2)
        if self.type == 'C':
            return rho_C
        else:
            return rho_C - T*exp(-r*T)*K

    def __Theta(self):
        S = self.spot
        K = self.strike
        r = self.r
        d = self.d
        T = self.T
        vol = self.implied_vol
        a1 = (log(S/K) + (r - d + 0.5*vol*vol)*T) / (vol*sqrt(T))
        a2 = (log(S/K) + (r - d - 0.5*vol*vol)*T) / (vol*sqrt(T))
        c1 = -0.5*log(S/K)/vol/T**(3/2) + 0.5*(r-d+0.5*vol*vol)/vol/sqrt(T)
        c2 = -0.5*log(S/K)/vol/T**(3/2) + 0.5*(r-d-0.5*vol*vol)/vol/sqrt(T)
        theta_C = d*exp(-d*T)*S*self.__phi(a1) - r*exp(-r*T)*K*self.__phi(a2) - exp(-d*T)*S*self.__phiprime(a1)*c1 + exp(-r*T)*K*self.__phiprime(a2)*c2
        if self.type == 'C':
            return theta_C
        else:
            return theta_C - d*exp(-d*T)*S + r*exp(-r*T)*K

    def snapshot(self):
        details = {"contract_type": 'call' if self.type == 'C' else 'put',
                   "exercise_style": self.exercise_style,
                   "expiration_date": self.expiry_date.strftime('%Y-%m-%d'),
                   "strike_price": self.strike,
                   "ticker": self.ticker}
        greeks = {"delta": self.__Delta(),
                  "gamma": self.__Gamma(),
                  "vega": self.__Vega(),
                  "rho": self.__Rho(),
                  "theta": self.__Theta()}
        return {"details": details, 
                "implied_volatility": self.implied_vol,
                 "greeks" : greeks }

opt1 = Option('AAPL', 200, 'C', dt.date(2023,8,18))
print(opt1.snapshot())