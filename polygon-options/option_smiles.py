import option
import datetime as dt
import numpy as np
import time

expiry = dt.date(2023,10,20)
underlying = 'SPY'
strikes = np.arange(420, 561, 5)
type = 'C'

def smile(underlying, strikes, expiry, type):
    query_count = 0
    f = open("SM_{0}{1}{2}.csv".format(underlying, type, expiry.strftime("%y%m%d")), "w")
    f.write('K,IV')
    for K in strikes:
        print(K)
        if query_count == 5:
            print('Query limit reached, waiting for 60 seconds.')
            time.sleep(60)
            query_count = 0
        opt = option.Option(underlying, K, type, expiry)
        output = str(K) + ',' + str(opt.implied_vol)
        f.write('\n' + output)
        query_count += 1
    f.close()

smile(underlying, strikes, expiry, type)

 