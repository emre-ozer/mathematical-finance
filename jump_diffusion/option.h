#include <iostream>
#include <math.h>
#include <string>

using namespace std;

int factorial(int n) {
    if (n == 0) { return 1; }
    else {
        int res = 1;
        for (int i = 1; i <= n; i++) {
            res *= i;
        }
        return res;
    }
}

double phi(double x) {
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

double phi_prime(double x) {
    return exp(-x*x/2) / sqrt(2*M_PI);
}

class Option {
    
    public:
    double S0, K, T, r, d, vol;
    string contract, exercise;
    Option(double=100.0,double=100.0,double=1.0,double=0.05,double=0.0,double=0.2,string="call",string="european");
    double analytic_price();
    double payoff(double);
    double update_rule(double,double,double);
    double analytic_jump_diffusion(double,double,double,int=10);
    double analytic_vega();
    double implied_volatility(double,double=0.2);
};

Option::Option(double spot, double strike, double expiry, double interest, double dividend, double volatility, string contract_type, string exercise_type) {
    /*
    Constructor for the Option class.
    Default is a European call option with strike and spot 100, expiry 1, interest 0.05, dividend 0.0 and volatility 0.2.
    */
    S0 = spot;
    K = strike;
    T = expiry;
    r = interest;
    d = dividend;
    vol = volatility;
    contract = contract_type;
    exercise = exercise_type;
}

double Option::payoff(double S) {
    try {
        if (contract == "call") {
            return max(S - K, 0.0); 
        }
        else if (contract == "put") {
            return max(K - S, 0.0); 
        }
        else { throw "Unsupported contract_type"; }
    }
    catch (...) {
        cout << "Option::payoff supports only options of type 'call' or 'put'." << endl;
        exit(-1);
    }
}

double Option::update_rule(double node_up, double node_down, double S) {
    try {
        if (exercise == "european") {
            return 0.5 * (node_up + node_down);
        }
        else if (exercise == "american") {
            return max(0.5 * (node_up + node_down), payoff(S));
        }
        else { throw "Unsupported exercise_type"; }
    }
    catch (...) {
        cout << "Option::update_rule supports only options of exercise 'european' or 'american'." << endl;
        exit(-1);
    }
}

double Option::analytic_price() {
    double a1 = (log(S0/K) + (r-d+0.5*vol*vol)*T) / (vol*sqrt(T));
    double a2 = (log(S0/K) + (r-d-0.5*vol*vol)*T) / (vol*sqrt(T));
    try {
        if (exercise == "european" && contract == "call") {
            return exp(-d*T)*S0*phi(a1) - exp(-r*T)*K*phi(a2);
        }
        else if (exercise == "european" && contract == "put") {
            return -exp(-d*T)*S0*phi(-a1) + exp(-r*T)*K*phi(-a2);
        }
        else { throw "Unsupported option type"; }
    }
    catch (...) {
        cout << "Option::analytic_price supports only European put and call options." << endl;
        exit(-1);
    }
}

double Option::analytic_jump_diffusion(double m, double nu, double lam, int cutoff) {
    /*
    Returns the analytic jump-diffusion price for log-normal jump distributions. 
    The jumps are distributed as J \sim m exp(-1/2 nu^2 + nu Z) with Z random normal variable.
    Only works for options for which there exists an analytic price (currently, only European vanilla call and put are implemented.)
    Cutoff determines after which term the infinite series is terminated. By default, it is set to 10.
    */
    double S_i = S0; // initial spot
    double vol_i = vol; // initial volatility
    double mu = r - lam * (m - 1.0); // risk-free rate
    double sum = 0;
    for (int n = 0; n <= cutoff; n++) {
        S0 = S_i * pow(m,n) * exp((mu-r)*T);
        vol = sqrt(vol_i*vol_i + n*nu*nu/T);
        sum += exp(-lam*T) * pow(lam*T, n) / factorial(n) * analytic_price();
    }
    // restoring initial values
    S0 = S_i;
    vol = vol_i;
    return sum;
}

double Option::analytic_vega() {
    double a1 = (log(S0/K) + (r-d+0.5*vol*vol)*T) / (vol*sqrt(T));
    double a2 = (log(S0/K) + (r-d-0.5*vol*vol)*T) / (vol*sqrt(T));
    double b1 = -a1/vol + sqrt(T);
    double b2 = -a2/vol - sqrt(T);
    try {
        if (exercise == "european" && contract == "call") {
            return exp(-d*T)*S0*phi_prime(a1)*b1 - exp(-r*T)*K*phi_prime(a2)*b2;
        }
        else if (exercise == "european" && contract == "put") {
            return exp(-d*T)*S0*phi_prime(a1)*b1 - exp(-r*T)*K*phi_prime(a2)*b2;
        }
        else { throw "Unsupported option type"; }
    }
    catch (...) {
        cout << "Option::analytic_vega supports only European put and call options." << endl;
        exit(-1);
    }
}

double Option::implied_volatility(double market_price, double guess) {
    double vol_i = vol; // store initial volatility
    int n_guess = 0;
    double err;
    do {
        vol = guess;
        err = (analytic_price()-market_price) / analytic_vega();
        guess -= err;
        n_guess++;
    } while(abs(err) > 1e-3 || n_guess < 15);
    vol = vol_i; // restore initial volatility
    return guess;
}