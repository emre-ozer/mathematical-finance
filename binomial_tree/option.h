/*
IS THIS THE DEFINITIVE OPTION CLASS????
*/
#include <iostream>
#include <math.h>
#include <string>

using namespace std;

class Option {
    public:
    double S0, K, T, r, d, vol;
    string contract, exercise;
    Option(double=100.0,double=100.0,double=1.0,double=0.05,double=0.0,double=0.2,string="call",string="european");
    double payoff(double);
    double update_rule(double,double,double);
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
        cout << "Option::BT_update_rule supports only options of exercise 'european' or 'american'." << endl;
        exit(-1);
    }
}