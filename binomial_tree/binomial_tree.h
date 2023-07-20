#include "option.h"

using namespace std;

class BinomialTree {
        Option p;
        double S0, T, r, d, vol, dt;
        int N;
    public:
        double price;
        BinomialTree (Option, int);
        void set_option(Option opt) {
            p = opt;
            S0 = opt.S0;
            T = opt.T;
            r = opt.r;
            d = opt.d;
            vol = opt.vol;
            dt = T / N;
        }
        Option get_option() {
            return p;
        }
        void set_N(int m) {
            N = m;
            dt = (double) T / N;
        }
        void set_T(double t) {
            T = t;
            dt = T / N;
        }
        double compute_price();
        double compute_averaged_price();
};

BinomialTree::BinomialTree (Option opt, int n) {
    N = n;
    dt = (double) opt.T / N;
    set_option(opt);
}

double BinomialTree::compute_price() {
    /*
    Constructs the final layer using the option payoff.
    Applies the update rule recursively to the first node.
    Stores the initial node value as price, and outputs it.
    */
    double S_min, S;
    double * nodes = new double[N+1];
    // compute final layer
    S_min = S0 * exp( (r - d - 0.5*vol*vol)*T - N*vol*sqrt(dt) );
    for (int i = 0; i < N+1; i++) {
        S = S_min * exp(2*i*vol*sqrt(dt));
        nodes[i] = p.payoff(S);
    }
    // apply update rule backwards
    for (int m = N; m > 0; m--) {
        S_min = S0 * exp( (r - d - 0.5*vol*vol)*(m-1)*dt - (m-1)*vol*sqrt(dt) );
        for (int i = 0; i < m; i++) {
            S = S_min * exp(2*i*vol*sqrt(dt));
            nodes[i] = p.update_rule(exp(-r*dt) * nodes[i], exp(-r*dt) * nodes[i+1], S);
        }
    }
    price = nodes[0];
    delete[] nodes;
    return price;
}

double BinomialTree::compute_averaged_price() {
    /*
    Computes the even-odd averaged price, for a given N outputs the averaged value of N and N+1 layers.
    */
    double price1 = compute_price();
    set_N(N+1);
    double price2 = compute_price();
    set_N(N-1);
    price = 0.5 * (price1 + price2); 
    return price;
}