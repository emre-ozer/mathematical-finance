#include "option.h"

using namespace std;

class BinomialTree {
        Option p;
        double S0, T, r, d, vol, dt;
        int N;
        double *final_layer;
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
        }
        void set_N(int m) {
            N = m;
            dt = T / N;
        }
        void compute_final_layer();
        double compute_price(bool = true);
        double compute_averaged_price();
};

BinomialTree::BinomialTree (Option opt, int n) {
    set_option(opt);
    N = n;
    dt = T / N;
}

void BinomialTree::compute_final_layer() {
    /*
    Uses the option payoff to initialize the final layer of the tree. This is stored in the final_layer array.
    */
    double S;
    double nodes[N+1];
    double S_min = S0 * exp( (r - d - 0.5*vol*vol)*T - N*vol*sqrt(dt) );
    for (int i = 0; i < N+1; i++) {
        S = S_min * exp(2*i*vol*sqrt(dt));
        nodes[i] = p.payoff(S);
    }
    final_layer = nodes;
}

double BinomialTree::compute_price(bool update_price) {
    /*
    Constructs the final layer using the compute_final_layer function.
    Applies the update rule recursively to the first node.
    Stores the initial node value as price, and outputs it.
    */
    compute_final_layer();
    double *nodes = final_layer;
    double S_min, S;
    for (int m = N; m > 0; m--) {
        S_min = S0 * exp( (r - d - 0.5*vol*vol)*(m-1)*dt - (m-1)*vol*sqrt(dt) );
        for (int i = 0; i < m; i++) {
            S = S_min * exp(2*i*vol*sqrt(dt));
            nodes[i] = p.update_rule(exp(-r*dt) * nodes[i], exp(-r*dt) * nodes[i+1], S);
        }
    }
    if (update_price) { price = nodes[0]; }
    return nodes[0];
}

double BinomialTree::compute_averaged_price() {
    /*
    Computes the even-odd averaged price, for a given N outputs the averaged value of N and N+1 layers.
    */
    double price1 = compute_price(false);
    set_N(N+1);
    double price2 = compute_price(false);
    set_N(N-1);
    price = 0.5 * (price1 + price2); 
    return price;
}
