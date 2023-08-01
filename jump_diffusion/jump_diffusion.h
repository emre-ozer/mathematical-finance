#include "option.h"
#include <random>

int poisson_draw(double lam) {
    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<double> uni(0.0, 1.0);
    double L = exp(-lam), p = 1.0;
    int k = 0;
    do {
        k += 1;
        p *= uni(gen);
    } 
    while (p > L);
    return k-1;
}

vector<int> poisson_sample(double lam, int N) {
    vector<int> sample;
    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<double> uni(0.0, 1.0);
    double L = exp(-lam);
    double p;
    int k;
    for (int i = 0; i < N; i++) {
        p = 1.0;
        k = 0;
        do {
            k += 1;
            p *= uni(gen);
        } 
        while (p > L);
        sample.push_back(k-1);
    }
    return sample;
}

double moro(double x) {
    /*
    Inverse cumulative normal function is defined for x in [0,1]. Below, x denotes the argument of the function. I know the implementation is slow, I will make it faster if it comes to it.
    */
    const double a[] =
    {   2.50662823884,
        -18.61500062529,
        41.39119773534,
        -25.44106049637
    };
    const double b[] =
    {   -8.47351093090,
        23.08336743743,
        -21.06224101826,
        3.13082909833
    };
    const double c[] =
    {   0.3374754822726147,
        0.9761690190917186,
        0.1607979714918209,
        0.0276438810333863,
        0.0038405729373609,
        0.0003951896511919,
        0.0000321767881768,
        0.0000002888167364,
        0.0000003960315187
    };
    double y = x-0.5;
    double r;
    if (fabs(y) < 0.42) {
        r = y*y;
        double A = 0;
        double B = 0;
        for (int j = 0; j<4; j++) {
            A += a[j] * pow(r,j);
            B += b[j] * pow(r,j+1);
        }
        return y*A / (B + 1.0);
    }
    else {
        if (y < 0) { r = x; }
        else { r = 1-x; }
        double s = log(-log(r));
        double t = 0;
        for (int j = 0; j < 9;j++) {
            t += c[j]*pow(s,j);
        }
        if (x > 0.5) { return t; }
        else { return -t; }
    }
}

vector<double> normal_sample(int N) {
    /*
    Given sample size N, generates an N dimensional vector of normally distributed random numbers. Uses the inverse cumulative distribution, approximated by the Moro method.
    */
    vector<double> sample;
    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<double> uni(0.0, 1.0);

    for (int i = 0; i < N; i++) {
        sample.push_back(moro(uni(gen)));
    }
    return sample;
}

double jump_diffusion_MC(Option opt, double m, double nu, double lam, int N) {
    vector<int> P = poisson_sample(lam*opt.T, N);
    vector<double> W = normal_sample(N);
    double mu = opt.r - lam*(m-1);
    double logST;
    double price = 0;
    for (int i = 0; i < N; i++) {
        logST = log(opt.S0);
        logST += mu*opt.T;
        logST += -0.5*opt.vol*opt.vol*opt.T;
        logST += P[i]*log(m);
        logST += -0.5*nu*nu*P[i];
        logST += sqrt(opt.vol*opt.vol*opt.T + P[i]*nu*nu)*W[i];
        price += exp(-opt.r*opt.T)*opt.payoff(exp(logST));
    }
    return price / N;
}