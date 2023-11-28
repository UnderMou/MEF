#include <cmath>

using namespace std;

// Reference : https://pomax.github.io/bezierinfo/legendre-gauss.html

vector<vector<double>> shl_PG(int nen, int nint, double xe1, double xe2, double epsilon){
    vector<double> pt(nint);
    vector<double> w(nint);

    vector<vector<double>> shg(nen, vector<double>(nint));
    
    
    pt[0] = 0.0;
    pt[1] = -0.5384693101056831;
    pt[2] = 0.5384693101056831;
    pt[3] = -0.9061798459386640;
    pt[4] = 0.9061798459386640;

    double tx;
    for (int l = 0; l < nint; l++){
        tx = ((xe2 - xe1)/2)*pt[l] + 0.5*(xe2 + xe1);

        shg[0][l] = -(exp((tx-xe2)/sqrt(epsilon)) - exp((xe2 - tx)/sqrt(epsilon))) / (exp((xe2 - xe1)/sqrt(epsilon)) - exp((xe1 - xe2)/sqrt(epsilon)));
        shg[1][l] = (exp((tx-xe1)/sqrt(epsilon)) - exp((xe1 - tx)/sqrt(epsilon))) / (exp((xe2 - xe1)/sqrt(epsilon)) - exp((xe1 - xe2)/sqrt(epsilon)));

    }

    return shg;

}