#include <cmath>

using namespace std;

// Reference : https://pomax.github.io/bezierinfo/legendre-gauss.html

vector<vector<double>> dshl(int nen, int nint){
    vector<double> pt(nint);
    vector<double> w(nint);

    vector<vector<double>> dshg(nen, vector<double>(nint));
    
    pt[0] = 0.0;
    pt[1] = -0.5384693101056831;
    pt[2] = 0.5384693101056831;
    pt[3] = -0.9061798459386640;
    pt[4] = 0.9061798459386640;

    double t;
    for (int l = 0; l < nint; l++){
        t = pt[l];
        
        dshg[0][l] = -1.0/2.0;
        dshg[1][l] = 1.0/2.0;
        
    }

    return dshg;

}