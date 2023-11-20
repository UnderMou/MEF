#include <cmath>

using namespace std;

// Reference : https://pomax.github.io/bezierinfo/legendre-gauss.html

vector<vector<double>> shl(int nen, int nint){
    vector<double> pt(nint);
    vector<double> w(nint);

    vector<vector<double>> shg(nint, vector<double>(nen));
    
    if(nint == 5){
        pt[0] = 0.0;
        pt[1] = -0.5384693101056831;
        pt[2] = 0.5384693101056831;
        pt[3] = -0.9061798459386640;
        pt[4] = 0.9061798459386640;
    }

    double t;
    for (int l = 0; l < nint; l++){
        t = pt[l];

        if (nen == 5){
            shg[0][l] = (2.0/3.0)*(t+1.0/2.0)*t*(t-1.0/2.0)*(t-1.0);
            shg[1][l] = (-8.0/3.0)*(t+1.0)*t*(t-1.0/2.0)*(t-1.0);
            shg[2][l] = 4.0*(t+1.0)*(t+1.0/2.0)*(t-1.0/2.0)*(t-1.0);
            shg[3][l] = (-8.0/3.0)*(t+1.0)*t*(t+1.0/2.0)*(t-1.0);
            shg[4][l] = (2.0/3.0)*(t+1.0)*(t+1.0/2.0)*t*(t-1.0/2.0);
        }

        // TODO: Implementar casos quadráticos, cúbicos, etc. ...

    }

    return shg;

}