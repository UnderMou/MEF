// function sh = shl(nen,nint)
//     if(nint == 2)
//         pt(1) = -sqrt(3.)/3.;
//         pt(2) = sqrt(3.)/3.;
//         w(1) = 1.;
//         w(2) = 1.;
//     end

//     for l=1:nint
//         t=pt(l);
//         if(nen==2)
//             sh(1,l) = (1.0-t)/2.0;
//             sh(2,l) = (1.0+t)/2.0;
//         end
//     end
// end


#include <cmath>

using namespace std;

// Reference : https://pomax.github.io/bezierinfo/legendre-gauss.html

vector<vector<float>> shl(int nen, int nint){
    vector<float> pt(nint);
    vector<float> w(nint);

    vector<vector<float>> shg(nint, vector<float>(nen));
    
    if(nint == 2){
        pt[0] = -0.5773502691896257;
        pt[1] = 0.5773502691896257;
        w[0] = 1.0;
        w[1] = 1.0;
    }

    if(nint == 3){
        pt[0] = 0.0;
        pt[1] = -0.7745966692414834;
        pt[2] = 0.7745966692414834;
        w[0] = 0.8888888888888888;
        w[1] = 0.5555555555555556;
        w[2] = 0.5555555555555556;
    }

    if(nint == 4){
        pt[0] = -0.3399810435848563;
        pt[1] = 0.3399810435848563;
        pt[2] = -0.8611363115940526;
        pt[3] = 0.8611363115940526;
        w[0] = 0.6521451548625461;
        w[1] = 0.6521451548625461;
        w[2] = 0.3478548451374538;
        w[3] = 0.3478548451374538;
    }

    if(nint == 5){
        pt[0] = 0.0;
        pt[1] = -0.5384693101056831;
        pt[2] = 0.5384693101056831;
        pt[3] = -0.9061798459386640;
        pt[4] = 0.9061798459386640;
        w[0] = 0.5688888888888889;
        w[1] = 0.4786286704993665;
        w[2] = 0.4786286704993665;
        w[3] = 0.2369268850561891;
        w[4] = 0.2369268850561891;
    }

    //     for l=1:nint
    //         t=pt(l);
    //         if(nen==2)
    //             sh(1,l) = (1.0-t)/2.0;
    //             sh(2,l) = (1.0+t)/2.0;
    //         end
    //     end

    float t;
    for (int l = 0; l < nint; l++){
        t = pt[l];

        if (nen == 2){
            shg[0][l] = (1.0-t)/2.0;
            shg[1][l] = (1.0+t)/2.0;
        }

        // TODO: Implementar casos quadráticos, cúbicos, etc. ...

    }

    return shg;

}