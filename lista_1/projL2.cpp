#include <iostream>
#include <vector>
#include "shl.cpp"
#include "we.cpp"

using namespace std;

vector<vector<float>> init_Matrix(int n1, int n2){
    vector<vector<float>> M(n1, vector<float>(n2));
    for (int i = 0; i < n1; i++){
        for (int j = 0; j < n2; j++){
            M[i][j] = 0.0;
        }
    }
    return M;
}

vector<float> init_F(int n){
    vector<float> F(n);
    for (int i = 0; i < n; i++) F[i] = 0.0;
    return F;
}

void print_Matrix(vector<vector<float>> M, int dim){
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            cout << M[i][j] << " ";
        }
        cout << "\n";
    }
}

void print_Vector(vector<float> F, int dim){
    for (int i = 0; i < dim; i++){
        cout << F[i] << "\n";
    }
}

int main(){

    // Defining the domain
    float a = -2.0;
    float b = 2.0;

    int nel = 4;        // number of elements
    int k = 1;          // polynomial degree
    int np = k*nel+1;   // mesh total nodes

    int nen = k+1;      // number of element nodes
    int nint = k+1;     // number of integration points

    float h = (b-a)/nel;    // element length

    vector<float> xl(np,0.0);
    
    xl[0] = a;
    for (int i = 1; i < xl.size(); i++) xl[i] = xl[i-1] + h;
    // for (int i = 0; i < xl.size(); i++) cout << xl[i] << "\n";


    // Global matrix and source vector
    vector<vector<float>> M(np, vector<float>(np));
    M = init_Matrix(np,np);
    print_Matrix(M,np);

    vector<float> F(np);
    F = init_F(np);
    print_Vector(F,np);

    vector<vector<float>> Me(nen, vector<float>(nen));
    Me = init_Matrix(nen,nen);

    vector<vector<float>> shg(nen, vector<float>(nint));
    shg = init_Matrix(nen,nint);

    vector<float> Fe(nen);
    Fe = init_F(nen);

    vector<float> w(nint);
    w = init_F(nint);

    float xx = 0.0;

    // Global problem construction
    for (int n = 0; n < nel; n++){
        
        Me = init_Matrix(nen,nen);
        Fe = init_F(nen);

        shg = shl(nen,nint);
        // print_Matrix(shg,nen);
        w = we(nint);
        // print_Vector(w,nint);

        for (int l = 0; l < nint; l++){
            xx = 0.0;

            for (int i = 0; i < nen; i++){
                //xx=xx+shg[i][l]*xl[n+i-1]; //funciona apenas para linear (k=1)
            }
            for (int j = 0; j < nen; j++){
                //Fe(j) = Fe(j) + f(xx)*shg(j,l)*w(l)*h/2.; 
                for (int i = 0; i < nen; i++){
                    // Me(i,j) = Me(i,j) + shg(i,l)*shg(j,l)*w(l)*h/2.;
                }
            }
        }

        // for (int j = 0; j < nen; j++){
        //     // F(n+j-1) = F(n+j-1) + Fe(j); %funciona apenas para linear (k=1)
        //     for (int i = 0; i < nen; i++){
        //         // M(n+i-1,n+j-1) = M(n+i-1,n+j-1) + Me(i,j); %funciona apenas para linear (k=1)
        //     }
        // }
    }
    // print_Matrix(M,np);
    // print_Vector(F,np);


    return 0;
}