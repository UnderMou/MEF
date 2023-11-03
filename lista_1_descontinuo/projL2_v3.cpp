#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

#include <Eigen/Dense>

#include "shl_v2.cpp"
#include "we_v2.cpp"
#include "pts_v2.cpp"

using namespace std;

vector<vector<double>> init_Matrix(int n1, int n2){
    vector<vector<double>> M(n1, vector<double>(n2));
    for (int i = 0; i < n1; i++){
        for (int j = 0; j < n2; j++){
            M[i][j] = 0.0;
        }
    }
    return M;
}

vector<double> init_F(int n){
    vector<double> F(n);
    for (int i = 0; i < n; i++) F[i] = 0.0;
    return F;
}

void print_Matrix(vector<vector<double>> M, int dim){
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            cout << fixed << setprecision(2) << M[i][j] << "\t";
        }
        cout << endl;
    }
}

void print_Vector(vector<double> F, int dim){
    for (int i = 0; i < dim; i++){
        cout << F[i] << "\n";
    }
}

double f(double xx){
    return sin(M_PI * xx)*sin(M_PI * xx);
}

Eigen::MatrixXd convert_matrix(vector<vector<double>> M, int dim){
    Eigen::MatrixXd M_eigen(dim,dim);

    for (int i = 0; i < dim; i++){
            for (int j = 0; j < dim; j++){
                M_eigen(i,j) = M[i][j];
            }
    }

    return M_eigen;
}

Eigen::MatrixXd convert_vector(vector<double> F, int dim){
    Eigen::MatrixXd F_eigen(dim,1);

    for (int j = 0; j < dim; j++){
        F_eigen(j,0) = F[j];
    }

    return F_eigen;
}

void projL2_analisys(int number_el){
       
    int nel = number_el;  // number of elements

    // Defining the domain
    double a = -2.0;
    double b = 2.0;
    
    int k = 4;          // polynomial degree
    int np = k*nel+nel;   // mesh total nodes

    int nen = k+1;      // number of element nodes
    int nint = k+1;     // number of integration points

    double h = (b-a)/nel;    // element length

    vector<double> xl(np,0.0);
    
    xl[0] = a;
    for (int i = 1; i < xl.size(); i++) {
        xl[i] = (i) % nen == 0 ? xl[i-1] : xl[i-1] + h/(nen-1);
    }
    // for (int i = 0; i < xl.size(); i++) cout << xl[i] << "\n";


    // Global matrix and source vector
    // vector<vector<double>> M(np, vector<double>(np));
    // M = init_Matrix(np,np);
    // print_Matrix(M,np);

    // vector<double> F(np);
    // F = init_F(np);
    // print_Vector(F,np);

    vector<vector<double>> Me(nen, vector<double>(nen));
    Me = init_Matrix(nen,nen);

    vector<vector<double>> shg(nint, vector<double>(nint));
    shg = init_Matrix(nint,nint);

    vector<double> Fe(nen);
    Fe = init_F(nen);

    vector<double> w(nint);
    w = init_F(nint);

    vector<double> pt(nint);
    pt = init_F(nint);


    vector<double> u_global(np);
    u_global = init_F(np);

    double xx = 0.0;

    // Global problem construction
    for (int n = 0; n < nel; n++){
        
        Me = init_Matrix(nen,nen);
        Fe = init_F(nen);

        shg = shl(nen,nint);
        // print_Matrix(shg,nen);
        w = we(nint);
        // print_Vector(w,nint);
        pt = pts(nint);
        // print_Vector(pt,nint);

        // cout << "Element " << n << endl;
        
        for (int l = 0; l < nint; l++){
            
            xx = h/2*pt[l] + 0.5*(xl[n*(nen) + nen-1] + xl[n*(nen)]);

            for (int j = 0; j < nen; j++){
                Fe[j] = Fe[j] + f(xx)*shg[j][l]*w[l]*h/2; 

                for (int i = 0; i < nen; i++){
                    Me[i][j] = Me[i][j] + shg[i][l]*shg[j][l]*w[l]*h/2;
                }

            }
        }
        // print_Matrix(Me,nen);
        // print_Vector(Fe,nen);
        // cout << "\n\n";

        // for (int j = 0; j < nen; j++){
        //     // if (j == nen-1) F[n+j] += F[j];
        //     // else F[n+j] = F[j];
        //     F[n*(nen-1)+j] += Fe[j];

        //     for (int i = 0; i < nen; i++){
        //         // if ((i == nen-1) && (j == nen-1) && (n!=nel-1)) {M[n+i][n+j] += Me[i][j]; cout << "OK\n";}
        //         // else {M[n+i][n+j] = Me[i][j];}
        //         M[n*(nen-1)+i][n*(nen-1)+j] += Me[i][j];
        //     }
        // }


        Eigen::MatrixXd M_eigen(nen,nen);
        M_eigen = convert_matrix(Me, nen);
        // std::cout << "Here is the matrix M_eigen:\n" << M_eigen << std::endl;

        Eigen::VectorXd F_eigen(nen);
        F_eigen = convert_vector(Fe, nen);
        // std::cout << "Here is the vector F_eigen:\n" << F_eigen << std::endl;

        Eigen::VectorXd u_eigen = M_eigen.colPivHouseholderQr().solve(F_eigen);
        // std::cout << "The solution is:\n" << u_eigen << std::endl;

        for (int i = 0; i < nen; i++){
            u_global[n*(nen) + i] = u_eigen(i);
        }

    }
    // print_Matrix(M,np);
    // print_Vector(F,np);

    
    


    stringstream ss;
    ss << "data_" << std::setw(4) << std::setfill('0') << nel << ".csv";
    string FileName = ss.str();

    ofstream csvFile(FileName);
    if (!csvFile.is_open()) {
        std::cerr << "Error opening the new CSV file." << std::endl;
        //return 1; // Return an error code
    }

    // Write data to the CSV file
    for (int i = 0; i < np; ++i) {
        csvFile << xl[i];
        if (i < np - 1) {
            csvFile << ","; // Use a comma as a delimiter
        }
    }
    csvFile << endl;
    for (int i = 0; i < np; i++) {
        csvFile << u_global[i];
        if (i < np - 1) {
            csvFile << ","; // Use a comma as a delimiter
        }
    }
    csvFile << endl;
    for (int i = 0; i < np; ++i) {
        csvFile << f(xl[i]);
        if (i < np - 1) {
            csvFile << ","; // Use a comma as a delimiter
        }
    }
    csvFile << "\n"; // Add a newline character to separate rows

    // Close the CSV file
    csvFile.close();

    std::cout << "CSV file written successfully." << std::endl;
}

int main(){

    const int size = 10; 
    int numb_el;

    // Initialize the array (optional)
    for (int i = 2; i < size; ++i) {
        numb_el = pow(2,i);
        cout << "Runing nel = " << numb_el << "\n";
        projL2_analisys(numb_el);
    }    

    return 0;
}