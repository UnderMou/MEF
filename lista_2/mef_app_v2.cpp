#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

#include <Eigen/Dense>

#include "shl_v2.cpp"
#include "dshl_v2.cpp"
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
            cout << fixed << setprecision(6) << M[i][j] << "\t";
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
    return M_PI*M_PI*sin(M_PI * xx);
    // return 1.0;
}

double u_exact(double xx){
    return sin(M_PI * xx);
    // return 1.0;
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

void mef_analisys(int number_el){
       
    int nel = number_el;  // number of elements

    // Defining the domain
    double a = 0.0;
    double b = 1.5;

    // boundary conditions
    // Dirichlet
    double BC_h = sin(M_PI*b);
    // Neumann
    double BC_g = -M_PI*cos(M_PI*a);
    cout << "Dirichlet h = " << BC_h << " | Neumann g = " << BC_g << endl;
    
    int k = 1;          // polynomial degree
    int np = k*nel+1;   // mesh total nodes

    int nen = k+1;      // number of element nodes
    int nint = k+1;     // number of integration points

    double h = (b-a)/nel;    // element length

    vector<double> xl(np,0.0);
    
    xl[0] = a;
    for (int i = 1; i < xl.size(); i++) xl[i] = xl[i-1] + h/(nen-1);
    // for (int i = 0; i < xl.size(); i++) cout << xl[i] << "\n";


    // Global K matrix
    vector<vector<double>> K(np, vector<double>(np));
    K = init_Matrix(np,np);
    // print_Matrix(M,np);

    // Global source vector
    vector<double> F(np);
    F = init_F(np);
    // print_Vector(F,np);

    vector<vector<double>> Ke(nen, vector<double>(nen));
    Ke = init_Matrix(nen,nen);

    vector<vector<double>> shg(nint, vector<double>(nint));
    shg = init_Matrix(nint,nint);

    vector<vector<double>> dshg(nint, vector<double>(nint));
    dshg = init_Matrix(nint,nint);

    vector<double> Fe(nen);
    Fe = init_F(nen);

    vector<double> w(nint);
    w = init_F(nint);

    vector<double> pt(nint);
    pt = init_F(nint);

    double xx = 0.0;

    // Global problem construction
    for (int n = 0; n < nel; n++){
        
        Ke = init_Matrix(nen,nen);
        Fe = init_F(nen);

        shg = shl(nen,nint);
        // print_Matrix(shg,nen);
        dshg = dshl(nen,nint);
        // print_Matrix(dshg,nen);
        w = we(nint);
        // print_Vector(w,nint);
        pt = pts(nint);
        // print_Vector(pt,nint);

        // cout << "Element " << n << endl;
        
        for (int l = 0; l < nint; l++){
            
            xx = h/2*pt[l] + 0.5*(xl[n*(nen-1) + nen-1] + xl[n*(nen-1)]);

            for (int j = 0; j < nen; j++){
                Fe[j] = Fe[j] + f(xx)*shg[j][l]*w[l]*h/2.0; 

                for (int i = 0; i < nen; i++){
                    Ke[i][j] = Ke[i][j] + (dshg[i][l]*2.0/h)*(dshg[j][l]*2.0/h)*w[l]*h/2.0;
                }

            }
        }
        // print_Matrix(Ke,nen);
        // print_Vector(Fe,nen);
        // cout << "\n\n";

        for (int j = 0; j < nen; j++){
            // if (j == nen-1) F[n+j] += F[j];
            // else F[n+j] = F[j];
            F[n*(nen-1)+j] += Fe[j];

            for (int i = 0; i < nen; i++){
                // if ((i == nen-1) && (j == nen-1) && (n!=nel-1)) {M[n+i][n+j] += Me[i][j]; cout << "OK\n";}
                // else {M[n+i][n+j] = Me[i][j];}
                K[n*(nen-1)+i][n*(nen-1)+j] += Ke[i][j];
            }
        }
    }
    // print_Matrix(K,np);
    // print_Vector(F,np);


    // Boundary conditions
    // Neumann
    F[0] += BC_g;

    // Dirichlet
    int idx = np - nint;
    F[np-1] = BC_h;
    K[np-1][np-1] = 1.0;
    for(int i = 0; i < nint-1; i++){
        idx += i;
        //cout << idx << endl;
        //cout << K[idx][np-1]<< endl;
        F[idx] += -BC_h*K[idx][np-1]; 
        K[np-1][idx] = 0.0;
        K[idx][np-1] = 0.0;
    }


    // print_Matrix(K,np);
    // print_Vector(F,np);
    
    Eigen::MatrixXd K_eigen(np,np);
    K_eigen = convert_matrix(K, np);
    // std::cout << "Here is the matrix M_eigen:\n" << M_eigen << std::endl;

    Eigen::VectorXd F_eigen(np);
    F_eigen = convert_vector(F, np);
    // std::cout << "Here is the vector F_eigen:\n" << F_eigen << std::endl;

    Eigen::VectorXd u_eigen = K_eigen.colPivHouseholderQr().solve(F_eigen);
    // std::cout << "The solution is:\n" << u_eigen << std::endl;


    // stringstream ss;
    // ss << "data_" << std::setw(4) << std::setfill('0') << nel << ".csv";
    // string FileName = ss.str();

    // ofstream csvFile(FileName);
    // if (!csvFile.is_open()) {
    //     std::cerr << "Error opening the new CSV file." << std::endl;
    //     //return 1; // Return an error code
    // }

    // // Set precision to output all decimals
    // // csvFile << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    // csvFile << std::fixed << std::setprecision(16);  // Adjust precision as needed


    // // Write data to the CSV file
    // for (int i = 0; i < np; ++i) {
    //     csvFile << xl[i];
    //     if (i < np - 1) {
    //         csvFile << ","; // Use a comma as a delimiter
    //     }
    // }
    // csvFile << endl;
    // for (int i = 0; i < np; ++i) {
    //     csvFile << u_eigen(i);
    //     if (i < np - 1) {
    //         csvFile << ","; // Use a comma as a delimiter
    //     }
    // }
    // csvFile << endl;
    // for (int i = 0; i < np; ++i) {
    //     csvFile << u_exact(xl[i]);
    //     if (i < np - 1) {
    //         csvFile << ","; // Use a comma as a delimiter
    //     }
    // }
    // csvFile << "\n"; // Add a newline character to separate rows

    // // Close the CSV file
    // csvFile.close();

    // std::cout << "CSV file written successfully." << std::endl;
    // cout << k << endl;

}


double interpolate_uh(double xx, vector<double> xl, Eigen::VectorXd u_eigen){

    int i = 0;
    while (xx > xl[i]){
        i++;
    }

    double x1, x2, y1, y2;
    x1 = xl[i-1];
    x2 = xl[i];
    cout << " xx = " << xx << " | [" << x1 << "," << x2 << "] == ";

    y1 = u_eigen(i-1);
    y2 = u_eigen(i);
    double uh =  y1 + ((xx-x1) * (y2-y1))/(x2-x1);

    cout << " uh = " << uh << " | [" << y1 << "," << y2 << "]";    

    return uh;
}

int main(){

    int first = 0;

    const int size = 2; 
    vector<double> erros(size-first);

    

    // Initialize the array (optional)
    for (int kk = first; kk < size; ++kk) {

        int nel = pow(2,kk); // number of elements
        cout << "Runing nel = " << nel << "\n";

        // Defining the domain
        double a = 0.0;
        double b = 1.5;

        // boundary conditions
        // Dirichlet
        double BC_h = sin(M_PI*b);
        // Neumann
        double BC_g = -M_PI*cos(M_PI*a);
        cout << "Dirichlet h = " << BC_h << " | Neumann g = " << BC_g << endl;
        
        int k = 4;          // polynomial degree
        int np = k*nel+1;   // mesh total nodes

        int nen = k+1;      // number of element nodes
        int nint = k+1;     // number of integration points

        double h = (b-a)/nel;    // element length

        // cout << "h = " << h << endl;

        vector<double> xl(np,0.0);
        
        xl[0] = a;
        for (int ii = 1; ii < xl.size(); ii++) xl[ii] = xl[ii-1] + h/(nen-1);
        // for (int i = 0; i < xl.size(); i++) cout << xl[i] << "\n";


        // Global K matrix
        vector<vector<double>> K(np, vector<double>(np));
        K = init_Matrix(np,np);
        // print_Matrix(M,np);

        // Global source vector
        vector<double> F(np);
        F = init_F(np);
        // print_Vector(F,np);

        vector<vector<double>> Ke(nen, vector<double>(nen));
        Ke = init_Matrix(nen,nen);

        vector<vector<double>> shg(nint, vector<double>(nint));
        shg = init_Matrix(nint,nint);

        vector<vector<double>> dshg(nint, vector<double>(nint));
        dshg = init_Matrix(nint,nint);

        vector<double> Fe(nen);
        Fe = init_F(nen);

        vector<double> w(nint);
        w = init_F(nint);

        vector<double> pt(nint);
        pt = init_F(nint);

        double xx = 0.0;

        // Global problem construction
        for (int n = 0; n < nel; n++){
            
            Ke = init_Matrix(nen,nen);
            Fe = init_F(nen);

            shg = shl(nen,nint);
            // print_Matrix(shg,nen);
            dshg = dshl(nen,nint);
            // print_Matrix(dshg,nen);
            w = we(nint);
            // print_Vector(w,nint);
            pt = pts(nint);
            // print_Vector(pt,nint);

            // cout << "Element " << n << endl;
            
            for (int l = 0; l < nint; l++){
                
                //xx = h/2*pt[l] + 0.5*(xl[n*(nen-1) + nen-1] + xl[n*(nen-1)]);
                xx = 0.0;

                for (int i = 0; i < nen; i++){
                    xx += shg[i][l]*xl[n*(nen-1) + i];
                }

                for (int j = 0; j < nen; j++){
                    Fe[j] = Fe[j] + f(xx)*shg[j][l]*w[l]*h/2.0; 

                    for (int i = 0; i < nen; i++){
                        Ke[i][j] = Ke[i][j] + (dshg[i][l]*2.0/h)*(dshg[j][l]*2.0/h)*w[l]*h/2.0;
                    }

                }
            }
            // print_Matrix(Ke,nen);
            // print_Vector(Fe,nen);
            // cout << "\n\n";

            for (int j = 0; j < nen; j++){
                // if (j == nen-1) F[n+j] += F[j];
                // else F[n+j] = F[j];
                F[n*(nen-1)+j] += Fe[j];

                for (int i = 0; i < nen; i++){
                    // if ((i == nen-1) && (j == nen-1) && (n!=nel-1)) {M[n+i][n+j] += Me[i][j]; cout << "OK\n";}
                    // else {M[n+i][n+j] = Me[i][j];}
                    K[n*(nen-1)+i][n*(nen-1)+j] += Ke[i][j];
                }
            }
        }
        // print_Matrix(K,np);
        // print_Vector(F,np);


        // Boundary conditions
        // // Neumann
        // F[0] += BC_g;

        // // Dirichlet
        // int idx = np - nint;
        // F[np-1] = BC_h;
        // K[np-1][np-1] = 1.0;
        // for(int i = 0; i < nint-1; i++){
        //     idx += i;
        //     //cout << idx << endl;
        //     //cout << K[idx][np-1]<< endl;
        //     F[idx] += -BC_h*K[idx][np-1]; 
        //     K[np-1][idx] = 0.0;
        //     K[idx][np-1] = 0.0;
        // }

        double kappa_a = 0.0;
        double kappa_b = 1e9;
        double g_a = 0.0;
        double g_b = BC_h;
        double q_a = BC_g;
        double q_b = 0.0;

        K[0][0] += kappa_a;
        K[np-1][np-1] += kappa_b;

        F[0] += kappa_a*g_a + q_a;
        F[np-1] += kappa_b*g_b + q_b;


        // print_Matrix(K,np);
        // print_Vector(F,np);
        
        Eigen::MatrixXd K_eigen(np,np);
        K_eigen = convert_matrix(K, np);
        // std::cout << "Here is the matrix M_eigen:\n" << M_eigen << std::endl;

        Eigen::VectorXd F_eigen(np);
        F_eigen = convert_vector(F, np);
        // std::cout << "Here is the vector F_eigen:\n" << F_eigen << std::endl;

        Eigen::VectorXd u_eigen = K_eigen.colPivHouseholderQr().solve(F_eigen);
        // std::cout << "The solution is:\n" << u_eigen << std::endl;







        stringstream ss;
        ss << "data_" << std::setw(4) << std::setfill('0') << nel << ".csv";
        string FileName = ss.str();

        ofstream csvFile(FileName);
        if (!csvFile.is_open()) {
            std::cerr << "Error opening the new CSV file." << std::endl;
            //return 1; // Return an error code
        }

        // Set precision to output all decimals
        // csvFile << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        csvFile << std::fixed << std::setprecision(16);  // Adjust precision as needed


        // Write data to the CSV file
        for (int i = 0; i < np; ++i) {
            csvFile << xl[i];
            if (i < np - 1) {
                csvFile << ","; // Use a comma as a delimiter
            }
        }
        csvFile << endl;
        for (int i = 0; i < np; ++i) {
            csvFile << u_eigen(i);
            if (i < np - 1) {
                csvFile << ","; // Use a comma as a delimiter
            }
        }
        csvFile << endl;
        for (int i = 0; i < np; ++i) {
            csvFile << u_exact(xl[i]);
            if (i < np - 1) {
                csvFile << ","; // Use a comma as a delimiter
            }
        }
        csvFile << "\n"; // Add a newline character to separate rows

        // Close the CSV file
        csvFile.close();

        std::cout << "CSV file written successfully." << std::endl;
        cout << k << endl;






        // ERRO NORMA L2

        // double erul2 = 0.0;

        // double uh;
        
        // for (int j = 0; j < nel; j++) {
        //     double eru = 0.0;

        //     for (int l = 0; l < nint; l++) {

        //         xx = h/2*pt[l] + 0.5*(xl[j*(nen-1) + nen-1] + xl[j*(nen-1)]);

        //         uh = interpolate_uh(xx, xl, u_eigen);
                
        //         cout << " | u_exact(xx) = " << u_exact(xx) << endl;
                
        //         eru = eru + pow(u_exact(xx) - uh ,2)*w[l]*h/2;
        //     }
        //     cout << endl;

        //     erul2 = erul2 + eru;
        // }
        // erul2 = sqrt(erul2);







        double erul2 = 0.0;

        cout << nen << endl;
        cout << nint << endl;
        cout << nel << endl;
        
        for (int j = 0; j < nel; j++) {
            double eru = 0.0;

            for (int l = 0; l < nint; l++) {

                double uh = 0.0;

                xx = 0.0;

                for (int i = 0; i < nen; i++){
                    uh += shg[i][l]*u_eigen(j*(nen-1)  +i);
                    xx += shg[i][l]*xl[j*(nen-1) + i];
                }
                
                eru = eru + pow(u_exact(xx) - uh ,2)*w[l]*h/2;
            }

            erul2 = erul2 + eru;
        }
        erul2 = sqrt(erul2);





        erros[kk-first] = erul2;


    } 

    // // Set precision for printing
    // std::cout << std::fixed << std::setprecision(16);

    // for (int i = 0; i < size-2; i++){
    //     cout << erros[i] << ", ";
    // } 
    // cout << endl;


    // ESCREVE ERRO 

    stringstream ss;
    ss << "errorL2.csv";
    string FileName2 = ss.str();

    ofstream csvFile(FileName2);
    if (!csvFile.is_open()) {
        std::cerr << "Error opening the new CSV file." << std::endl;
        //return 1; // Return an error code
    }

    // csvFile << std::fixed << std::setprecision(16);

    for (int i = 0; i < size-first; i++) {
        csvFile << erros[i];
        if (i < (size - first - 1)) {
            csvFile << ","; // Use a comma as a delimiter
        } 
    }
    csvFile.close();

    return 0;
}