#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

#include <Eigen/Dense>

#include "shl.cpp"
#include "we.cpp"
#include "pts.cpp"

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
            cout << fixed << setprecision(2) << M[i][j] << "\t";
        }
        cout << endl;
    }
}

void print_Vector(vector<float> F, int dim){
    for (int i = 0; i < dim; i++){
        cout << F[i] << "\n";
    }
}

float f(float xx){
    return sin(M_PI * xx);//*sin(M_PI * xx);
}

Eigen::MatrixXd convert_matrix(vector<vector<float>> M, int dim){
    Eigen::MatrixXd M_eigen(dim,dim);

    for (int i = 0; i < dim; i++){
            for (int j = 0; j < dim; j++){
                M_eigen(i,j) = M[i][j];
            }
    }

    return M_eigen;
}

Eigen::MatrixXd convert_vector(vector<float> F, int dim){
    Eigen::MatrixXd F_eigen(dim,1);

    for (int j = 0; j < dim; j++){
        F_eigen(j,0) = F[j];
    }

    return F_eigen;
}


int main(){

    // Defining the domain
    float a = -2.0;
    float b = 2.0;

    const int size = 9; 
    int numb_el;


    // Initialize the array (optional)
    for (int i = 2; i < size; ++i) {
        numb_el = pow(2,i); 
        cout << "Runing nel = " << numb_el << "\n";
        

        int nel = numb_el;  // number of elements
        
        int k = 2;          // polynomial degree
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
        // print_Matrix(M,np);

        vector<float> F(np);
        F = init_F(np);
        // print_Vector(F,np);

        vector<vector<float>> Me(nen, vector<float>(nen));
        Me = init_Matrix(nen,nen);

        vector<vector<float>> shg(nen, vector<float>(nint));
        shg = init_Matrix(nen,nint);

        vector<float> Fe(nen);
        Fe = init_F(nen);

        vector<float> w(nint);
        w = init_F(nint);

        vector<float> pt(nint);
        pt = init_F(nint);

        float xx = 0.0;

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
                
                xx = h/2*pt[l] + 0.5*(xl[n*(nen-1) + 1] + xl[n*(nen-1)]);

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

            for (int j = 0; j < nen; j++){
                // if (j == nen-1) F[n+j] += F[j];
                // else F[n+j] = F[j];
                F[n*(nen-1)+j] += Fe[j];

                for (int i = 0; i < nen; i++){
                    // if ((i == nen-1) && (j == nen-1) && (n!=nel-1)) {M[n+i][n+j] += Me[i][j]; cout << "OK\n";}
                    // else {M[n+i][n+j] = Me[i][j];}
                    M[n*(nen-1)+i][n*(nen-1)+j] += Me[i][j];
                }
            }
        }
        // print_Matrix(M,np);
        // print_Vector(F,np);

        
        Eigen::MatrixXd M_eigen(np,np);
        M_eigen = convert_matrix(M, np);
        // std::cout << "Here is the matrix M_eigen:\n" << M_eigen << std::endl;

        Eigen::VectorXd F_eigen(np);
        F_eigen = convert_vector(F, np);
        // std::cout << "Here is the vector F_eigen:\n" << F_eigen << std::endl;

        Eigen::VectorXd u_eigen = M_eigen.colPivHouseholderQr().solve(F_eigen);
        // std::cout << "The solution is:\n" << u_eigen << std::endl;


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
        for (int i = 0; i < np; ++i) {
            csvFile << u_eigen(i);
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

    
















    for (int i = 2; i < size; ++i) {

        int numb_el = pow(2,i); 
        int nel = numb_el;  // number of elements
        
        int k = 2;          // polynomial degree
        int np = k*nel+1;   // mesh total nodes

        int nen = k+1;      // number of element nodes
        int nint = k+1;     // number of integration points

        float h = (b-a)/nel;    // element length

        vector<float> w(nint);
        w = init_F(nint);
        w = we(nint);

        stringstream ss;
        ss << "data_" << std::setw(4) << std::setfill('0') << numb_el << ".csv";
        string FileName = ss.str();
        ifstream file(FileName);

        if (!file.is_open()) {
            std::cerr << "Failed to open the CSV file." << std::endl;
            return 1;
        }

        
        std::vector<float> fields_x(pow(2,i)+1);
        std::string line_x;
        std::getline(file, line_x);
        std::stringstream ss_x(line_x);
        std::string field_x;

        int j = 0;
        while (std::getline(ss_x, field_x, ',')) {
            fields_x[j] = std::stof(field_x);
            j++;
        }

        // for (int j = 0; j < pow(2,i)+1; j++) {
        //     std::cout << fields_x[j] << " ";
        // }
        // std::cout << std::endl;







        std::vector<float> fields_u(pow(2,i)+1);
        std::string line_u;
        std::getline(file, line_u);
        std::stringstream ss_u(line_u);
        std::string field_u;

        j = 0;
        while (std::getline(ss_u, field_u, ',')) {
            fields_u[j] = std::stof(field_u);
            j++;
        }


        // for (int j = 0; j < pow(2,i)+1; j++) {
        //     std::cout << fields_u[j] << " ";
        // }
        // std::cout << std::endl;

        






        std::vector<float> fields_uTrue(pow(2,i)+1);
        std::string line_uTrue;
        std::getline(file, line_uTrue);
        std::stringstream ss_uTrue(line_uTrue);
        std::string field_uTrue;

        j = 0;
        while (std::getline(ss_uTrue, field_uTrue, ',')) {
            fields_uTrue[j] = std::stof(field_uTrue);
            j++;
        }

        // for (int j = 0; j < pow(2,i)+1; j++) {
        //     std::cout << fields_uTrue[j] << " ";
        // }
        // std::cout << std::endl;

        file.close();



        float erul2 = 0.0;
        for (int j = 0; j < nel-1; j++) {
            float eru = 0.0;

            for (int l = 0; l < nint-1; l++) {
                // std::cout << "f(xx) - uh = " << fields_uTrue[j] << " - " << fields_u[j] << endl;
                // cout << w[l] << endl;
                eru = eru + pow( fields_uTrue[j+l] - fields_u[j+l] ,2)*w[l]*h/2;
            }

            erul2 = erul2 + eru;
        }
        erul2 = sqrt(erul2);
        cout << "Error L2 norm = " << erul2 << endl;

        
        
    } 

    return 0;
}