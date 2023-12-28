#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

// #include <Eigen/Dense>

#include "shl_v2.cpp"
#include "shl_wh.cpp"
#include "dshl_v2.cpp"
#include "dshl2_v2.cpp"
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
    return 0.0;
}

double u_exact(double xx, double tt, double kappa, double epsilon){
    double A = 1.0 / (sqrt(4*tt + 1.0));
    double B = pow(xx - kappa*tt - 0.5,2);
    double C = epsilon*(4*tt + 1.0);
    return A * exp(-B/C);
}

double K_func(double xx, double epsilon){
    return epsilon;
}

double gamma_func(double xx){
    return 0.0;
}

// Eigen::MatrixXd convert_matrix(vector<vector<double>> M, int dim){
//     Eigen::MatrixXd M_eigen(dim,dim);

//     for (int i = 0; i < dim; i++){
//             for (int j = 0; j < dim; j++){
//                 M_eigen(i,j) = M[i][j];
//             }
//     }

//     return M_eigen;
// }

// Eigen::MatrixXd convert_vector(vector<double> F, int dim){
//     Eigen::MatrixXd F_eigen(dim,1);

//     for (int j = 0; j < dim; j++){
//         F_eigen(j,0) = F[j];
//     }

//     return F_eigen;
// }


// Function to perform Gaussian elimination and back substitution
void solveLinearSystem(vector<vector<double>>& A, vector<double>& b, vector<double>& x) {
    int n = A.size();

    // Forward elimination
    for (int i = 0; i < n; ++i) {
        // Make the diagonal element 1
        double divisor = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= divisor;
        }
        b[i] /= divisor;

        // Eliminate other elements in the current column
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
    }

    // Back substitution
    for (int i = 0; i < n; ++i) {
        x[i] = b[i];
    }
}

// Function to multiply a matrix and a vector
std::vector<double> matrixVectorMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector) {
    // Check if the matrix and vector dimensions are compatible for multiplication
    size_t numRows = matrix.size();
    size_t numCols = matrix[0].size();
    size_t vectorSize = vector.size();

    if (numCols != vectorSize) {
        // Dimensions are not compatible for multiplication
        throw std::invalid_argument("Matrix and vector dimensions are not compatible for multiplication");
    }

    // Resulting vector
    std::vector<double> result(numRows, 0.0);

    // Perform matrix-vector multiplication
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

// Function to add two matrices
std::vector<std::vector<double>> addMatrices(const std::vector<std::vector<double>>& matrix1, const std::vector<std::vector<double>>& matrix2) {
    // Check if the matrices have the same dimensions
    size_t numRows1 = matrix1.size();
    size_t numCols1 = matrix1[0].size();
    size_t numRows2 = matrix2.size();
    size_t numCols2 = matrix2[0].size();

    if (numRows1 != numRows2 || numCols1 != numCols2) {
        // Matrices have different dimensions, cannot perform addition
        throw std::invalid_argument("Matrices have different dimensions, cannot perform addition");
    }

    // Resulting matrix
    std::vector<std::vector<double>> result(numRows1, std::vector<double>(numCols1, 0.0));

    // Perform matrix addition
    for (size_t i = 0; i < numRows1; ++i) {
        for (size_t j = 0; j < numCols1; ++j) {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }

    return result;
}

// Function to add two vectors
std::vector<double> addVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    // Check if the vectors have the same size
    size_t size1 = vector1.size();
    size_t size2 = vector2.size();

    if (size1 != size2) {
        // Vectors have different sizes, cannot perform addition
        throw std::invalid_argument("Vectors have different sizes, cannot perform addition");
    }

    // Resulting vector
    std::vector<double> result(size1, 0.0);

    // Perform vector addition
    for (size_t i = 0; i < size1; ++i) {
        result[i] = vector1[i] + vector2[i];
    }

    return result;
}

bool isMultiple(double t, double dt_save, double tolerance = 1e-4) {
    // Check if the absolute difference is within the specified tolerance
    return std::abs(t - dt_save * std::round(t / dt_save)) < tolerance;
}

int main(){

    double epsilon = 1e-2;
    cout << "epsilon = " << epsilon << endl;

    double kappa = 1.0;
    cout << "kappa = " << kappa << endl;

    

    int nel = 100;
    cout << "Runing nel = " << nel << "\n";
    // Defining the domain
    double a = 0.0;
    double b = 2.0;
    double h = (b-a)/nel;    // element length
    cout << "h = " << h << endl;
 
    double Pe_h = abs(kappa)*h/(2*epsilon);
    cout << "Pe_h = " << Pe_h << endl;

    double beta = max(0.0,1.0 - 1.0/Pe_h); // beta >= este valor
    cout << "beta = " << beta << endl;

    double tau = beta*h/(2.0*abs(kappa));


    double t = 0.0;
    double T = 3.0;
    double dt = h*h;
    double dt_save = 0.25;
    cout << "dt = " << dt << endl;


    double CFL = kappa*dt/h;
    cout << "CFL = " << CFL << endl;


    int k = 1;          // polynomial degree
    int np = k*nel+1;   // mesh total nodes
    int nen = k+1;      // number of element nodes
    int nint = k+1;     // number of integration points



    // Initializing the mesh and matrices
    vector<double> xl(np,0.0);
        
    xl[0] = a;
    for (int ii = 1; ii < xl.size(); ii++) xl[ii] = xl[ii-1] + h/(nen-1);
    // for (int i = 0; i < xl.size(); i++) cout << xl[i] << "\n";


    // Global K matrix
    vector<vector<double>> K(np, vector<double>(np));
    K = init_Matrix(np,np);
    // print_Matrix(K,np);

    // Global C matrix
    vector<vector<double>> C(np, vector<double>(np));
    C = init_Matrix(np,np);
    // print_Matrix(C,np);

    // Global M matrix
    vector<vector<double>> M(np, vector<double>(np));
    M = init_Matrix(np,np);
    // print_Matrix(M,np);

    // Global source vector
    vector<double> F(np);
    F = init_F(np);
    // print_Vector(F,np);

    // Solution vector
    vector<double> u_solution_n(np);
    u_solution_n = init_F(np);
    for (int ii = 1; ii < u_solution_n.size(); ii++) u_solution_n[ii] = u_exact(xl[ii], 0.0, kappa, epsilon);  
    // for (int ii = 0; ii < u_solution_n.size(); ii++) cout << u_solution_n[ii] << ",";
    // cout << endl;

    vector<double> u_solution_np1(np);
    u_solution_np1 = init_F(np);

    vector<double> un_M(np);
    un_M = init_F(np);



    // Matrizes e vetores locais (elemento)
    vector<vector<double>> Ke(nen, vector<double>(nen));
    Ke = init_Matrix(nen,nen);

    vector<vector<double>> Me(nen, vector<double>(nen));
    Me = init_Matrix(nen,nen);

    vector<vector<double>> Ce(nen, vector<double>(nen));
    Ce = init_Matrix(nen,nen);

    vector<double> Fe(nen);
    Fe = init_F(nen);

    
    



    // Funçoes de base e integração numérica
    vector<vector<double>> shg(nint, vector<double>(nint));
    shg = init_Matrix(nint,nint);

    vector<vector<double>> shg_wh(nint, vector<double>(nint));
    shg_wh = init_Matrix(nint,nint);

    vector<vector<double>> dshg(nint, vector<double>(nint));
    dshg = init_Matrix(nint,nint);

    vector<vector<double>> dshg2(nint, vector<double>(nint));
    dshg2 = init_Matrix(nint,nint);

    vector<double> w(nint);
    w = init_F(nint);

    vector<double> pt(nint);
    pt = init_F(nint);
        



    vector<double> erros;

 
    while (t <= T) { 
        // cout << "t = " << t << endl;

        double xx = 0.0;

        M = init_Matrix(np,np);
        K = init_Matrix(np,np);
        C = init_Matrix(np,np);
        F = init_F(np);

        // Global problem construction
        for (int n = 0; n < nel; n++){
            
            Ke = init_Matrix(nen,nen);
            Me = init_Matrix(nen,nen);
            Ce = init_Matrix(nen,nen);
            Fe = init_F(nen);

            shg = shl(nen,nint);
            // print_Matrix(shg,nen);
            shg_wh = shl_wh(nen,nint,beta);
            // print_Matrix(shg_wh,nen);
            dshg = dshl(nen,nint);
            // print_Matrix(dshg,nen);
            dshg2 = dshl2(nen,nint);
            // print_Matrix(dshg2,nen);
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
                    Fe[j] = Fe[j] + f(xx)*shg_wh[j][l]*w[l]*h/2.0; 

                    for (int i = 0; i < nen; i++){
                        Ke[i][j] = Ke[i][j] + K_func(xx, epsilon)*(dshg[i][l]*2.0/h)*(dshg[j][l]*2.0/h)*w[l]*h/2.0
                                   + gamma_func(xx)*shg[i][l]*shg_wh[j][l]*w[l]*h/2;
                    
                        Me[i][j] = Me[i][j] + (1.0/dt)*shg[i][l]*shg[j][l]*w[l]*h/2;

                        Ce[i][j] = Ce[i][j] + kappa*(dshg[j][l]*2.0/h)*(shg_wh[i][l])*w[l]*h/2.0
                                   - K_func(xx, epsilon)*(dshg2[i][l]*2.0/h)*kappa*tau*(dshg[j][l]*2.0/h)*w[l]*h/2.0;
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
                    C[n*(nen-1)+i][n*(nen-1)+j] += Ce[i][j];
                    M[n*(nen-1)+i][n*(nen-1)+j] += Me[i][j];
                }
            }
        }

        un_M = matrixVectorMultiply(M,u_solution_n);
        F = addVectors(F,un_M);

        K = addMatrices(K,M);
        K = addMatrices(K,C);

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

        double kappa_a = 1e9;
        double kappa_b = 1e9;
        double g_a = u_exact(0.0, t, kappa, epsilon);
        double g_b = u_exact(2.0, t, kappa, epsilon);
        double q_a = 0.0;
        double q_b = 0.0;

        K[0][0] += kappa_a;
        K[np-1][np-1] += kappa_b;

        F[0] += kappa_a*g_a + q_a;
        F[np-1] += kappa_b*g_b + q_b;


        // print_Matrix(K,np);
        // print_Vector(F,np);
        
        // Eigen::MatrixXd K_eigen(np,np);
        // K_eigen = convert_matrix(K, np);
        // // std::cout << "Here is the matrix M_eigen:\n" << M_eigen << std::endl;

        // Eigen::VectorXd F_eigen(np);
        // F_eigen = convert_vector(F, np);
        // // std::cout << "Here is the vector F_eigen:\n" << F_eigen << std::endl;

        // Eigen::VectorXd u_eigen = K_eigen.colPivHouseholderQr().solve(F_eigen);
        // // std::cout << "The solution is:\n" << u_eigen << std::endl;

        solveLinearSystem(K, F, u_solution_np1);


        if (isMultiple(t,dt_save)){

            cout << "t = " << t << endl;

            stringstream ss;
            ss << "data_" << std::setw(4) << std::setfill('0') << nel << "t_" << std::fixed << std::setprecision(4) << t <<".csv";
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
                csvFile << u_solution_n[i];
                if (i < np - 1) {
                    csvFile << ","; // Use a comma as a delimiter
                }
            }
            csvFile << endl;
            for (int i = 0; i < np; ++i) {
                csvFile << u_exact(xl[i], t, kappa, epsilon);
                if (i < np - 1) {
                    csvFile << ","; // Use a comma as a delimiter
                }
            }
            csvFile << "\n"; // Add a newline character to separate rows

            // Close the CSV file
            csvFile.close();

            std::cout << "CSV file written successfully." << std::endl;
            cout << k << endl;
        
        }





        // // ERRO NORMA L2

        // double erul2 = 0.0;

        // cout << nen << endl;
        // cout << nint << endl;
        // cout << nel << endl;
        
        // for (int j = 0; j < nel; j++) {
        //     double eru = 0.0;

        //     for (int l = 0; l < nint; l++) {

        //         double uh = 0.0;

        //         xx = 0.0;

        //         for (int i = 0; i < nen; i++){
        //             uh += shg[i][l]*u_solution[j*(nen-1) + i];
        //             xx += shg[i][l]*xl[j*(nen-1) + i];
        //         }
                
        //         eru = eru + pow(u_exact(xx, epsilon, kappa) - uh ,2)*w[l]*h/2;
        //     }

        //     erul2 = erul2 + eru;
        // }
        // erul2 = sqrt(erul2);



        for (int ii = 0; ii < u_solution_n.size(); ii++) u_solution_n[ii] = u_solution_np1[ii];

        //erros[kk-first] = erul2;
        t += dt;

    } 



    // // ESCREVE ERRO 

    // stringstream ss;
    // ss << "errorL2.csv";
    // string FileName2 = ss.str();

    // ofstream csvFile(FileName2);
    // if (!csvFile.is_open()) {
    //     std::cerr << "Error opening the new CSV file." << std::endl;
    //     //return 1; // Return an error code
    // }

    // csvFile << std::fixed << std::setprecision(16);

    // for (int i = 0; i < size-first; i++) {
    //     csvFile << erros[i];
    //     if (i < (size - first - 1)) {
    //         csvFile << ","; // Use a comma as a delimiter
    //     } 
    // }
    // csvFile.close();

    return 0;
}