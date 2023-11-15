#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <cfloat>

#include "shl_v2.cpp"
#include "we_v2.cpp"
#include "pts_v2.cpp"

using namespace std;

vector<double> init_F(int n){
    vector<double> F(n);
    for (int i = 0; i < n; i++) F[i] = 0.0;
    return F;
}

int main(){
    const int size = 6; 

    vector<double> erros(size-2);
    vector<double> derros(size-2);


    for (int i = 2; i < size; ++i) {

        double a = 0.0;
        double b = 1.5;

        int numb_el = pow(2,i); 
        int nel = numb_el;  // number of elements
        
        int k = 1;          // polynomial degree
        int np = k*nel+1;   // mesh total nodes

        int nen = k+1;      // number of element nodes
        int nint = k+1;     // number of integration points

        double h = (b-a)/nel;    // element length

        vector<double> w(nint);
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

        
        std::vector<double> fields_x(np);
        std::string line_x;
        std::getline(file, line_x);
        std::stringstream ss_x(line_x);
        std::string field_x;

        int j = 0;
        while (std::getline(ss_x, field_x, ',')) {
            fields_x[j] = std::stof(field_x);
            j++;
        }

        // for (int j = 0; j < np; j++) {
        //     std::cout << fields_x[j] << " ";
        // }
        // std::cout << std::endl;







        std::vector<double> fields_u(np);
        std::string line_u;
        std::getline(file, line_u);
        std::stringstream ss_u(line_u);
        std::string field_u;

        j = 0;
        while (std::getline(ss_u, field_u, ',')) {
            fields_u[j] = std::stof(field_u);
            j++;
        }


        // for (int j = 0; j < np; j++) {
        //     std::cout << fields_u[j] << " ";
        // }
        // std::cout << std::endl;

        






        std::vector<double> fields_uTrue(np);
        std::string line_uTrue;
        std::getline(file, line_uTrue);
        std::stringstream ss_uTrue(line_uTrue);
        std::string field_uTrue;

        j = 0;
        while (std::getline(ss_uTrue, field_uTrue, ',')) {
            fields_uTrue[j] = std::stof(field_uTrue);
            j++;
        }

        // for (int j = 0; j < np; j++) {
        //     std::cout << fields_uTrue[j] << " ";
        // }
        // std::cout << std::endl;

        file.close();


        // ERRO NORMA L2

        double erul2 = 0.0;
        for (int j = 0; j < nel; j++) {
            double eru = 0.0;

            for (int l = 0; l < nint-1; l++) {
                // std::cout << "f(xx) - uh = " << fields_uTrue[j] << " - " << fields_u[j] << endl;
                // cout << w[l] << endl;
                eru = eru + pow( fields_uTrue[j*(nen-1)+l] - fields_u[j*(nen-1)+l] ,2)*w[l]*h/2;
            }

            erul2 = erul2 + eru;
        }
        erul2 = sqrt(erul2);



        // // ERRO NORMA INFTY

        // double erul2 = 0.0;
        // double eru = 0.0;

        // for (int j = 0; j < np; j++) {
            
        //     eru = abs(fields_uTrue[j] - fields_u[j]);
        //     if(eru > erul2){
        //         erul2 = eru;
        //     }
            
        // }


        // // ERRO NORMA L2 - DERIVADA

        // double derul2 = 0.0;
        // for (int j = 0; j < nel; j++) {
        //     double deru = 0.0;

        //     for (int l = 0; l < nint-1; l++) {
        //         // std::cout << "f(xx) - uh = " << fields_uTrue[j] << " - " << fields_u[j] << endl;
        //         // cout << w[l] << endl;
        //         deru = deru + pow( fields_uTrue[j*(nen-1)+l] - fields_u[j*(nen-1)+l] ,2)*w[l]*h/2;
        //     }

        //     derul2 = derul2 + deru;
        // }
        // derul2 = sqrt(derul2);
        


        // SALVAR ERRO
        erros[i-2] = erul2;
        // derros[i-2] = derul2;
        
    } 

    // ESCREVE ERRO 

    stringstream ss;
    ss << "errorL2.csv";
    string FileName2 = ss.str();

    ofstream csvFile(FileName2);
    if (!csvFile.is_open()) {
        std::cerr << "Error opening the new CSV file." << std::endl;
        //return 1; // Return an error code
    }

    for (int i = 0; i < size-2; i++) {
        csvFile << erros[i];
        if (i < size - 3) {
            csvFile << ","; // Use a comma as a delimiter
        } 
    }
    csvFile.close();

    // ESCREVE ERRO DERIVADA
    
    // stringstream ss;
    // ss << "d_errorL2.csv";
    // string FileName2 = ss.str();

    // ofstream csvFile(FileName2);
    // if (!csvFile.is_open()) {
    //     std::cerr << "Error opening the new CSV file." << std::endl;
    //     //return 1; // Return an error code
    // }

    // for (int i = 0; i < size-2; i++) {
    //     csvFile << erros[i];
    //     if (i < size - 3) {
    //         csvFile << ","; // Use a comma as a delimiter
    //     } 
    // }
    // csvFile.close();

    return 0;
}