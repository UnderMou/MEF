

int main(){
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
}