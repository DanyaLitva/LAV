//#include "LAV.h"
#include "LAV_without_segmentation.h"
using namespace std;



int main()
{
    int rows_count = 8;
    int cols_count = 8;
    int count_el_in_matrix = 19;
    
    
    
    
    int rows = rows_count;
    int cols = cols_count;
    int count_el = count_el_in_matrix;
    auto dense_matrix = create_random_matrix(rows, cols, count_el);

    dense_matrix = {
        { 1,0,0,2,0,3,0,0 },
        {0,4,0,5,0,0,0,0},
        {6,0,7,0,8,0,0,0},
        {0,9,0,0,0,0,0,0},
        {0,10,0,11,0,12,0,0},
        {0,13,0,14,0,0,15,0},
        {0,0,16,0,0,0,0,17},
        {0,18,0,0,0,19,0,0}
    };

    cout << "Original Matrix:" << endl;
    print_matrix(dense_matrix);

    CSRMatrix csr_matrix = dense_to_csr(dense_matrix);

    cout << "\nCSR Matrix:" << endl;

    print_csr_matrix(csr_matrix);
    /*
    auto reconstructed_matrix = csr_to_dense(csr_matrix, rows, cols);

    cout << "\nReconstructed Matrix:" << endl;
    print_matrix(reconstructed_matrix);

    cout << endl << "Correct? ";
    if (matrix_comprasion(dense_matrix, reconstructed_matrix)) cout << "yes"; else cout << "no"; cout << endl;*/

    //cout << endl << "CSR to LAV convert" << endl;
    LAVMatrix lav_matrix;

    CSR_to_LAV(csr_matrix, lav_matrix, rows, cols);

    //print_lav_matrix(lav_matrix);
    //print_lav_matrix_w_letters_and_shift(lav_matrix);


    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    


    //start = std::chrono::steady_clock::now();
    //end = std::chrono::steady_clock::now();
    //elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    //cout << "last it time: " << elapsed / 1000. << " seconds" << endl;

    cout << "Mult: \n";

    cols = rows = 10001;
    dense_matrix = create_random_matrix(cols, rows,1000000);


    vector<double> x(cols);
    vector<double> y(rows,0.0);
    vector<double> temp_y(rows);

    x = vector<double>(cols, 1);

    

    cout << "dense: \n\t";
    start = std::chrono::steady_clock::now();
    MV_dense(dense_matrix, x, y, rows, cols);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "it time: " << elapsed << " miliseconds" << endl;

    /*for (double i : y) cout << i << " ";
    cout << endl;*/
    temp_y = y;

    start = std::chrono::steady_clock::now();
    csr_matrix = dense_to_csr(dense_matrix);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "\ndense to csr time: " << elapsed << " miliseconds" << endl << endl;


    cout << "CSR: \n\t";
    start = std::chrono::steady_clock::now();
    SpMV_CSR(csr_matrix, x, y, rows, cols);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "it time: " << elapsed << " miliseconds" << endl;
    /*for (double i : y) cout << i << " ";
    cout << endl;*/

    cout << "\tIs correct?: ";
    if (vector_comprasion(y, temp_y)) cout << "yes";
    else cout << "no";
    cout << endl;


    start = std::chrono::steady_clock::now();
    CSR_to_LAV(csr_matrix, lav_matrix, rows, cols);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "\ncsr to lav time: " << elapsed << " miliseconds" << endl<<endl;

    cout << "LAV: \n\t";
    start = std::chrono::steady_clock::now();
    SpMV_LAV(lav_matrix, x, y, rows, cols);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "it time: " << elapsed << " miliseconds" << endl;
    /*for (double i : y) cout << i << " ";
    cout << endl;*/

    cout << "\tIs correct?: ";
    if (vector_comprasion(y, temp_y)) cout << "yes";
    else cout << "no";
    cout << endl;
}






