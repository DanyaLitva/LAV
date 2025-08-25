#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <iomanip>
#include <bitset>
#include <chrono>
using namespace std;


const double T = 0.7;
const int SIMD_Lanes = 4;
const double MinVal = 0.1;
const double MaxVal = 10.0;
double error_rate = 1e-10;
const size_t segment_len = 10;


struct CSRMatrix {
    vector<double> vals;
    vector<int> col_ind;
    vector<int> row_ptr;
};

struct LAVSegment {
    vector<vector<int>> out_order;
    vector<int> chunk_offsets = vector<int>(1, 0);
    vector<vector< bitset<SIMD_Lanes> >> mask;
    vector<vector<vector<double>>> vals;
    vector<vector<vector<int>>> col_id;
};

struct LAVMatrix {
    vector<LAVSegment> segments;

    CSRMatrix sparse_part;
};

vector<vector<double>> create_random_matrix(int rows, int cols, int non_zero) {
    vector<vector<double>> matrix(rows, vector<double>(cols, 0.0));
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis_row(0, rows - 1);
    uniform_int_distribution<> dis_col(0, cols - 1);
    uniform_real_distribution<> dis_val(MinVal, MaxVal);

    int count = 0;
    while (count < non_zero) {
        int i = dis_row(gen);
        int j = dis_col(gen);
        if (matrix[i][j] == 0.0) {
            matrix[i][j] = dis_val(gen);
            count++;
        }
    }
    return matrix;
}

void print_matrix(const vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    cout << "    ";
    for (int j = 0; j < cols; j++) {
        cout << setw(5) << j << " ";
    }
    cout << endl;

    cout << "     ";
    for (int j = 0; j < cols; j++) {
        cout << "------";
    }
    cout << endl;

    for (int i = 0; i < rows; i++) {
        cout << setw(3) << i << " |";
        for (int j = 0; j < cols; j++) {
            cout << fixed << setw(5) << setprecision(2) << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

CSRMatrix dense_to_csr(const vector<vector<double>>& dense) {
    CSRMatrix csr;
    int rows = dense.size();
    int cols = dense[0].size();
    csr.row_ptr.push_back(0);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (dense[i][j] != 0.0) {
                csr.vals.push_back(dense[i][j]);
                csr.col_ind.push_back(j);
            }
        }
        csr.row_ptr.push_back(csr.vals.size());
    }
    return csr;
}

vector<vector<double>> csr_to_dense(const CSRMatrix& csr, int rows, int cols) {
    vector<vector<double>> dense(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        int start = csr.row_ptr[i];
        int end = csr.row_ptr[i + 1];

        for (int j = start; j < end; j++) {
            int col_index = csr.col_ind[j];
            dense[i][col_index] = csr.vals[j];
        }
    }
    return dense;
}

bool matrix_comprasion(vector<vector<double>> left, vector<vector<double>> right) {
    if (left.size() != right.size()) return false;
    if (left.front().size() != right.front().size()) return false;
    for (size_t i = 0; i < left.size(); ++i) {
        for (size_t j = 0; j < left.front().size(); ++j) {
            if (abs(left[i][j] - right[i][j]) > error_rate) return false;
        }
    }

    return true;
}

void print_csr_matrix(const CSRMatrix& csr) {
    cout << "CSR vals: ";
    for (double val : csr.vals) {
        cout << fixed << setprecision(2) << val << " ";
    }
    cout << endl;

    cout << "CSR Col Indices: ";
    for (int col : csr.col_ind) {
        cout << col << " ";
    }
    cout << endl;

    cout << "CSR Row Pointers: ";
    for (int ptr : csr.row_ptr) {
        cout << ptr << " ";
    }
    cout << endl;
}

void print_lav_matrix(const LAVMatrix& lav) {

}

void CSR_to_LAV(CSRMatrix& csr_matrix, LAVMatrix& lav_matrix, int rows, int cols) {
    lav_matrix = LAVMatrix();

    int count_el = csr_matrix.row_ptr.back() - csr_matrix.row_ptr.front();

    vector<pair<int, int>> lav_col(cols);

    for (size_t it = 0; it < csr_matrix.col_ind.size(); ++it) lav_col[csr_matrix.col_ind[it]].second++;
    for (size_t it = 0; it < cols; ++it) lav_col[it].first = it;

    std::sort(lav_col.begin(), lav_col.end(), [](const std::pair<int, int>& left, const std::pair<int, int>& right) {
        return left.second > right.second;
        });

    int separation = 0;
    int temp = 0;
    while ((double(temp) / count_el) < T) {
        temp += lav_col[separation].second;
        separation++;
    }

    //
    //SPARSE PART OF LAV MATRIX
    CSRMatrix sparse_part;
    sparse_part.row_ptr.push_back(0);
    
    vector<int> old_to_sparse(cols, -1);
    int sparse_col_count = cols - separation;
    for (int i = separation; i < cols; i++) {
        old_to_sparse[lav_col[i].first] = i - separation;
    }

    for (int i = 0; i < rows; i++) {
        int start = csr_matrix.row_ptr[i];
        int end = csr_matrix.row_ptr[i + 1];
        int count_in_row = 0;

        for (int j = start; j < end; j++) {
            int old_col = csr_matrix.col_ind[j];
            if (old_to_sparse[old_col] != -1) {
                sparse_part.vals.push_back(csr_matrix.vals[j]);
                sparse_part.col_ind.push_back(old_col);
                count_in_row++;
            }
        }
        sparse_part.row_ptr.push_back(sparse_part.row_ptr.back() + count_in_row);
    }

    lav_matrix.sparse_part = sparse_part;

    //
    //DENSE PART OF LAV MATRIX
    
    size_t segment_count = separation / segment_len;
    if ((separation % segment_len) != 0) segment_count++;

    for (size_t segment_num = 0; segment_num < segment_count; segment_num++) {

    }
}

void MV_dense(vector<vector<double>> dense_matrix, vector<double> x, vector<double>& res, int rows, int cols) {
    for (size_t i = 0; i < rows; ++i) res[i] = 0.0;

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            res[i] += dense_matrix[i][j] * x[j];
        }
    }
}

void SpMV_CSR(const CSRMatrix& mat, const vector<double>& vec, vector<double>& result, int rows, int cols) {
    result.clear();
    result.resize(rows);
    for (int i = 0; i < rows; ++i) {
        double sum = 0.0;
        int start = mat.row_ptr[i];
        int end = mat.row_ptr[i + 1];

        for (int j = start; j < end; ++j) {
            sum += mat.vals[j] * vec[mat.col_ind[j]];
        }
        result[i] = sum;
    }
}

bool vector_comprasion(vector<double> left, vector<double> right) {
    if (left.size() != right.size()) return false;
    for (size_t i = 0; i < left.size(); ++i) {
            if (abs(left[i] - right[i]) > error_rate) return false;
    }

    return true;
}



void SpMV_LAV(const LAVMatrix& lav_matrix, const vector<double>& vec, vector<double>& result, int rows, int cols) {
    result.clear();
    result.resize(rows,0.0);

   


    vector<double> temp_res(rows,0.0);
    SpMV_CSR(lav_matrix.sparse_part, vec, temp_res, rows, cols);
    for (size_t i = 0; i < temp_res.size(); ++i) result[i] += temp_res[i];
}
