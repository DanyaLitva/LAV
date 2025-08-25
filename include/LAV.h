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
    lav_matrix.segments.resize(segment_count);

    vector<int> old_to_new(cols, -1);
    for (int i = 0; i < cols; i++) {
        old_to_new[lav_col[i].first] = i;
    }

    for (size_t segment_num = 0; segment_num < segment_count; segment_num++) {
        size_t start_seg = segment_num * segment_len;
        size_t end_seg = start_seg + segment_len;
        if (segment_num == (segment_count - 1)) end_seg = separation;

        vector<pair<int, int>> row_nonzeros(rows);
        for (int i = 0; i < rows; i++) {
            row_nonzeros[i] = { i, 0 };

            int start = csr_matrix.row_ptr[i];
            int end = csr_matrix.row_ptr[i + 1];
            for (int j = start; j < end; j++) {
                int old_col = csr_matrix.col_ind[j];
                int new_col = old_to_new[old_col];
                if ((start_seg<= new_col) && (new_col < end_seg)) {
                    row_nonzeros[i].second++;
                }
            }
        }

        sort(row_nonzeros.begin(), row_nonzeros.end(),
            [](const pair<int, int>& a, const pair<int, int>& b) {
                return a.second > b.second;
            });

        vector<int> old_row_to_new(rows);
        for (int i = 0; i < rows; i++) {
            old_row_to_new[row_nonzeros[i].first] = i;
        }

        std::vector<std::vector<std::vector<double>>> chunks(
            int(std::ceil(double(rows) / SIMD_Lanes)),
            std::vector<std::vector<double>>(
                SIMD_Lanes,
                std::vector<double>(end_seg - start_seg, 0.0)
            )
        );

        vector<int> max_len_in_chunk;
        for (size_t i = 0; i < rows; ++i) {
            if (i % SIMD_Lanes == 0) max_len_in_chunk.push_back(row_nonzeros[i].second);
            else {
                if (max_len_in_chunk.back() < row_nonzeros[i].second) max_len_in_chunk.back() = row_nonzeros[i].second;

            }
        }

        int ind = 0;
        for (int i : max_len_in_chunk) {
            lav_matrix.segments[segment_num].chunk_offsets.push_back(lav_matrix.segments[segment_num].chunk_offsets.back() + i);
            lav_matrix.segments[segment_num].vals.push_back(vector<vector<double>>());
            lav_matrix.segments[segment_num].col_id.push_back(vector<vector<int>>());
            for (int j = 0; j < SIMD_Lanes; ++j) {
                lav_matrix.segments[segment_num].vals[ind].push_back(vector<double>(i));
                lav_matrix.segments[segment_num].col_id[ind].push_back(vector<int>(i, -1));
            }
            ++ind;
        }

        if (rows % SIMD_Lanes != 0) {
            lav_matrix.segments[segment_num].vals.back() = vector<vector<double>>(rows % SIMD_Lanes, vector<double>(max_len_in_chunk.back(), 0.0));
            lav_matrix.segments[segment_num].col_id.back() = vector<vector<int>>(rows % SIMD_Lanes, vector<int>(max_len_in_chunk.back(), -1));
        }

        vector<pair<double, int>> actual_lane_vals;
        vector<pair<int, int>> actual_lane_cols;

        for (int i = 0; i < rows; i++) {
            int start = csr_matrix.row_ptr[i];
            int end = csr_matrix.row_ptr[i + 1];
            int new_row = old_row_to_new[i];

            actual_lane_vals.clear();
            actual_lane_cols.clear();
            for (int j = start; j < end; j++) {
                int old_col = csr_matrix.col_ind[j];
                int new_col = old_to_new[old_col];
                if ((start_seg <= new_col) && (new_col < end_seg)) {
                    int chunk_index = new_row / SIMD_Lanes;
                    int row_in_chunk = new_row % SIMD_Lanes;
                    chunks[chunk_index][row_in_chunk][new_col] = csr_matrix.vals[j];
                    actual_lane_vals.push_back(make_pair(csr_matrix.vals[j], new_col));
                    actual_lane_cols.push_back(make_pair(old_col, new_col));
                }
            }

            int temp = 0;
            sort(actual_lane_vals.begin(), actual_lane_vals.end(),
                [](const pair<double, int>& a, const pair<double, int>& b) {
                    return a.second < b.second;
                });

            sort(actual_lane_cols.begin(), actual_lane_cols.end(),
                [](const pair<int, int>& a, const pair<int, int>& b) {
                    return a.second < b.second;
                });
            for (int j = 0; j < actual_lane_vals.size(); ++j) {
                lav_matrix.segments[segment_num].vals[new_row / SIMD_Lanes][new_row % SIMD_Lanes][j] = actual_lane_vals[j].first;
                lav_matrix.segments[segment_num].col_id[new_row / SIMD_Lanes][new_row % SIMD_Lanes][j] = actual_lane_cols[j].first;
            }
        }


        vector<vector<int>> out_order(int(std::ceil(double(rows) / SIMD_Lanes)));
        for (size_t i = 0; i < rows; ++i) {
            out_order[i / SIMD_Lanes].push_back(row_nonzeros[i].first);
        }

        lav_matrix.segments[segment_num].out_order = out_order;
        lav_matrix.segments[segment_num].mask.resize(lav_matrix.segments[segment_num].chunk_offsets.size() - 1);
        bitset<SIMD_Lanes> actual_mask;
        for (size_t i = 0; i < lav_matrix.segments[segment_num].chunk_offsets.size() - 1; ++i) {
            for (size_t j = 0; j < lav_matrix.segments[segment_num].chunk_offsets[i + 1] - lav_matrix.segments[segment_num].chunk_offsets[i]; ++j) {
                actual_mask = 0;
                for (size_t k = 0; k < lav_matrix.segments[segment_num].col_id[i].size(); ++k) {
                    if (lav_matrix.segments[segment_num].col_id[i][k][j] != -1) actual_mask |= (1 << k);
                }
                lav_matrix.segments[segment_num].mask[i].push_back(actual_mask);
            }
        }


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
    result.resize(rows, 0.0);

    for (size_t s = 0; s < lav_matrix.segments.size(); ++s) {
        const LAVSegment& seg = lav_matrix.segments[s];
        const size_t seg_nChunks = seg.vals.size();
        if (seg_nChunks == 0) continue;

        for (size_t c = 0; c < seg_nChunks; ++c) {
            int pos_begin = seg.chunk_offsets[c];
            int pos_end = seg.chunk_offsets[c + 1];
            int num_pos = pos_end - pos_begin;
            if (num_pos <= 0) continue;

            int num_lanes = static_cast<int>(seg.vals[c].size());
            if (num_lanes <= 0) continue;

            uint32_t rowmask_bits;
            if (SIMD_Lanes >= 32) rowmask_bits = 0xFFFFFFFFu;
            else rowmask_bits = ((1u << SIMD_Lanes) - 1u);

            if (c == seg_nChunks - 1) {
                int rows_in_last = rows - static_cast<int>(c) * SIMD_Lanes;
                if (rows_in_last < 0) rows_in_last = 0;
                if (rows_in_last < SIMD_Lanes) {
                    int shift = SIMD_Lanes - rows_in_last;
                    if (shift >= 0 && shift < SIMD_Lanes) rowmask_bits = rowmask_bits >> shift;
                    else if (rows_in_last == 0) rowmask_bits = 0;
                }
            }

            std::bitset<SIMD_Lanes> rowmask_bitset;
            for (int b = 0; b < SIMD_Lanes; ++b) {
                if (rowmask_bits & (1u << b)) rowmask_bitset.set(b);
                else rowmask_bitset.reset(b);
            }

            for (int p = 0; p < num_pos; ++p) {
                std::bitset<SIMD_Lanes> active = seg.mask[c][p] & rowmask_bitset;
                if (active.none()) continue;

                int common_col = seg.col_id[c][0][p];
                bool can_use_common = (common_col >= 0);
                if (can_use_common) {
                    for (int lane = 0; lane < num_lanes; ++lane) {
                        if (!active.test(lane)) continue;
                        if (seg.col_id[c][lane][p] != common_col) { can_use_common = false; break; }
                    }
                }

                double xval_common = 0.0;
                if (can_use_common) {
                    if (static_cast<size_t>(common_col) < vec.size()) xval_common = vec[common_col];
                    else can_use_common = false;
                }

                for (int lane = 0; lane < num_lanes; ++lane) {
                    if (!active.test(lane)) continue;
                    int col = seg.col_id[c][lane][p];
                    if (col < 0) continue;
                    if (static_cast<size_t>(col) >= vec.size()) continue;
                    double xval = can_use_common ? xval_common : vec[col];
                    int row = seg.out_order[c][lane];
                    if (row < 0 || static_cast<size_t>(row) >= result.size()) continue;
                    result[row] += seg.vals[c][lane][p] * xval;
                }
            }
        }
    }

    vector<double> temp_res(rows, 0.0);
    SpMV_CSR(lav_matrix.sparse_part, vec, temp_res, rows, cols);
    for (int i = 0; i < rows; ++i) result[i] += temp_res[i];
}