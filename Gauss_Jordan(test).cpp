#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>

template <typename scale_type>
class Matrix {
public:
    Matrix() = delete;
    Matrix(std::vector<std::vector<scale_type>> array, int row, int col);

    // get matrix length
    int row() const;
    int col() const;

    // operator[]
    std::vector<scale_type>& operator[](int i);
    const std::vector<scale_type>& operator[](int i) const;

private:
    int _row;
    int _col;
    std::vector<std::vector<scale_type>> _array;
};

template <typename scale_type>
class Gauss_Jordan {
public:
    Gauss_Jordan() = delete;
    Gauss_Jordan(std::vector<std::vector<scale_type>> &M, int row, int col);
    void solve(std::vector<std::vector<scale_type>> &array);
    void print_inverse_matrix() const;

private:
    Matrix<scale_type> _M;
};


// Matrix
template <typename scale_type>
Matrix<scale_type>::Matrix(std::vector<std::vector<scale_type>> array, int row, int col)
                         : _row{row}, _col{col}, _array{array} {}

template <typename scale_type>
int Matrix<scale_type>::row() const { return _row; }

template <typename scale_type>
int Matrix<scale_type>::col() const { return _col; }

template <typename scale_type>
std::vector<scale_type>& Matrix<scale_type>::operator[](int i) {
    assert(i < _row);
    return _array[i];
}

template <typename scale_type>
const std::vector<scale_type>& Matrix<scale_type>::operator[](int i) const {
    assert(i < _row);
    return _array[i];
}

template <typename scale_type>
std::ostream& operator<<(std::ostream &os, const Matrix<scale_type> &mat) {
    int row = mat.row();
    int col = mat.col();
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nThe inverse Matrix is:\n";
    for (int i = 0; i < row; i++) {
        for (int j = row; j < col; j++) {
            os << std::setw(10) << mat[i][j];
        }
        os << std::endl;
    }
    return os;
}

// Gauss_Jordan
template <typename scale_type>
Gauss_Jordan<scale_type>::Gauss_Jordan(std::vector<std::vector<scale_type>> &M, int row, int col) : _M{M, row, col} {};

template <typename scale_type>
void Gauss_Jordan<scale_type>::solve(std::vector<std::vector<scale_type>> &array) {
    int n = _M.row();
    int col = _M.col();
    for (int i = 0; i < n; i++) {

        for (int j = 0; j < n; j++) {
            _M[i][j] = array[i][j];

            if (i == j)
                _M[i][j+n] = (scale_type)1.0;
        }
    }

    for (int i = 0; i < n; i++) {

        // Find pivot element in column i
        int pivot_row = i;
        for (int j = i+1; j < n; j++) {
            if (std::fabs(_M[j][i]) > std::fabs(_M[pivot_row][i]))
                pivot_row = j;
        }

        // Swap current row with pivot row if necessary
        if (pivot_row != i) {
            for (int k = 0; k < col; k++) {
                double tmp = _M[i][k];
                _M[i][k] = _M[pivot_row][k];
                _M[pivot_row][k] = tmp;
            }
        }

        // Scale ith row so that pivot element becomes 1
        double pivot = _M[i][i];
        for (int k = 0; k < col; k++) {
            _M[i][k] = _M[i][k]/pivot;
        }

        // Subtract a multiple of ith row from all other rows to eliminate column i
        for (int j = 0; j < n; j++) {

            if (j != i) {
                double factor = _M[j][i];

                for (int k = 0; k < col; k++) {
                    _M[j][k] = _M[j][k] - factor*_M[i][k];
                }
            }
        }
    }
}

template <typename scale_type>
void Gauss_Jordan<scale_type>::print_inverse_matrix() const {
    std::cout << _M;
}

int main() {
    std::cout << "Input your Matrix(n x n) size\n";
    int n = 0;
    std::cout << "n = ";
    std::cin >> n;
    int row = n;
    int col = n*2;
    std::vector<std::vector<double>> array(n, std::vector<double>(n, 0.0));

    // std::vector<std::vector<double>> array {{1.0, 1.0,  3.0}, 
    //                                         {1.0, 3.0, -3.0},
    //                                         {-2.0, -4.0, -4.0}};                                        

    std::cout << "\nInput your Matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << "A[" << i << "][" << j << "] = ";
            std::cin >> array[i][j];
        }
    }

    std::vector<std::vector<double>> M(row, std::vector<double>(col, 0.0));
    Gauss_Jordan<double> gd(M, row, col);
    gd.solve(array);
    gd.print_inverse_matrix();
    std::cout << "\nDone.\n";
    system("pause");
    return 0;
}




// template <typename scale_type>
// class Gauss_Jordan {
// public:
//     Gauss_Jordan() = delete;
//     Gauss_Jordan(std::vector<std::vector<scale_type>> Ao,
//                  std::vector<std::vector<scale_type>> array,
//                  std::vector<std::vector<scale_type>> Identity, int row, int col)
//                  : _row{row}, _col{col}, _A{Ao, row, col} {

//         // Augment Identity Matrix of Order n to Matrix A           
//         for (int i = 0; i < _row; i++) {
//             for (int j = 0; j < _row; j++) {
//                 _A[i][j] = array[i][j];

//                 if (i == j) {
//                     _A[i][j+_row] = Identity[i][i];
//                 }
//             }
//         }

//         // Apply Gauss Jordan Elimination on Augmented Matrix (A)
//         for (int i = 0; i < _row; i++) {
//             if (_A[i][i] == 0.0) {
//                 std::cout << "Mathematical Error!";
//                 continue;
//             }
//             for (int j = 0; j < _row; j++) {
//                 if (i != j) {
//                     double ratio = _A[j][i]/_A[i][i];
//                     for (int k = 0; k < _col; k++) {
//                         _A[j][k] = _A[j][k] - ratio*_A[i][k];
//                     }
//                 }
//             }
//         }

//         // Row Operation to Convert Principal Diagonal to 1
//          for (int i = 0; i < _row; i++) {
//             for (int j = _row; j < _col; j++) {
//                 _A[i][j] = _A[i][j]/_A[i][i];
//             }
//         }
//     }

//     inline void print_inverse_matrix() const {
//         for (int i = 0; i < _row; i++) {
//             for (int j = _row; j < _col; j++) {
//                 std::cout << _A[i][j] << " ";
//             }
//             std::cout << std::endl;
//         }
//     }

// private:
//     int _row;
//     int _col;
//     Matrix<scale_type> _A;
// };