#include <iostream>
#include <vector>
#include <cmath>

std::vector<double> Gauss_Elimination(std::vector<std::vector<double>> &_A, const int &n);
void print_solution(std::vector<double> &x);

int main() {
    const int n {3};
    std::vector<std::vector<double>> _A {{3.0,  1.0,  2.0, 0.0}, 
                                         {2.0, -3.0, -1.0, 0.0},
                                         {1.0, -2.0,  1.0, 1.0}};
    std::vector<double> x = Gauss_Elimination(_A, n);
    print_solution(x);
    return 0;
}

std::vector<double> Gauss_Elimination(std::vector<std::vector<double>> &_A, const int &n) {

    std::vector<double> x(n, 0.0);

    for (int k = 0; k < n-1; k++) {

        int m = k;
        for (int j = k+1; j < n; j++) {
            if (std::fabs(_A[m][k]) < std::fabs(_A[j][k]))
                m = j;
        }

        if (_A[m][k] == 0.0) {
            std::cout << "No unique solution exists.\n";
            return x;
        }
        else {
            for (int j = 0; j < n+1; j++) {
                double tmp = _A[m][j];
                _A[m][j] = _A[k][j];
                _A[k][j] = tmp;
            }
        }

        if (_A[n-1][n-1] == 0.0) {
            std::cout << "No unique solution exists.\n";
            return x;
        }
        else {
            for (int j = k+1; j < n; j++) {
                double mjk = _A[j][k]/_A[k][k];
                for (int p = k+1; p < n+1; p++) {
                    _A[j][p] = _A[j][p] - mjk*_A[k][p];
                }
            }
        }
    }  

    x[n-1] = _A[n-1][n]/_A[n-1][n-1];
    for (int i = n-2; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i+1; j < n; j++) {
            sum += _A[i][j]*x[j];
        }
        x[i] = (1.0/_A[i][i])*(_A[i][n]-sum);    
    }
    return x;
}

void print_solution(std::vector<double> &x) {
    int j = 0;
    for (double &value : x) {
        std::cout << "x[" << j << "] = " << value << std::endl;
        j++;
    } 
}

// Gauss Jordan Elimination Method from chatGPT
// 1. Let A be the n x n input matrix.
// 2. Set I to the n x n identity matrix.
// 3. Concatenate A and I to form the n x 2n augmented matrix [A|I].
// 4. For k = 1 to n:
//      a. If the kth diagonal entry of [A|I] is zero, swap the kth row with a lower row to make the entry nonzero. If all entries in the kth column below the kth row are zero, terminate the algorithm and report that A is not invertible.
//      b. Divide the kth row of [A|I] by the kth diagonal entry to make the diagonal entry equal to one.
//      c. For each nonzero entry in the kth column of [A|I] below the kth row, subtract a multiple of the kth row from that row to make the entry zero.
// 5. [A|I] is now in row echelon form. For k = n to 1:
//      a. For each nonzero entry in the kth column of [A|I] above the kth row, subtract a multiple of the kth row from that row to make the entry zero.
// 6. [A|I] is now in reduced row echelon form. The right half of [A|I] is the inverse of A.
