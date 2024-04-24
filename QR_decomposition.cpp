
/*
学籍番号: 1211201118
氏名: 林 優花
*/

//Add if necessary
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;
const int n = 5;

template<typename A, size_t N, typename T>
void Fill(A (&array)[N], const T &val){
    std::fill( (T*)array, (T*)(array+N), val );
}

void show_matrix(double matrix[n][n], int N) {
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			printf("%10f",matrix[i][j]);
  	}
		cout << endl;
	}
	cout << endl;
}

double inner_product(double A[n][n], int a, double E[n][n], int e){//Aのa列目とEのe列目の内積を求める
  double sum = 0;
  for(int i = 0; i < n; i++){
		sum += A[i][a] * E[i][e];
  }
  return sum;
}

double size(double A[n][n], int a){//Aのa列目の長さを求める
  double sum = 0;
  for(int i = 0; i < n; i++) {
		sum += A[i][a] * A[i][a];
	}
  return sqrt(sum);
}

double sigma(double matrix_A[n][n], double matrix_q[n][n],int i, int k) {
	double sigma = 0;
	for (int j = 0; j < k; j++) {
		sigma += inner_product(matrix_A, k, matrix_q, j) / inner_product(matrix_q, j, matrix_q, j) * matrix_q[i][j]; 
	}
	//cout << "sigma" << sigma << endl;
	return sigma;
}

int main () {
    
    /* input data for QR decomposition */
    double input_mat[][5] = {
        {2.3, -2.1, 8.0, 5.0, -4.1},
        {2.8, 1.3, 1.0, -9.2, -1.0},
        {1.1, 2.9, 3.1, -2.1, 5.1},
        {-2.1, -6.0, -2.1, 2.0, 3.0},
        {7.0, 0.5, -2.1, 2.0, -5.0},
    };


    
    
    /* print the input_mat before calling gramSchmidt_QRdecomposition */
    std::cout << "input_mat = " << std::endl;
    std::cout << std::endl;
    // write here...
	  show_matrix(input_mat, 5);




    // Gram-Schmidt & QR decomposition
    // write here...
	  int i, j;
	  //matrix_Q
    double matrix_Q[5][5], matrix_q[5][5];
	  for (j = 0; j < n; j++) {
			for (i = 0; i < n; i++) {
				matrix_q[i][j] = input_mat[i][j] - sigma(input_mat, matrix_q, i, j);
		  }
		}	
	  for (j = 0; j < n; j++) {
			for (i = 0; i < n; i++) {
        matrix_Q[i][j] = matrix_q[i][j] / size(matrix_q, j);
		  }
		}
    
    //matrix_R
		double matrix_R[5][5];
	  Fill(matrix_R, 0);
	  for (j = 0; j < n; j++) {
			for (i = 0; i < j+1; i++) {
				matrix_R[i][j] = inner_product(input_mat, j, matrix_Q, i);
			}	
		}		

	  //matrix_QR
	  double matrix_QR[5][5], matrix_A_QR[5][5], term;
		for (i = 0; i < n; i++) {
	    for (j = 0; j < n; j++) {
			  matrix_QR[i][j] = matrix_Q[i][j] * matrix_R[i][j];
				for (i = 0; i < n; i++) {
	        for (j = 0; j < n; j++) {
						term = 0;
						for(int k = 0; k < n; k++) {
              term = term + matrix_Q[i][k] * matrix_R[k][j];
						}
            matrix_QR[i][j] = term;
				  }		
			  }
		  }
		}
		for (i = 0; i < n; i++) {
	    for (j = 0; j < n; j++) {
				matrix_A_QR[i][j] = input_mat[i][j] - matrix_QR[i][j];
			}
		}
	
    /* print the matrix Q resulting from gramSchmidt_QRdecomposition */
    std::cout << "Q = " << std::endl;
    std::cout << std::endl;
    // write here...
    show_matrix(matrix_Q, 5);




    /* print the matrix R resulting from gramSchmidt_QRdecomposition */
    std::cout << "R = " << std::endl;
    std::cout << std::endl;
    // write here...
    show_matrix(matrix_R, 5);



    
    /* print the matrix Q * R */
    std::cout << "Q * R = ?" << std::endl;
    std::cout << std::endl;
    // write here...
    show_matrix(matrix_QR, 5);




    
    /* print Input_mat - (Q * R) */
    std::cout << "Input_mat - (Q * R) = ?" << std::endl;
    std::cout << std::endl;
    // write here...
    show_matrix(matrix_A_QR, 5);
    
    


    return 0;
}
