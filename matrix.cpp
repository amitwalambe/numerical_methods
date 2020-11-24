/*
Amit Walambe

Program to demonstrate matrix and numerical methods libraries
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>

#include "matrix.h"
#include "nm.h"

using namespace std;

void get_matrix(matrix &m, string s)
{
	cout << endl << s;
    cin >> m;
    cout << m;
}

int condition_number()
{
    matrix A, X1;

	get_matrix(A, "Original matrix");
	get_matrix(A, "Guess Vector");
    
    float max = power(A, X1);
    float min = inv_power(A, X1);

    cout << endl << "Max. eigen value : " << max << endl 
	<< "Min. eigen value : " << min;
    cout << endl << "Condition Number of given matrix is : " 
	<< max / min << endl;
}

/* Function to apply Cholesky method. */
void cholesky()
{
    matrix A;
	get_matrix(A, "Original matrix");

	matrix X = cholesky(A);
    cout << X << endl;
}

/* Function to apply Jacobi method. */
void jacobi()
{
    matrix A;
	get_matrix(A, "Original matrix");
	A = Jacobi(A);
    cout << "After applying Jacobi" << endl << A << endl;
}

/* Function to obtain triangular matrix using Gaussian Elimination method. */
void ge()
{
    matrix A;
	get_matrix(A, "Original matrix");

    gauss_elimination(A);
    cout << "Upper triangular matrix is :" << endl;
    cout << A;
    matrix X = solve_equation(A);
    cout << X;
}

void inverse_power()
{
    matrix A, X;
	get_matrix(A, "Original matrix");
	get_matrix(X, "Guess Vector");
    cout << endl << "Minimum value of lambada : " << inv_power(A, X) << endl;
}

void shifted_power()
{
	matrix A, X;
	get_matrix(A, "Original matrix");
	get_matrix(X, "Guess Vector");
	shifted_power(A, X);
}

/* Function to obtain max. eigen value using power method */
void max_eigen_value()
{
    matrix A, X;
	get_matrix(A, "Original matrix");
	get_matrix(X, "Guess Vector");

    cout << endl << "Maximum value of lambada : " << power(A, X) << endl;
}

/* Function to obtain eigen values using shifted power method */
void shifted_power_function()
{
    matrix A, X;
	get_matrix(A, "Original matrix");
	get_matrix(X, "Guess Vector");

    shifted_power(A, X);
}

int main()
{
    long i;
    int op;
    struct timeval t1, t2;

    fd_set set;
    FD_ZERO(&set);
    FD_SET(fileno(stdin), &set);

    while (1) {
		system("clear");
		printf("1 : Gaussian Elimination\n"
			"2 : Jacobi\n"
			"3 : Cholesky\n"
			"4 : Shifted Power\n"
			"5 : Inverse Power\n"
			"6 : Condition Number\n"
			"7 : Max Eigen Value\n"
			"8 : Exit\n"
			"Select option : ");
		scanf("%d", &op);
		if (op == 8) {
			return 0;
		}
		if (op < 1 || op > 8) {
			printf("Select correct option.");
		} else {
			gettimeofday(&t1, NULL);
			switch (op) {
			case 1:
				ge();
				break;
			case 2:
				jacobi();
				break;
			case 3:
				cholesky();
				break;
			case 4:
				shifted_power();
				break;
			case 5:
				inverse_power();
				break;
			case 6:
				condition_number();
				break;
			case 7:
				max_eigen_value();
				break;
			}
		}
		gettimeofday(&t2, NULL);
		printf("Time for %d operations of given function is : "
			"%ld seconds %ld msec %ld usec.\n",
			op,
			t2.tv_sec - t1.tv_sec,
			(t2.tv_usec - t1.tv_usec) / 1000,
			(t2.tv_usec - t1.tv_usec) % 1000);

		/* wait for user to press any key */
		select(fileno(stdin) + 1, &set, NULL, NULL, &t1);
    }
}
