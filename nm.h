/*
Amit Walambe

functionality :
	Find L by Cholesky method (Applicable for symmetric matrix)
	Gaussian elimination method
	inverse power method
	jacobi method
	power method
	solve_equation for solving upper triangular matrix
*/

#include <math.h>

#define PRECISION	0.00001
#define MAX_ITER	200
#define PI		3.142

#define square(x) x*x;

matrix cholesky(const matrix &);
void find_disks(float *&, int &, const matrix &);
int find_max_in_row(const matrix &, int);
void gauss_elimination(matrix &);
float inv_power(const matrix &, matrix &);
matrix Jacobi(const matrix & A);
float power(const matrix &, matrix);
matrix rotation(const matrix &, int, int);
void shifted_power(const matrix &, matrix );
matrix solve_equation(matrix &);
float sum(const matrix &, int);
float sum(const matrix &, int, int);

matrix cholesky(const matrix & A)
{
    int i = 0, j = 0;
    matrix L(A.rows, A.cols);
    for (j = 0; j < A.cols; j++)
	for (i = 0; i < A.rows; i++) {
	    if (i == j) {			//calculate l[i][i]
		*(L.A + i * L.rows + j) =
		    sqrt(fabs(*(A.A + (i * A.rows) + i) - sum(L, i)));
	    } else if (i > j) {		//calculate l[i][j]
			if (*(L.A + j * L.rows + j) != 0) {
				*(L.A + i * L.rows + j) =
				(*(L.A + i * L.rows + j) -
				sum(L, i, j)) / (*(L.A + i * L.rows + j));
			} else {
				cout << "L[" << j << "][" << j << "] = 0";
				exit(0);
			}
	    } else {				//i<j
			*(L.A + i * L.rows + j) = 0;
		}
	}
    return L;
}

/* end of cholesky */

void find_disks(float *&disk, int &no_of_disks, const matrix & A)
{
    int i, j, k;
    float radius;
    for (i = 0; i < A.rows; i++) {
		radius = 0;
		for (j = 0; j < A.cols; j++) {
			if (i != j) {
				radius += *(A.A + (i * A.cols) + j);
			}
		}
		*(disk + (i * 2)) = *(A.A + (i * A.cols) + i) - radius;
		*(disk + (i * 2) + 1) = *(A.A + (i * A.cols) + i) + radius;
    }
    no_of_disks = A.rows;
	
    //check for duplicate disks
    for (i = 0; i <= no_of_disks; i++) {
		for (j = i + 1; j <= no_of_disks; j++) {
			//if found
			if (*(disk + (i * 2)) == *(disk + (j * 2))
			&& *(disk + (i * 2) + 1) == *(disk + (j * 2) + 1)) {
				//remove the row
				for (k = j; k < no_of_disks; k++) {
					*(disk + (k * 2)) = *(disk + (k + 1) * 2);
					*(disk + (k * 2) + 1) = *(disk + ((k + 1) * 2) + 1);
				}
				no_of_disks--;
				i = -1;
			}
		}
	}

    //print disks
    cout << endl << no_of_disks << " disk(s) found: " << endl;
    for (i = 0; i < no_of_disks; i++) {
		cout << *(disk + (i * 2)) << " " << *(disk + (i * 2) + 1) << endl;
	}
}

/* end of find_disks */

int find_max_in_row(const matrix & A, int i)
{
    float max;
    int pos = -1;
    max = 0;
    for (int j = i + 1; j < A.cols; j++) {
		float val = fabs(*(A.A + (i * A.rows) + j));
		if (val > max) {
			max = val;
			pos = j;
		}
    }
    if (max <= PRECISION) {
		pos = -1;
	}
    return pos;
}

/* end of find_max_in_row */

void gauss_elimination(matrix & A)
{
    float div, element;
    int i, j, k;

    for (i = 0; i < A.rows; i++) {
		//if diagonal element is zero
		if (*(A.A + (i * A.cols) + i) == 0) {
			for (k = i + 1; k < A.rows; k++)
			//if row with non-zero element found
			if (*(A.A + (k * A.cols) + i) != 0) {
				//row found ! swap rows
				for (j = i; j < A.cols; j++) {
				element = *(A.A + (i * A.cols) + j);
				*(A.A + (i * A.cols) + j) =
					*(A.A + (k * A.cols) + j);
				*(A.A + (k * A.cols) + j) = element;
				}
			}
			if (k == A.rows) {
				cout << "All diagonal elements can't be made non-zero."
					<< endl << "Method not applicable" << endl;
				exit(1);
			}
		} else {
			div = *(A.A + (i * A.cols) + i);
			for (j = i; j < A.cols; j++) {
				*(A.A + (i * A.cols) + j) =
					float (*(A.A + (i * A.cols) + j) / div);
			}
			for (k = i + 1; k < A.rows; k++) {
				element = *(A.A + (k * A.cols) + i);
				for (j = i; j < A.cols; j++)
					*(A.A + (k * A.cols) + j) =
					(*(A.A + (k * A.cols) + j)) -
					(*(A.A + (i * A.cols) + j)) * element;
			}
		}
    }
}

/* end of gaussian elimination */

float inv_power(const matrix & A, matrix & X)
{
    int i, j = 0, f = 0, m, n;
    float max1 = 1;		/* just to make initial condition true */
    float max2 = 0;
    matrix B;

    while ((fabs(max1 - max2)) > PRECISION) {
		B.augment(A, X);
		// cout << endl << "After augmentation" << endl << B;

		gauss_elimination(B);
		// cout << endl << "After GE" << endl << B;

		X = solve_equation(B);
		// cout << endl << "X" << endl << X;

		max1 = max2;		//saving old value of max
		max2 = *X.A;
		j = 0;
		//for finding maximum
		for (i = 1; i < X.rows; i++)
			if (fabs(*(X.A + i)) > fabs(max2))
			j = i;

		//dividing by max
		max2 = *(X.A + j);
		X.divide_by_max(max2);

		f++;
		if (f > MAX_ITER) {
			cout << endl << "Ill-conditioned matrix";
			return PRECISION;
		}
    }
    return 1 / max2;
}

matrix Jacobi(const matrix & A)
{
    int iterations = 0;
	int i = 0;
	matrix tmp = A;

    for (i = 0; i < tmp.rows; i++) {
		int j = -1;
		j = find_max_in_row(tmp, i);	//to find maximum of the row
		while (j != -1) {
			iterations++;
			if (iterations > MAX_ITER) {
				break;
			}
			tmp = rotation(tmp, i, j);
			j = find_max_in_row(tmp, i);
		}
		// cout << "Number of iterations " << iterations << endl;
		iterations = 0;
    }
    return tmp;
}

float power(const matrix & A, matrix X1)
{
    int i, j;

    float max = 0, new_max;
    matrix X2;
    while (1) {
		X2 = matrix(A.cols, 1);
		X2 = *(A * X1);
		new_max = X2.find_max();
		X2.divide_by_max(new_max);
		X1 = X2;
		if (fabs(max - new_max) < PRECISION) {
			break;
		} else {
			max = new_max;
		}
    }
    return max;
}

matrix rotation(const matrix & A, int i, int j)
{
    matrix R(A.rows, 'i');	// rotation matrix
	matrix RT;    			// transpose of rotation matrix

	float cos, sin, tan2, cos2;

    if (*(A.A + (i * A.rows) + i) == *(A.A + (j * A.rows) + j)) {
		sin = 1 / sqrt(2);
		cos = sin;
    } else {			// find the angle of rotation
	//aij                 aii                   -  ajj
		tan2 =
			2 * (*(A.A + (i * A.rows) + j)) /
			((*(A.A + (i * A.rows) + i)) - (*(A.A + (j * A.rows) + j)));
		cos2 = 1 / (sqrt(1 + tan2 * tan2));
		cos = sqrt((1 + cos2) / 2);
		sin = tan2 * cos2 / (2 * cos);
    }
    //form rotation matrix
    *(R.A + i * A.rows + i) = cos;
    *(R.A + i * A.rows + j) = (-1) * sin;
    *(R.A + j * A.rows + i) = sin;
    *(R.A + j * A.rows + j) = cos;
    //find transpose of the rotation matrix
   	RT.transpose(R);
    matrix b;
	b = *(RT * A);
    b = *(b * R);
    return b;
}

void shifted_power(const matrix & A, matrix X)
{
    float *disk, step = 0.05;
    int i, j, k;
    int no_of_disks;

    //find disks
    disk = new float[A.rows * 2];
    find_disks(disk, no_of_disks, A);

    //applying shifted power
    float s;
    float *lambda;
    int count = 0;
    float max;
    lambda = new float[A.rows];

	matrix B;
    for (i = 0; i < no_of_disks; i++) {
		while (count < A.rows) {
			if (*(disk + (i * 2)) == *(disk + (i * 2) + 1)) {
				// both values are same
				*(lambda + count) = *(disk + (i * 2));
				count++;
			} else {
				s = *(disk + (i * 2));
				while (s <= *(disk + (i * 2) + 1)) {
					B = A;
					//calc. B=A-sI
					for (j = 0; j < B.rows; j++) {
						*(B.A + (j * B.cols) + j) -= s;
					}
					max = inv_power(B, X);
					// (max<*(disk+(i*2)) || max>*(disk+(i*2)+1) ))
					if (fabs(max - PRECISION) >= PRECISION) {
						int k, flag = 0;
						for (k = 0; k < count; k++) {
							if (fabs(*(lambda + k) - (max + s)) < 1) {
								flag = 1;
							}
						}
						if (flag == 0) {
							break;
						}
					}
					s += step;
				}
				*(lambda + count) = max + s;
				count++;
			}
		}
    }

	cout << "Lambda values: " << endl;
	fflush(stdout);
    for (i = 0; i < count; i++) {
		cout << *(lambda + i) << endl;
	}
	cout << endl;
	delete disk;
	delete lambda;
}

matrix solve_equation(matrix & mat)
{
	matrix  tmp_mat = mat;

    int i;
    int rows = mat.rows;
    rows--;
    *(tmp_mat.A + rows) = *(mat.A + (rows * mat.cols) + mat.cols - 1);
    for (i = (rows - 1); i >= 0; i--) {
		*(tmp_mat.A + i) = *(mat.A + (i * mat.cols) + mat.cols - 1);
		for (int j = i + 1; j < mat.cols - 1; j++) {
			*(tmp_mat.A + i) =
			*(tmp_mat.A + i) -
			((*(mat.A + (i * mat.cols) + j)) * (*(tmp_mat.A + j)));
		}
    }

    return tmp_mat;
}

float sum(const matrix & L, int i)
{
    int j;
    float sum = 0.0;
    for (j = 0; j < i - 1; j++) {
		sum += *(L.A + (L.rows * i) + j) * (*(L.A + (L.rows * i) + j));
	}
    return sum;
}

float sum(const matrix & L, int i, int j)
{
    int k;
    float sum = 0.0;
    for (k = 0; k < j - 1; k++) {
		sum += *(L.A + L.rows * i + k) * (*(L.A + L.rows * j + k));
	}
    return sum;
}
