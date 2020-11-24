/*
Amit Walambe

functionality :
1.to perform foll. operations on matrix
	- read matrix from a file into array
	- display matrix
	- free memory of a matrix
2.manipulator used for display
*/

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <cstring>

using namespace std;

class matrix;
// class matrix *matrix_ptr;

ostream & setwidth(ostream &);
istream & operator>>(istream &, matrix &);
ostream & operator<<(ostream &, matrix &);

/* Definition of class matrix */
class matrix {
  private:
    float *A;
    int rows, cols;

     public:
	//methods provided by matrix class
     matrix();
     matrix(const matrix &);
     matrix(int, int);
     matrix(int, char);
    ~matrix();
    void mem_alloc();
    void mem_free();

    matrix *operator*(const matrix &) const;
    void operator =(const matrix &temp);
    void transpose(const matrix &);
    void augment(const matrix &, const matrix &);
    void divide_by_max(float max);
    float find_max();
    float *element(int row, int col) const;

    //friend functions of class matrix

    //      overloaded istream and ostream operators
    //      for accepting an displaying matrix
    //      definitions of foll. functions is in same file
    friend istream & operator>>(istream &, matrix &);
    friend ostream & operator<<(ostream &, matrix &);

    //      definitions of foll. functions is in nm.h
    friend matrix cholesky(const matrix &);
    friend void find_disks(float *&, int &, const matrix &);
    friend int find_max_in_row(const matrix &, int);
    friend void gauss_elimination(matrix &);
    friend float inv_power(const matrix &, matrix &);
    friend matrix Jacobi(const matrix &);
    friend float power(const matrix &, matrix);
    friend matrix rotation(const matrix &, int, int);
    friend void shifted_power(const matrix &, matrix);
    friend matrix solve_equation(matrix &);
    friend float sum(const matrix &, int);
    friend float sum(const matrix &, int, int);

    //modules currently under development
    friend void calc_sign_change(float val, int &sign_change,
				 const matrix &);
    friend void find_eigen(const matrix &, const float, const float,
			   float *&, int &);
    friend matrix givens_method(const matrix &);
    friend matrix tridiagonalise(const matrix &, int, int);
};

void matrix::mem_alloc()
{
    A = new float[rows * cols];
    memset(A, 0, (sizeof(*A) * rows * cols));
    // cout << "Allocating " << A << endl;
    // fflush(stdout);
}

void matrix::mem_free()
{
    // cout << "Deleting " << A << endl;
    // fflush(stdout);
    delete A;
}

// constructors of class matrix
matrix::matrix()		//default constructor
{
    A = NULL;
    rows = 0;
    cols = 0;
}

matrix::matrix(const matrix &mat)		//copy constructor
{
    rows = mat.rows;
    cols = mat.cols;
    mem_alloc();
    memcpy(A, mat.A, (sizeof(*A) * rows * cols));
}

matrix::matrix(int m, int n)	//when dimentions are specified
{
    rows = m;
    cols = n;
    mem_alloc();
}

matrix::matrix(int m, char ch)	// to generate identity matrix
{
    int i, j;
    if (toupper(ch) != 'I') {
        cout << "Invalid matrix type. Currently only identity "
            "matrix supported." << endl;
    }

    rows = m;
    cols = m;
    mem_alloc();
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            if (i == j) {
                *element(i, j) = 1;
            } else {
                *element(i, j) = 0;
            }
        }
    }
}
// end of definitions of constructors of class matrix

float *matrix::element(int row, int col) const
{
    return (A + (row * cols) + col);
}

matrix::~matrix()		// destructor
{
    // cout << "~matrix" << endl;
    if (A == NULL) {
	    return;
    } else {
	    mem_free();
    }
}

// definitions of methods of class matrix

matrix * matrix::operator *(const matrix & B) const
{
    if (cols != B.rows) {
        cout << endl << rows << "x" << cols <<
            " marix can't be multiplied by " << B.rows << "x" << B.
            cols << " matrix.";
        exit(1);
    }

    matrix *mat;
    mat = new matrix(rows, B.cols);

    int i, j, k;
    for (i = 0; i < rows; i++) {
        for (k = 0; k < B.cols; k++) {
            *(mat->element(i, k)) = 0;
            for (j = 0; j < cols; j++) {
                *(mat->element(i, k)) +=
                    (*element(i, j)) * (*B.element(j, k));
            }
        }
    }
    return mat;
}

void matrix::operator =(const matrix &temp)
{
    int i, j;
    if (A != NULL) {
    	mem_free();
    }
    rows = temp.rows;
    cols = temp.cols;
    mem_alloc();
    memcpy(A, temp.A, (sizeof(*A) * rows * cols));
    // for (i = 0; i < rows; i++) {
    //     for (j = 0; j < cols; j++) {
    //         *element(i, j) = *(temp.element(i, j));
    //     }
    // }
}

ostream & setwidth(ostream & o)
{
    o.width(12);
    o.precision(2);
    return (o);
}

/* overloaded istream and ostream operators
 * for accepting an displaying matrix */
istream & operator>>(istream & cin, matrix & mat)
{
    float val;
    char ch;
    string file_name;
    ifstream fp;

    if(mat.A != NULL) {
        delete mat.A;
    }

    do {
        cout << "\nEnter name file containing matrix : ";
        cin >> file_name;
        fp.open(file_name);
        if (fp.rdstate() != 0x00) {
            cout << endl << "There was error opening file. "
            "Do you want to try again (y/n) ? ";
            cin >> ch;
            if (toupper(ch) == 'N') {
                exit(1);
            }
        } else {
            break;
        }
    } while (1);

    /* File opened successfully.
     * Now, count number of columns. */
    while (!fp.eof()) {
        do {
            ch = fp.get();
        } while (ch != ' ' && ch != '\n');
        mat.cols++;
        if (ch == '\n') {
            mat.rows = 1;
            break;		/* done counting columns */
        }
    }
    /* count rows now */
    while (!fp.eof()) {
        do {
            ch = fp.get();
        } while (ch != '\n' && !fp.eof());
        if (fp.eof())
            break;		/* don't recount the last row */
        mat.rows++;
    }

    // allocating new size
    mat.mem_alloc();

    fp.close();
    fp.open(file_name);

    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            fp >> *(mat.A + (i * mat.cols) + j);
        }
    }
    return (cin);
}

ostream & operator<<(ostream & cout, matrix & A)
{
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            cout << setwidth << *A.element(i, j);
        }
        cout << endl;
    }
    return cout;
}


void matrix::augment(const matrix & a, const matrix & b)
{
    if (a.rows != b.rows) {
	    cout << endl <<
	    "Number of rows not equal. Augmentation not possible.";
	    exit(1);
    }

    if (A != NULL) {
	    mem_free();
    }
    rows = a.rows;
    cols = a.cols + 1;
    mem_alloc();

    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < a.cols; j++) {
            *element(i, j) = *a.element(i, j);
        }
        *element(i, j) = *b.element(i, 0);
    }
}

float matrix::find_max()
{
    int i, j;
    float max = 0;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if (*(A + (i * cols) + j) > max) {
                max = *(A + (i * cols) + j);
            }
        }
    }
    return max;
}

void matrix::divide_by_max(float max)
{
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            *(A + (i * cols) + j) = float (*(A + (i * cols) + j) / max);
        }
    }
}

void matrix::transpose(const matrix & a)
{
    if (A != NULL) {
	    mem_free();
    }
    rows = a.cols;
    cols = a.rows;
    mem_alloc();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            *(A + (j * rows) + i) = *(a.A + (i * a.rows) + j);
        }
    }
}

/* End of definition of class matrix */
