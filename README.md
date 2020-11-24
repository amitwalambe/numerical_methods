# numerical_methods
This code was written circa 2001-2; during my masters days. Instead of boiler-plating matrix related operations for numerical methods (and partially taking fancy to being able to multiply matrices using '*' operator etc.), I implemented the code in C++.
Looking back, so many things are cringeworthy (generous use of friend functions, making copies of matrices which will be inefficient when operating on large matrices). But a good trip down the memory lane nonetheless.

$  g++ matrix.cpp -o nm
$ ./nm
Use any of the files in test_matrices folder as input.


