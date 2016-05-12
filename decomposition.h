bool lu(double* a, int* pivot, int n);

bool guass(double const* lu, int const* p, double* b, int n);

void qr(double* a, double* d, int n);

bool hshld(double const* qr, double const* d, double* b, int n);

void swap_row(double* a, int i, int j, int n);//交换矩阵中i行和j行

const int forward_replace(double const* lu, double* b, int n);//U矩阵，向前替换法

const int backward_replace(double const* lu, double* b, int n);//L矩阵，向后替换法