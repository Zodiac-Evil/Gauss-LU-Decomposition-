#include "decomposition.h"
#include <iostream>
#include <cmath>

using namespace std;

void swap_row(double* a, int i, int j, int n){
	for(int k = 0; k < n; ++k){
		double temp = 0.0;
		temp = *(a + i * n + k);
		*(a + i * n + k) = *(a + j * n + k);
		*(a + j * n + k) = temp;
	}
}

const int forward_replace(double const* lu, double* b, int n){

	for(int i = n - 1; i >= 0; --i){
		if(i == (n - 1))
			*(b + i) = *(b + i) / *(lu + i * n + i);
		else{
			double temp_total = 0.0;
			for(int j = i + 1; j < n; ++j){
				temp_total += *(b + j) * *(lu + i * n + j);
			}
			*(b + i) = (*(b + i) - temp_total) / *(lu + i * n + i);
		}
	}
	return 1;
}

const int backward_replace(double const* lu, double* b, int n){

	for(int i = 0; i < n; ++i){
		if(i == 0)
			*(b + i) = *(b + i);
		else{
			double temp_total = 0.0;
			for(int j = 0; j <= i - 1; ++j){
				temp_total += *(b + j) * *(lu + i * n + j);
			}
			*(b + i) = *(b + i) - temp_total;
		}
	}
	return 1;
}

bool lu(double* a, int* pivot, int n){
	double max = 0.0;

	for(int i = 0; i < n; ++i){
		*(pivot + i) = i;//默认是对角线上
	}

	for(int i = 0; i < n; ++i){//高斯列选主元法第i次迭代这个矩阵，到达第i列
		max = *(a + i * n + i);
		int m = 0, r = 0;
		for(m = i; m < n; ++m){//j为找到主元所在的行数
			if(max < *(a + m * n + i)){
				max = *(a + m * n + i);
				r = m;
			}
		}

		if(*(pivot + i) != i){
			swap_row(a, i, *(pivot + i), n);
			*(pivot + i) = r;//更新主元位置排列
			*(pivot + r) = i;
		}

		if(*(a + i * n + i) != 0){
			for(int j = i + 1; j < n; ++j){
				*(a + j * n + i) = *(a + j * n + i) / *(a + i * n + i);//就地填充L矩阵
			}
		}
		else{
			cout << i << endl;
			return true;
		}

		for(int j = i + 1; j < n; ++j){
			for(int k = i + 1; k < n; ++k){
				*(a + j * n + k) = *(a + j * n + k) - *(a + j * n + i) * *(a + i * n + k);//更新剩下的矩阵,也是最后的到U矩阵
			}
		}
	}
	return false;
}

bool guass(double const* lu, int const* p, double* b, int n){

	for(int i = 0; i < n; ++i){//对常数向量同样左乘置换矩阵
		if(*(p + i) != i){
			double temp = 0.0;
			temp = *(b + i);
			*(b + i) = *(b + *(p + i));
			*(b + *(p + i)) = temp;
		}
	}

	//先向后代替，然后向前代替
	const int i = backward_replace(lu, b, n);
	cout << "" << endl;
	const int j = forward_replace(lu, b, n);


	if(i == 1 && j == 1)
		return false;
	else
		return true;
}

void qr(double* a, double* d, int n){
	int i = 0, j = 0, k = 0, l = 0;
	double tempro = 0.0, norm = 0.0, *temp = new double[n];

	for(i = 0; i < n - 1; ++i){
		norm = 0.0;
		for(j = i; j < n; ++j){
			norm += *(a + n * j + i) * *(a + n * j + i);
		}

		if(*(a + n * i + i) > 0){
			norm = 0.0 - sqrt(norm);
		}
		else{
			norm = sqrt(norm);
		}

		tempro = 0.0;
		*(d + i) = norm;
		*(a + n * i + i) = *(a + n * i + i) - norm;

		for(j = i; j <= n - 1; ++j){
			tempro += *(a + n * j + i) * *(a + n * j + i);
		}
		tempro = sqrt(tempro);

		for(j = i; j <= n - 1; ++j){
			*(a + n * j + i) /= tempro;
		}

		for(k = i + 1; k < n; ++k){
			for(j = i; j < n; j++){
				tempro = 0.0;
				for(l = i; l < n; ++l){
					tempro += *(a + n * j + i) * *(a + n * l + i) * *(a + n * l + k);
				}
				*(temp + j) = *(a + n * j + k) - 2 * tempro;
			}
			for(j = i; j < n; ++j){
				*(a + n * j + k) = *(temp + j);
			}
		}
	}
	*(d + n - 1) = *(a + n * (n - 1) + n - 1);
	delete [] temp;
}

bool hshld(double const* qr, double const* d, double* b, int n){
	int i = 0, j = 0, l = 0;
	double buffer = 0.0, *temp = new double[n];

	for(i = 0; i < n - 1; ++i){
		for(j = i; j < n; ++j){
			buffer = 0.0;
			for(l = i; l < n; ++l){
				buffer += *(qr + n * l + i) * *(qr + n * j + i) * *(b + l);
			}
			*(temp + j) = *(b + j) - 2 * buffer;
		}

		for(j = i; j < n; ++j){
			*(b + j) = *(temp + j);
		}
	}

	for(j = n - 1; j > -1; --j){
		for(l = n - 1; l != j; --l){
			*(b + j) -= *(b + l) * *(qr + n * j + l);
		}
		if(*(d + j) == 0.0){
			return true;
		}
		*(b + j) /= *(d + j);
	}

	return false;
}
