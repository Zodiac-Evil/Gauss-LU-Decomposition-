/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.cc
 * Copyright (C) 2016 Java <java@java-Lenovo-Product>
 *
 * lu_decomposition is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * lu_decomposition is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include "decomposition.h"

using namespace std;

int main()
{
	ofstream out("solution.txt");

	int power = 0;
	out << "Hello world!" << endl;

	//第一题
	out << "----------第一题----------" << endl;
	enum {B1 = 10, B2 = 20, B3 = 30};

	double* matrix1 = new double[B1 * B1];
	double* matrix2 = new double[B2 * B2];
	double* matrix3 = new double[B3 * B3];

	double *b1 = new double[B1];
	double *b2 = new double[B2];
	double *b3 = new double[B3];

	for(int i = 0; i < B1; ++i){
		for(int j = 0; j < B1; ++j){
			*(matrix1 + i* B1 + j) = 1.0 / (i + 1 + j + 1 - 1);
		}

		double sum = 0.0;
		for(int j = i + 1; j <= (B1 + i); ++j){
			sum += 1.0 / j;
		}
		*(b1 + i) = sum;
	}

	for(int i = 0; i < B2; ++i){
		for(int j = 0; j < B2; ++j){
			*(matrix2 + i* B2 + j) = 1.0 / (i + 1 + j + 1 - 1);
		}

		double sum = 0.0;
		for(int j = i + 1; j <= (B2 + i); ++j){
			sum += 1.0 / j;
		}
		*(b2 + i) = sum;
	}

	for(int i = 0; i < B3; ++i){
		for(int j = 0; j < B3; ++j){
			*(matrix3 + i* B3 + j) = 1.0 / (i + 1 + j + 1 - 1);
		}

		double sum = 0.0;
		for(int j = i + 1; j <= (B3 + i); ++j){
			sum += 1.0 / j;
		}
		*(b3 + i) = sum;
	}

	double* matrix1_1 = new double[B1 * B1];
	double* matrix2_2 = new double[B2 * B2];
	double* matrix3_3 = new double[B3 * B3];

	double* b1_1 = new double[B1];
	double* b2_2 = new double[B2];
	double* b3_3 = new double[B3];

	for(int i = 0; i < B1; ++i){
		for(int j = 0; j < B1; ++j){
			*(matrix1_1 + i* B1 + j) = 1.0 / (i + 1 + j + 1 - 1);
		}

		double sum = 0.0;
		for(int j = i + 1; j <= (B1 + i); ++j){
			sum += 1.0 / j;
		}
		*(b1_1 + i) = sum;
	}

	for(int i = 0; i < B2; ++i){
		for(int j = 0; j < B2; ++j){
			*(matrix2_2 + i* B2 + j) = 1.0 / (i + 1 + j + 1 - 1);
		}

		double sum = 0.0;
		for(int j = i + 1; j <= (B2 + i); ++j){
			sum += 1.0 / j;
		}
		*(b2_2 + i) = sum;
	}

	for(int i = 0; i < B3; ++i){
		for(int j = 0; j < B3; ++j){
			*(matrix3_3 + i* B3 + j) = 1.0 / (i + 1 + j + 1 - 1);
		}

		double sum = 0.0;
		for(int j = i + 1; j <= (B3 + i); ++j){
			sum += 1.0 / j;
		}
		*(b3_3 + i) = sum;
	}

	int* pivot1 = new int[B1];

	lu(matrix1, pivot1, B1);
	guass(matrix1, pivot1, b1, B1);

	out << "LU分解的结果:" << B1 << endl;
	for(int i = 0; i < B1; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b1 + i) << endl;
	}
	double erfanshu = 0.0;
	for(int i = 0; i < B1; ++i){
		erfanshu += (*(b1 + i) - 1.0) * (*(b1 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "LU分解得到的数值解与理论解的二范数为：" << erfanshu << endl;
	delete [] matrix1;
	delete [] b1;
	delete [] pivot1;


	double* pivot1_1 = new double[B1];

	qr(matrix1_1, pivot1_1, B1);
	hshld(matrix1_1, pivot1_1, b1_1, B1);

	out << "QR分解的结果:" << B1 << endl;
	for(int i = 0; i < B1; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b1_1 + i) << endl;
	}
	erfanshu = 0.0;
	for(int i = 0; i < B1; ++i){
		erfanshu += (*(b1_1 + i) - 1.0) * (*(b1_1 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "QR分解得到的数值解与理论解的二范数为：" << erfanshu << endl;
	delete [] matrix1_1;
	delete [] b1_1;
	delete [] pivot1_1;

	int* pivot2 = new int[B2];
	lu(matrix2, pivot2, B2);
	guass(matrix2, pivot2, b2, B2);

	out << "LU分解的结果:" << B2 << endl;
	for(int i = 0; i < B2; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b2 + i) << endl;
	}

	erfanshu = 0;
	for(int i = 0; i < B2; ++i){
		erfanshu += (*(b2 + i) - 1.0) * (*(b2 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "LU分解得到的数值解与理论解的二范数为：" << erfanshu << endl;
	delete [] matrix2;
	delete [] b2;
	delete [] pivot2;


	double* pivot2_2 = new double[B2];

	qr(matrix2_2, pivot2_2, B2);
	hshld(matrix2_2, pivot2_2, b2_2, B2);

	out << "QR分解的结果:" << B2 << endl;
	for(int i = 0; i < B2; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b2_2 + i) << endl;
	}
	erfanshu = 0.0;
	for(int i = 0; i < B2; ++i){
		erfanshu += (*(b2_2 + i) - 1.0) * (*(b2_2 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "QR分解得到的数值解与理论解的二范数为：" << erfanshu << endl;
	delete [] matrix2_2;
	delete [] b2_2;
	delete [] pivot2_2;

	int* pivot3 = new int[B3];
	lu(matrix3, pivot3, B3);
	guass(matrix3, pivot3, b3, B3);

	out << "LU分解的结果:" << B3 << endl;
	for(int i = 0; i < B3; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b3 + i) << endl;
	}

	erfanshu = 0;
	for(int i = 0; i < B3; ++i){
		erfanshu += (*(b3 + i) - 1.0) * (*(b3 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "LU分解得到的数值解与理论解的二范数为：" << erfanshu << endl;
	delete [] matrix3;
	delete [] b3;
	delete [] pivot3;


	double* pivot3_3 = new double[B3];

	qr(matrix3_3, pivot3_3, B3);
	hshld(matrix3_3, pivot3_3, b3_3, B3);

	out << "QR分解的结果:" << B3 << endl;
	for(int i = 0; i < B3; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b3_3 + i) << endl;
	}
	erfanshu = 0.0;
	for(int i = 0; i < B3; ++i){
		erfanshu += (*(b3_3 + i) - 1.0) * (*(b3_3 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "QR分解得到的数值解与理论解的二范数为：" << erfanshu << endl;
	delete [] matrix3_3;
	delete [] b3_3;
	delete [] pivot3_3;

	//第二题
	out << "----------第二题----------" << endl;
	enum {A1 = 50, A2 = 100, A3 = 150};
	//悬空指针重利用
	matrix1 = new double[A1 * A1];
	matrix2 = new double[A2 * A2];
	matrix3 = new double[A3 * A3];

	b1 = new double[A1];
	b2 = new double[A2];
	b3 = new double[A3];

	pivot1 = new int[A1];
	pivot2 = new int[A2];
	pivot3 = new int[A3];

	matrix1_1 = new double[A1 * A1];
	matrix2_2 = new double[A2 * A2];
	matrix3_3 = new double[A3 * A3];

	b1_1 = new double[A1];
	b2_2 = new double[A2];
	b3_3 = new double[A3];

	pivot1_1 = new double[A1];
	pivot2_2 = new double[A2];
	pivot3_3 = new double[A3];

	for(int i = 0; i < A1; ++i){
		for(int j = 0; j < A1; ++j){
			if(i == j){
				*(matrix1 + i * A1 + j) = 1.0;
				*(matrix1_1 + i * A1 + j) = 1.0;
			}
			else if(i < j && j != A1 - 1){
				*(matrix1 + i * A1 + j) = 0.0;
				*(matrix1_1 + i * A1 + j) = 0.0;
			}
			else if(i > j){
				*(matrix1 + i * A1 + j) = 0.0 - 1.0;
				*(matrix1_1 + i * A1 + j) = 0.0 - 1.0;
			}
			else{
				*(matrix1 + i * A1 + j) = 1.0;
				*(matrix1_1 + i * A1 + j) = 1.0;
			}
		}

		if(i == 0){
			*(b1 + i) = 2;
			*(b1_1 + i) = 2;
		}
		else if( i > 0 && i < A1 - 1){
			*(b1 + i) = *(b1 + i - 1) - 1;
			*(b1_1 + i) = *(b1_1 + i - 1) - 1;
		}
		else{
			*(b1 + i) = 2 - A1;
			*(b1_1 + i) = 2 - A1;
		}
	}

	for(int i = 0; i < A2; ++i){
		for(int j = 0; j < A2; ++j){
			if(i == j){
				*(matrix2 + i * A2 + j) = 1.0;
				*(matrix2_2 + i * A2 + j) = 1.0;
			}
			else if(i < j && j != A2 - 1){
				*(matrix2 + i * A2 + j) = 0.0;
				*(matrix2_2 + i * A2 + j) = 0.0;
			}
			else if(i > j){
				*(matrix2 + i * A2 + j) = 0.0 - 1.0;
				*(matrix2_2 + i * A2 + j) = 0.0 - 1.0;
			}
			else{
				*(matrix2 + i * A2 + j) = 1.0;
				*(matrix2_2 + i * A2 + j) = 1.0;
			}
		}

		if(i == 0){
			*(b2 + i) = 2;
			*(b2_2 + i) = 2;
		}
		else if( i > 0 && i < A2 - 1){
			*(b2 + i) = *(b2 + i - 1) - 1;
			*(b2_2 + i) = *(b2_2 + i - 1) - 1;
		}
		else{
			*(b2 + i) = 2 - A2;
			*(b2_2 + i) = 2 - A2;
		}
	}

	for(int i = 0; i < A3; ++i){
		for(int j = 0; j < A3; ++j){
			if(i == j){
				*(matrix3 + i * A3 + j) = 1.0;
				*(matrix3_3 + i * A3 + j) = 1.0;
			}
			else if(i < j && j != A3 - 1){
				*(matrix3 + i * A3 + j) = 0.0;
				*(matrix3_3 + i * A3 + j) = 0.0;
			}
			else if(i > j){
				*(matrix3 + i * A3 + j) = 0.0 - 1.0;
				*(matrix3_3 + i * A3 + j) = 0.0 - 1.0;
			}
			else{
				*(matrix3 + i * A3 + j) = 1.0;
				*(matrix3_3 + i * A3 + j) = 1.0;
			}
		}

		if(i == 0){
			*(b3 + i) = 2;
			*(b3_3 + i) = 2;
		}
		else if( i > 0 && i < A3 - 1){
			*(b3_3 + i) = *(b3_3 + i - 1) - 1;
			*(b3_3 + i) = *(b3_3 + i - 1) - 1;
		}
		else{
			*(b3_3 + i) = 2 - A3;
			*(b3_3 + i) = 2 - A3;
		}
	}

	lu(matrix1, pivot1, A1);
	guass(matrix1, pivot1, b1, A1);

	out << "LU分解的结果:" << A1 << endl;
	for(int i = 0; i < A1; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b1 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A1; ++i){
		erfanshu += (*(b1 + i) - 1.0) * (*(b1 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "LU分解得到的数值解与理论解的二范数为：" << erfanshu << endl;

	qr(matrix1_1, pivot1_1, A1);
	hshld(matrix1_1, pivot1_1, b1_1, A1);

	out << "QR分解的结果:" << A1 << endl;
	for(int i = 0; i < A1; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b1_1 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A1; ++i){
		erfanshu += (*(b1_1 + i) - 1.0) * (*(b1_1 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "QR分解得到的数值解与理论解的二范数为：" << erfanshu << endl;

	lu(matrix2, pivot2, A2);
	guass(matrix2, pivot2, b2, A2);

	out << "LU分解的结果:" << A2 << endl;
	for(int i = 0; i < A2; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b2 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A2; ++i){
		erfanshu += (*(b2 + i) - 1.0) * (*(b2 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "LU分解得到的数值解与理论解的二范数为：" << erfanshu << endl;

	qr(matrix2_2, pivot2_2, A2);
	hshld(matrix2_2, pivot2_2, b2_2, A2);

	out << "QR分解的结果:" << A2 << endl;
	for(int i = 0; i < A2; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b2_2 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A2; ++i){
		erfanshu += (*(b2_2 + i) - 1.0) * (*(b2_2 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "QR分解得到的数值解与理论解的二范数为：" << erfanshu << endl;

	lu(matrix3, pivot3, A3);
	guass(matrix3, pivot3, b3, A3);

	out << "LU分解的结果:" << A3 << endl;
	for(int i = 0; i < A3; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b3 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A3; ++i){
		erfanshu += (*(b3 + i) - 1.0) * (*(b3 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "LU分解得到的数值解与理论解的二范数为：" << erfanshu << endl;

	qr(matrix3_3, pivot3_3, A3);
	hshld(matrix3_3, pivot3_3, b3_3, A3);

	out << "QR分解的结果:" << A3 << endl;
	for(int i = 0; i < A3; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b3_3 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A3; ++i){
		erfanshu += (*(b3_3 + i) - 1.0) * (*(b3_3 + i) - 1.0);
	}
	erfanshu = sqrt(erfanshu);
	out << "QR分解得到的数值解与理论解的二范数为：" << erfanshu << endl;

	//第二题的变形
	out << "----------第二题的变形----------" << endl;

	for(int i = 0; i < A1; ++i){
		for(int j = 0; j < A1; ++j){
			if(i == j){
				*(matrix1 + i * A1 + j) = 1.0;
				*(matrix1_1 + i * A1 + j) = 1.0;
			}
			else if(i < j && j != A1 - 1){
				*(matrix1 + i * A1 + j) = 0.0;
				*(matrix1_1 + i * A1 + j) = 0.0;
			}
			else if(i > j){
				*(matrix1 + i * A1 + j) = 0.0 - 1.0;
				*(matrix1_1 + i * A1 + j) = 0.0 - 1.0;
			}
			else{
				*(matrix1 + i * A1 + j) = 1.0;
				*(matrix1_1 + i * A1 + j) = 1.0;
			}
		}

		if(i < A1 - 1){
			*(b1 + i) = i + 2;
			*(b1_1 + i) = i + 2;
		}
		else{
			*(b1 + i) = A1;
			*(b1_1 + i) = A1;
		}
	}

	lu(matrix1, pivot1, A1);
	guass(matrix1, pivot1, b1, A1);

	out << "LU分解的结果:" << A1 << endl;
	for(int i = 0; i < A1; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b1 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A1; ++i){
		double temp = 0.0;

		for(int j = 0; j < i; ++j){
			temp += *(b1 + j);
		}

		if(i < A1 - 1){
			erfanshu += (*(b1 + A1 - 1) + *(b1 + i) - temp - (i + 2)) * (*(b1 + A1 - 1) + *(b1 + i) - temp - (i + 2));
		}
		else{
			erfanshu += (*(b1 + A1 - 1) - temp - (i + 1)) * (*(b1 + A1 - 1) - temp - (i + 1));
		}
	}
	erfanshu = sqrt(erfanshu);
	out << "LU分解得到的数值解的残量的二范数为：" << erfanshu << endl;

	qr(matrix1_1, pivot1_1, A1);
	hshld(matrix1_1, pivot1_1, b1_1, A1);

	out << "QR分解的结果:" << A1 << endl;
	for(int i = 0; i < A1; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b1_1 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A1; ++i){
		double temp = 0.0;

		for(int j = 0; j < i; ++j){
			temp += *(b1_1 + j);
		}

		if(i < A1 - 1){
			erfanshu += (*(b1_1 + A1 - 1) + *(b1_1 + i) - temp - (i + 2)) * (*(b1_1 + A1 - 1) + *(b1_1 + i) - temp - (i + 2));
		}
		else{
			erfanshu += (*(b1_1 + A1 - 1) - temp - (i + 1)) * (*(b1_1 + A1 - 1) - temp - (i + 1));
		}
	}
	erfanshu = sqrt(erfanshu);
	out << "QR分解得到的数值解的残量的二范数为：" << erfanshu << endl;

	for(int i = 0; i < A2; ++i){
		for(int j = 0; j < A2; ++j){
			if(i == j){
				*(matrix2 + i * A2 + j) = 1.0;
				*(matrix2_2 + i * A2 + j) = 1.0;
			}
			else if(i < j && j != A2 - 1){
				*(matrix2 + i * A2 + j) = 0.0;
				*(matrix2_2 + i * A2 + j) = 0.0;
			}
			else if(i > j){
				*(matrix2 + i * A2 + j) = 0.0 - 1.0;
				*(matrix2_2 + i * A2 + j) = 0.0 - 1.0;
			}
			else{
				*(matrix2 + i * A2 + j) = 1.0;
				*(matrix2_2 + i * A2 + j) = 1.0;
			}
		}

		if(i < A2 - 1){
			*(b2 + i) = i + 2;
			*(b2_2 + i) = i + 2;
		}
		else{
			*(b2 + i) = A2;
			*(b2_2 + i) = A2;
		}
	}

	lu(matrix2, pivot2, A2);
	guass(matrix2, pivot2, b2, A2);

	out << "LU分解的结果:" << A2 << endl;
	for(int i = 0; i < A2; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b2 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A2; ++i){
		double temp = 0.0;

		for(int j = 0; j < i; ++j){
			temp += *(b2 + j);
		}

		if(i < A2 - 1){
			erfanshu += (*(b2 + A2 - 1) + *(b2 + i) - temp - (i + 2)) * (*(b2 + A2 - 1) + *(b2 + i) - temp - (i + 2));
		}
		else{
			erfanshu += (*(b2 + A2 - 1) - temp - (i + 1)) * (*(b2 + A2 - 1) - temp - (i + 1));
		}
	}
	erfanshu = sqrt(erfanshu);
	out << "LU分解得到的数值解的残量的二范数为：" << erfanshu << endl;

	qr(matrix2_2, pivot2_2, A2);
	hshld(matrix2_2, pivot2_2, b2_2, A2);

	out << "QR分解的结果:" << A2 << endl;
	for(int i = 0; i < A2; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b2_2 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A2; ++i){
		double temp = 0.0;

		for(int j = 0; j < i; ++j){
			temp += *(b2_2 + j);
		}

		if(i < A2 - 1){
			erfanshu += (*(b2_2 + A2 - 1) + *(b2_2 + i) - temp - (i + 2)) * (*(b2_2 + A2 - 1) + *(b2_2 + i) - temp - (i + 2));
		}
		else{
			erfanshu += (*(b2_2 + A2 - 1) - temp - (i + 1)) * (*(b2_2 + A2 - 1) - temp - (i + 1));
		}
	}
	erfanshu = sqrt(erfanshu);
	out << "QR分解得到的数值解的残量的二范数为：" << erfanshu << endl;

	for(int i = 0; i < A3; ++i){
		for(int j = 0; j < A3; ++j){
			if(i == j){
				*(matrix3 + i * A3 + j) = 1.0;
				*(matrix3_3 + i * A3 + j) = 1.0;
			}
			else if(i < j && j != A3 - 1){
				*(matrix3 + i * A3 + j) = 0.0;
				*(matrix3_3 + i * A3 + j) = 0.0;
			}
			else if(i > j){
				*(matrix3 + i * A3 + j) = 0.0 - 1.0;
				*(matrix3_3 + i * A3 + j) = 0.0 - 1.0;
			}
			else{
				*(matrix3 + i * A3 + j) = 1.0;
				*(matrix3_3 + i * A3 + j) = 1.0;
			}
		}

		if(i < A3 - 1){
			*(b3 + i) = i + 2;
			*(b3_3 + i) = i + 2;
		}
		else{
			*(b3 + i) = A3;
			*(b3_3 + i) = A3;
		}
	}

	lu(matrix3, pivot3, A3);
	guass(matrix3, pivot3, b3, A3);

	out << "LU分解的结果:" << A3 << endl;
	for(int i = 0; i < A3; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b3 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A3; ++i){
		double temp = 0.0;

		for(int j = 0; j < i; ++j){
			temp += *(b3 + j);
		}

		if(i < A3 - 1){
			erfanshu += (*(b3 + A3 - 1) + *(b3 + i) - temp - (i + 2)) * (*(b3 + A3 - 1) + *(b3 + i) - temp - (i + 2));
		}
		else{
			erfanshu += (*(b3 + A3 - 1) - temp - (i + 1)) * (*(b3 + A3 - 1) - temp - (i + 1));
		}
	}
	erfanshu = sqrt(erfanshu);
	out << "LU分解得到的数值解的残量的二范数为：" << erfanshu << endl;

	qr(matrix3_3, pivot3_3, A3);
	hshld(matrix3_3, pivot3_3, b3_3, A3);

	out << "QR分解的结果:" << A3 << endl;
	for(int i = 0; i < A3; ++i){
		out.precision(20);
		out.setf(ios::fixed);
		out << "x" << (i + 1) << " = " << *(b3_3 + i) << endl;
	}

	erfanshu = 0.0;
	for(int i = 0; i < A3; ++i){
		double temp = 0.0;

		for(int j = 0; j < i; ++j){
			temp += *(b3_3 + j);
		}

		if(i < A3 - 1){
			erfanshu += (*(b3_3 + A3 - 1) + *(b3_3 + i) - temp - (i + 2)) * (*(b3_3 + A3 - 1) + *(b3_3 + i) - temp - (i + 2));
		}
		else{
			erfanshu += (*(b3_3 + A3 - 1) - temp - (i + 1)) * (*(b3_3 + A3 - 1) - temp - (i + 1));
		}
	}
	erfanshu = sqrt(erfanshu);
	out << "QR分解得到的数值解的残量的二范数为：" << erfanshu << endl;

	delete [] matrix1;
	delete [] matrix2;
	delete [] matrix3;

	delete [] b1;
	delete [] b2;
	delete [] b3;

	delete [] pivot1;
	delete [] pivot2;
	delete [] pivot3;

	delete [] matrix1_1;
	delete [] matrix2_2;
	delete [] matrix3_3;

	delete [] b1_1;
	delete [] b2_2;
	delete [] b3_3;

	delete [] pivot1_1;
	delete [] pivot2_2;
	delete [] pivot3_3;

	return 0;
}
