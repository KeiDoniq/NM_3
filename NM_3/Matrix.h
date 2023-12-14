#pragma once
#include <Windows.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
using uint = unsigned int;



static std::ofstream out("out.txt");
class Matrix {
	uint N, K_max, k;
	int iMax1, iMax2, iMax3;
	double** A, ** H;
	double* lambda_accur;//собственные значения конструируемой матрицы A
	double lambda, eps_l, eps_g, r;
	double *x;


	void mem_alloc();
	void mem_clear();

	void gen_vector(double* vec, int range = 10);
	double rand_number(int range = 100);
	void calc_H();
	void calc_A();

	void gen_omega(double* omega);
	void normalize(double* vector);
	void lambda_to_matrix(double** res);
	void multi_matrices(double** m1, double** m2, double** res);

	std::string ToString_matrix(double** matrix, std::string message = "");
	std::string ToString_vec(double* vec, std::string message = "");

	void calculate();


	void findMax();
	void calculateNextA();
	void calc_r();
public:
	~Matrix();
	Matrix(uint N, uint K, int lambda_range, double epsilon);
	std::string ToString_A();
	std::string ToString_H();
	std::string ToString_lambdas_accur();
	std::string ToString_X();
	void Solve();
	double get_r() const;
	double get_x_accur();
	double get_l_accur() const;
	uint get_k() const;
	double get_l() const;
	

};
