#pragma once
#include "Matrix.h"

void Matrix::mem_alloc()
{
	lambda_accur = new double[N] {};
	x = new double[N] {};

	A = new double*[N];
	H = new double*[N];
	for (int i = 0; i < N; ++i)
	{
		A[i] = new double[N] {};
		H[i] = new double[N] {};
	}
}

void Matrix::mem_clear()
{
	for (int i = 0; i < N; ++i)
	{
		delete[]A[i];
		delete[]H[i];
	}
	delete[]A;
	delete[]H;
	delete[]lambda_accur;
	delete[]x;
}

void Matrix::gen_vector(double* vec, int range)
{
	for (int i = 0; i < N; ++i)
		vec[i] = rand_number(range);
}


double Matrix::rand_number(int range)
{
	return -range + (rand() / (RAND_MAX / (2.0 * range)));
}


void Matrix::gen_omega(double* omega)
{
	double sum = 0;
	for (int i = 0; i < N; ++i)
	{
		omega[i] = rand_number();
		sum += omega[i] * omega[i];
	}
	sum = sqrt(sum);
	for (int i = 0; i < N; ++i)
		omega[i] = omega[i] / sum;
}

void Matrix::normalize(double* vector)
{
	double sum = 0;
	for (int i = 0; i < N; ++i)
		sum += vector[i]*vector[i];
	sum = sqrt(sum);
	for (int i = 0; i < N; ++i)
		vector[i] /= sum;
}

void Matrix::lambda_to_matrix(double** res)
{
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			res[i][j] = (i != j ? 0 : lambda_accur[i]);
}

void Matrix::multi_matrices(double** m1, double** m2, double** res)
{
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j) {
			res[i][j] = 0;
			for (int k = 0; k < N; ++k)
				res[i][j] += m1[i][k] * m2[k][j];
		}
}

std::string Matrix::ToString_matrix(double** matrix, std::string message)
{
	std::stringstream str(message);
	str << message;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
			str << std::setw(15) << matrix[i][j] << ' ';
		str << '\n';
	}
	return str.str();
}

std::string Matrix::ToString_vec(double* vec, std::string message)
{
	std::stringstream str;
	str <<  message << '\n';
	for (int i = 0; i < N; ++i)
		str << std::setw(15) << std::setprecision(25) << vec[i] << '\n';
	return str.str();
}

void Matrix::calculate()
{

	double* curr_nu = new double[N], * prev_nu = new double[N],
		* curr_x = new double[N], * next_x = new double[N],
		prev_sigma, curr_sigma;

	
	auto get_eps_g = [&prev_nu, &curr_nu, this]() -> double {
		double scalarMulti = 0, modulePrev = 0, moduleCurr = 0;
		for (int i = 0; i < N; i++)
		{
			modulePrev += /*H[i][iMax3] * H[i][iMax3];*/prev_nu[i] * prev_nu[i];
			moduleCurr += curr_nu[i] * curr_nu[i];
			scalarMulti += prev_nu[i] * curr_nu[i];// H[i][iMax3] * curr_nu[i];
		}
		double cos = scalarMulti / (sqrt(modulePrev) * sqrt(moduleCurr));

	/*	if (abs(cos) > 1) {
			cos = 1;
		}*/
		return std::acos(cos);
	};
	auto get_eps_l = [&prev_sigma, &curr_sigma, this]() -> double {
		return abs(/*lambda_accur[iMax3]*/prev_sigma - curr_sigma);
	};
	auto copy = [this](double* source, double* dest) {
		for (int i = 0; i < N; ++i)
			dest[i] = source[i];
	};
	auto calc_nu = [this, &copy, &curr_x, &curr_nu]() {
		copy(curr_x, curr_nu);
		normalize(curr_nu);
	};
	//auto multi_row_vec = [this](int i, double* dest, double* vec) {
	//	dest[i] = 0;
	//	for (int j = 0; j < N; ++j)
	//		dest[i] += A[i][j] * vec[j];
	//};
	auto calc_x = [this, &curr_nu, &next_x]() {
		for (int i = 0; i < N; ++i) {
			next_x[i] = 0;
			for (int j = 0; j < N; ++j)
				next_x[i] += A[i][j] * curr_nu[j];
		}
	};
	auto calc_sigma = [this, &curr_nu, &next_x]() {
		double res = 0;
		for (int i = 0; i < N; ++i)
			res += curr_nu[i] * next_x[i];
		return res;
	};

	//x0 -> 
	// nu0 -> x1 -> sigma0 
	// -> nu1 -> x2 -> sigma1
	//  -> nu2 -> x3 -> sigma2 -> ....
	gen_vector(curr_x, 100); //x0


	calc_nu();//nu0
	calc_x();//x1
	curr_sigma = calc_sigma();//sigma0
	do
	{
		copy(curr_nu, prev_nu);
		copy(next_x, curr_x);
		prev_sigma = curr_sigma;

	/*	std::cout << ToString_vec(prev_nu, "\nprev_nu:") << '\n';
		std::cout << "\nprevSigma: " << std::setprecision(25) << prev_sigma << '\n';
		std::cout << ToString_vec(curr_x, "\ncurr_x:") << '\n';*/

		calc_nu();//nu1
		calc_x();//x2
		curr_sigma = calc_sigma();//s1

		//std::cout << ToString_vec(next_x, "\nnext_x:") << '\n';
		//std::cout << ToString_vec(curr_nu, "\ncurr_nu:") << '\n';
		//std::cout << "currSigma: " << std::setprecision(25) << curr_sigma << '\n';
		++k;
	/*	if (k < 50) {
			std::cout <<  curr_sigma << '\n';
		}*/


//		out << "\ndiff nu: " << std::setprecision(25) << get_eps_g() << '\n';
//		out << "diff sigma: " << std::setprecision(25) << abs(curr_sigma-prev_sigma) 
//			<< '\n' << '\n' << '\n';
	} while (k <= K_max && (get_eps_g() >= eps_g || get_eps_l() >= eps_l));

	lambda = curr_sigma;
	copy(next_x, x);
	delete[] curr_nu;
	delete[] prev_nu;
	delete[] curr_x;
	delete[] next_x;
}

void Matrix::calc_H()
{
	double* omega = new double[N];
	gen_omega(omega);
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			H[i][j] = -2 * omega[i] * omega[j];
			if (i == j)
				H[i][j] += 1;
		}
	delete[]omega;
}

void Matrix::calc_A()
{
	double** L_lambda = new double* [N], **A_tmp = new double* [N];
	for (int i = 0; i < N; ++i) {
		L_lambda[i] = new double[N];
		A_tmp[i] = new double[N];
	}
	lambda_to_matrix(L_lambda);
	
	//std::cout << ToString_matrix(L_lambda, "got Lambda matrix:\n") << '\n';

	multi_matrices(H, L_lambda, A_tmp);
	//std::cout << ToString_matrix(A_tmp, "H*Lambda matrix:\n") << '\n';

	multi_matrices(A_tmp, H, A);
	//std::cout << ToString_matrix(A, "Prev*H_t:\n") << '\n';


	for (int i = 0; i < N; ++i)
	{
		delete[]L_lambda[i];
		delete[]A_tmp[i];
	}
	delete[]L_lambda;
	delete[]A_tmp;
}

void Matrix::findMax()
{
	iMax1 = -1, iMax2 = -1, iMax3 = -1;
	for (int i = 0; i < N; ++i) {
		if (abs(lambda_accur[iMax1]) <= abs(lambda_accur[i]) || iMax1 == -1) {
			iMax3 = iMax2;
			iMax2 = iMax1;
			iMax1 = i;
		}
		else if (abs(lambda_accur[iMax2]) <= abs(lambda_accur[i]) || iMax2 == -1) {
			iMax3 = iMax2;
			iMax2 = i;
		}
		else if (abs(lambda_accur[iMax3]) <= abs(lambda_accur[i]) || iMax3 == -1)
			iMax3 = i;
	}
	out << "\n iMAXes: " << iMax1 << ' ' << iMax2 << ' ' << iMax3 << '\n';
}

void Matrix::calculateNextA()
{
	auto multi = [this](double** matr, double value, uint iMax) {
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j) {
				matr[i][j] = value * H[i][iMax] * H[j][iMax];
			}
	};

	double **firstA = new double* [N], **secondA = new double* [N];
	for (int i = 0; i < N; ++i) {
		firstA[i] = new double[N] {};
		secondA[i] = new double[N] {};
	}

	multi(firstA, lambda_accur[iMax1], iMax1);
	multi(secondA, lambda_accur[iMax2], iMax2);

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			A[i][j] = A[i][j] - firstA[i][j] - secondA[i][j];


	out << "Calculated A^(2):\n" << ToString_A() << '\n';

	for (int i = 0; i < N; ++i) {
		delete[] firstA[i];
		delete[] secondA[i];
	}
	delete[] firstA;
	delete[] secondA;
}

void Matrix::calc_r()
{
	//n * n  n * 1 -> n * 1
	//A			x		vector column

	double* Ax = new double[N] {};
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j)
			Ax[i] += A[i][j] * x[j];
	}

	r = abs(Ax[0] - lambda * x[0]);
	for (int i = 1; i < N; i++)
	{
		double tmp = abs(Ax[i] - lambda * x[i]);
		if (r < tmp)
			r = tmp;
	}
	delete[] Ax;
}

Matrix::~Matrix()
{
	mem_clear();
}

Matrix::Matrix(uint N, uint K, int lambda_range, double epsilon): K_max(K), N(N), eps_g(epsilon), eps_l(epsilon)
{
	mem_alloc();
	gen_vector(lambda_accur, lambda_range);
	out << "Generated lambdas\n" << ToString_lambdas_accur() << '\n';
	calc_H();
	out << "Calculated H:\n" << ToString_H() << '\n';
	calc_A();
	out << "Calculated A:\n" << ToString_A() << '\n';

}

std::string Matrix::ToString_A()
{
	return ToString_matrix(A);
}
std::string Matrix::ToString_H()
{
	return ToString_matrix(H);
}
std::string Matrix::ToString_lambdas_accur()
{
	return ToString_vec(lambda_accur);
}

void Matrix::Solve()
{
	findMax();
	calculateNextA();
	calculate();
	calc_r();
	out << "K = " << k << '\n';
	out << "Lambda diff: " << get_l_accur() << '\n';
	out << "X diff: " << get_x_accur() << '\n';
	out << "r: " << get_r() << '\n';
	
}

double Matrix::get_r() const
{
	return r;
}

double Matrix::get_l_accur() const
{
	return abs(lambda - lambda_accur[iMax3]);
}

double Matrix::get_l() const
{
	return lambda;
}

uint Matrix::get_k() const
{
	return k;
}

double Matrix::get_x_accur() const
{
	double* Hdivided = new double[N], * xdivided = new double[N];
	double hdiv = 1 / H[0][iMax3], xdiv = 1 / x[0];
	Hdivided[0] = 1;
	xdivided[0] = 1;
	for (int i = 1; i < N; ++i) {
		Hdivided[i] = H[i][iMax3]*hdiv;
		xdivided[i] = x[i] * xdiv;
	}

	double max = abs(Hdivided[0] - xdivided[0]);
	for (int i = 1; i < N; ++i) {
		double tmp = abs(Hdivided[i] - xdivided[i]);
		if (max < tmp)
			max = tmp;
	}
	delete[] Hdivided;
	delete[] xdivided;
	return max;
}

std::string Matrix::ToString_X()
{
	return ToString_vec(x, "calculated x:\n");
}

