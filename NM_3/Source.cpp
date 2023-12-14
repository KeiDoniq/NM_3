
#include "Matrix.h"


int main()
{
	srand(GetTickCount64());
	int K = 1000;
	std::ofstream out2("out2.txt");
	for (int N = 10; N <= 50; N += 20)
	{
		for (int lambda_range = 2; lambda_range <= 50; lambda_range += 48)
			for (double eps = 1e-5; eps >= 1e-8; eps /= 1e3)
			{
				double av_acc_x(0), av_acc_l(0), av_r(0), av_k(0); int n_tests(50);
				for (int test = 0; test < n_tests; ++test)
				{
					Matrix M(N, K, lambda_range, eps);
					M.Solve();
					/*if (M.get_k() != 1001)
					{*/
						av_acc_x += M.get_x_accur();
						av_acc_l += M.get_l_accur();
						av_r += M.get_r();
						av_k += M.get_k();
					/*}
					else
						test -= 1;*/
				}
				av_acc_x /= n_tests;
				av_acc_l /= n_tests;
				av_r /= n_tests;
				av_k/= n_tests*1.0;
				out2 << "N = "<<N<<" lambda range = "<<lambda_range<<" eps = "<<eps<<"\nav_diff_lambda = " << av_acc_l << " av_diff_x = " << av_acc_x
				<< " av_r = " << av_r <<" av_k = " <<av_k << std::endl;
			}
	}
}