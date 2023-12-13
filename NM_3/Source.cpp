
#include "Matrix.h"


int main()
{
	srand(GetTickCount64());
	int K = 1000;
	std::ofstream out2("out2.txt");
	for (int i = 10; i <= 50; i += 20)
	{
		for (int j = 2; j <= 50; j += 48)
			for (double t = 1e-5; t >= 1e-8; t /= 1e3)
			{
				Matrix M(i, K, j, t);
				M.Solve();

				out2 << "\nr = " << M.get_r()
					<< "\ndiff lambda = " << M.get_l_accur() << ",  lambda = " << std::setprecision(10) << M.get_l()
					<< "\ndiff x= " << M.get_x_accur() << "\n" << M.ToString_X() 
					<< "\nk = " << M.get_k() << "\n\n\n" << std::endl;


			}
	}
}