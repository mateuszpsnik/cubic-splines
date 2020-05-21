#include <iostream>
#include <vector>
#include <exception>
#include <fstream>

struct DivisionByZeroException : public std::exception
{
	char const* what() const throw()
	{
		return "The denominator was equal to zero.\nRemember that the value of B0 (0,0) cannot be equal to zero.\n";
	}
};

struct NonPositiveSizeException : public std::exception
{
	char const* what() const throw()
	{
		return "The size of the matrix must be positive.\n";
	}
};

using vec_doub = std::vector<double>;
using std::cout;
using std::cin;
using std::endl;

// this function initalizes the array used to print user's input correctly
void initialize_array(double** arr, vec_doub const& A, vec_doub const& B,
	vec_doub const& C, int size)
{
	for (int i = 0; i < size; i++)
	{
		arr[i][i] = B[i];
		if (i > 0)
			arr[i][i - 1] = A[i];
		if (i < size - 1)
			arr[i][i + 1] = C[i];
	}
}

// this function prints user's input (the equation)
void print_eq(double** arr, vec_doub const& F, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << arr[i][j] << "\t";
		}
		cout << "\t" << "u" << i << "\t\t" << F[i] << endl;
	}
}

// this is the function where the equation is solved
void tridiag_eq(vec_doub& A, vec_doub& B, vec_doub& C, vec_doub& F, vec_doub& u)
{
	double denom = 0.0;
	int n = u.size();
	vec_doub newc(n, 0.0);

	denom = B[0];
	u[0] = F[0] / denom;

	if (denom == 0.0)
		throw DivisionByZeroException{};

	// forward substitution
	for (int i = 1; i < n; i++)
	{
		newc[i] = C[i - 1] / denom;
		denom = B[i] - A[i] * newc[i];

		if (denom == 0.0)
			throw DivisionByZeroException{};

		u[i] = (F[i] - A[i] * u[i - 1]) / denom;
	}

	// back substitution
	for (int i = n - 2; i >= 0; i--)
	{
		u[i] -= newc[i + 1] * u[i + 1];
	}
}

double function(double coef1, double coef2, double coef3, double coef4, double t_i, double t_i1, double x)
{
	return coef1 * std::pow(t_i1 - x, 3) + coef2 * std::pow(x - t_i, 3)
		+ coef3 * (x - t_i) + coef4 * (t_i1 - x);
}

void calculate_print(vec_doub t, vec_doub y, std::string filename)
{
	vec_doub h;
	vec_doub u;
	vec_doub b;
	vec_doub v;
	vec_doub h_eq_C;

	for (size_t i = 0; i < t.size() - 1; i++)
	{
		h.push_back(t[i + 1] - t[i]);
		b.push_back((6 / h[i]) * (y[i + 1] - y[i]));
	}

	for (size_t i = 1; i < h.size(); i++)
	{
		u.push_back(2 * (h[i - 1] + h[i]));
		v.push_back(b[i] - b[i - 1]);
	}

	vec_doub z(t.size(), 0.0);
	vec_doub z_eq(v.size(), 0.0);

	tridiag_eq(h, u, h, v, z_eq);

	for (size_t i = 1; i < z.size() - 1; i++)
	{
		z[i] = z_eq[i - 1];
	}

	cout << "Splines:\n";
	
	std::ofstream file{ filename };
	for (size_t i = 0; i < t.size() - 1; i++)
	{
		double coef_1 = z[i] / (6 * h[i]);
		double coef_2 = z[i + 1] / (6 * h[i]);
		double coef_3 = ((y[i + 1] / h[i]) - (z[i + 1] * h[i] / 6));
		double coef_4 = ((y[i] / h[i]) - (z[i] * h[i] / 6));

		cout << "from " << t[i] << " to " << t[i + 1] << ": " << coef_1 << "(" 
			<< t[i + 1] << " - x)^3 + " << coef_2 << "(x - " << t[i] << ")^3 + " 
			<< coef_3 << "(x - " << t[i] << ") + " << coef_4 << "(" << t[i + 1] << " - x)\n";

		for (double j = t[i]; j < t[i + 1]; j += 0.01)
		{
			file << j << "," << function(coef_1, coef_2, coef_3, coef_4, t[i], t[i + 1], j) << endl;
		}
	}
}

int main()
{
	cout << "Sine:\n";

	vec_doub t_sine{ -6, -5, -3, -2, 0, 1.5, 3, 4, 6, 7, 8 };

	vec_doub y_sine;

	for (size_t i = 0; i < t_sine.size(); i++)
	{
		y_sine.push_back(std::sin(t_sine[i]));
	}

	calculate_print(t_sine, y_sine, "data_sine.txt");

	cout << "My example (y = x^3):\n";

	vec_doub t_my{ -12, -10, -8, -7, -6, -5, -3, -2, 0, 1.5, 3, 4, 6, 7, 8, 10, 12 };
	vec_doub y_my{ -1728, -1000, -512, -343, -216, -125, -27, -8, 0, 3.375, 27, 64, 216, 343, 512, 1000, 1728 };

	calculate_print(t_my, y_my, "data_my.txt");
}