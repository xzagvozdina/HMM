#include "math.h"
#include <iostream>
#include <vector>
#include <fstream>
#include<ctime>
#include<cstdlib>
#include<time.h>
#include <random>

using namespace std;

double e = 2.71828182845904523536;

void array_output(double **a, int n1, int n2)
{
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << endl;
	}
}

void vector_output(double *a, int n)
{
	for (int i = 0; i < n; i++)
	{
		cout << a[i] << " ";
	}
	cout << endl;
}

void vector_output_int(int *a, int n)
{
	for (int i = 0; i < n; i++)
	{
		cout << a[i] << " ";
	}
	cout << endl;
}

void define_parameters(double **transition_matrix, double **emission_matrix, double *start_distr, int M, int K)
{
	fstream file_par;
	file_par.open("test.txt");

	for (int i = 0; i < M; i++)
		for (int j = 0; j < M; j++)
			file_par >> transition_matrix[i][j];

	for (int i = 0; i < M; i++)
		for (int j = 0; j < K; j++)
			file_par >> emission_matrix[i][j];

	for (int i = 0; i < M; i++)
		file_par >> start_distr[i];

	file_par.close();
}

int argmax(double **vector, int M, int i)
{
	int argmax_res = 0;
	double max_res = vector[i][0];

	for (int j = 1; j < M; j++)
		if (max_res < vector[i][j])
		{
			argmax_res = j;
			max_res = vector[i][j];
		}

	return argmax_res;
}

void generate_sequence(int *real_pi, int *x, double **transition_matrix, double **emission_matrix, double *start_distr, int M, int K, int L)
{
	/*GENERATING FIRST STATE*/
	double double_random = (double)(rand()) / (double)(RAND_MAX);
	real_pi[0] = 0;
	for (int i = 0; i < M; i++)
	{
		double_random -= start_distr[i];
		if (double_random < 0)
		{
			real_pi[0] = i;
			break;
		}
	}

	/*GENERATING FIRST SEQUENCE ELEMENT*/
	x[0] = argmax(emission_matrix, K, real_pi[0]);

	int flag_eqprob = 0;
	for (int j = 1; j < L; j++)
	{
		/*GENERATING NEW STATE*/
		double_random = (double)(rand()) / (double)(RAND_MAX);

		real_pi[j] = 0;
		for (int i = 0; i < M; i++)
		{
			double_random -= transition_matrix[real_pi[j - 1]][i];
			if (double_random < 0)
			{
				real_pi[j] = i;
				break;
			}
		}

		/*GENERATING NEW SEQUENCE ELEMENT*/
		for (int i = 1; i < K; i++)
			if (emission_matrix[real_pi[j]][0] = emission_matrix[real_pi[j]][i])
				flag_eqprob += 1;
		if (flag_eqprob = K-1)
		{
			double_random = (double)(rand()) / (double)(RAND_MAX);
			for (int i = 0; i < K; i++)
			{
				double_random -= emission_matrix[real_pi[j]][i];
				if (double_random < 0)
				{
					x[j] = i;
					break;
				}
			}
		}
		else
			x[j] = argmax(emission_matrix, K, real_pi[j]);

		flag_eqprob = 0;

		/*x[j] = argmax(emission_matrix, K, real_pi[j]);*/
	}
}

void viterbi_pi(int *vit_pi, int *real_pi, int *x, double **transition_matrix, double **emission_matrix, double *start_distr, int M, int K, int L)
{
	cout << "Viterbi: ";
	
	double **v = new double*[M];
	for (int i = 0; i < M; i++)
	{
		v[i] = new double[L + 1];
	}
	int **ptr = new int*[L];
	for (int i = 0; i < L; i++)
	{
		ptr[i] = new int[M];
	}
	/*INITIALISATION*/
	for (int i = 0; i < L + 1; i++)
		for (int j = 0; j < M; j++)
		{
			v[j][i] = 0.;
		}
	v[0][0] = 1.;
	for (int i = 0; i < L; i++)
		for (int j = 0; j < M; j++)
		{
			ptr[i][j] = 0.;
		}

	/*RECURSION*/
	double max_v = 0.;
	int a_max = 0;
	double *help_arr = new double[M];
	for (int i = 1; i < L + 1; i++)
	{
		for (int l = 0; l < M; l++)
		{
			max_v = v[0][i-1]*transition_matrix[0][l];
			a_max = 0;
			for (int k = 1; k < M; k++)
			{
				if (max_v < v[k][i - 1]*transition_matrix[k][l])
				{
					max_v = v[k][i - 1]*transition_matrix[k][l];
					a_max = k;
				}
			}
			v[l][i] = max_v*emission_matrix[l][x[i - 1]];
			/*v[l][i] = log(v[a_max][i - 1]) + log(transition_matrix[a_max][l]) + log(emission_matrix[l][x[i - 1]]);*/
			ptr[i - 1][l] = a_max;
		}
	}

	/*TERMINATION*/
	double vit_prob = v[0][L];
	int vit_pi_L = 0;
	for (int k = 1; k < M; k++)
		if (vit_prob < v[k][L])
		{
			vit_prob = v[k][L];
			vit_pi_L = k;
		}

	/*TRACEBACK*/
	vit_pi[L - 1] = ptr[L - 1][vit_pi_L];
	for (int i = L-1; i > 0; i--)
		vit_pi[i - 1] = ptr[i - 1][vit_pi[i]];
}

void for_back_ward_algorithm(double **prob, int *real_pi, int *x, double **transition_matrix, double **emission_matrix, double *start_distr, int M, int K, int L)
{
	cout << "Forward algorithm: ";

	double **f = new double*[M];
	for (int i = 0; i < M; i++)
	{
		f[i] = new double[L + 1];
	}

	/*INITIALISATION*/
	for (int i = 0; i < L + 1; i++)
	{
		for (int j = 0; j < M; j++)
		{
			f[j][i] = 0.;
		}
	}
	f[0][0] = 1.;

	/*RECURSION*/
	for (int i = 1; i < L + 1; i++)
	{
		for (int l = 0; l < M; l++)
		{
			for (int k = 0; k < M; k++)
				f[l][i] += f[k][i-1] * transition_matrix[k][l];
			f[l][i] *= emission_matrix[l][x[i-1]];
		}
	}

	/*TERMINATION*/
	double forward_prob = 0.;
	for (int k = 0; k < M; k++)
		forward_prob += f[k][L];
	cout << forward_prob << endl;

	cout << "Backward algorithm: ";

	double **b = new double*[M];
	for (int i = 0; i < M; i++)
	{
		b[i] = new double[L];
	}

	/*INITIALISATION*/
	for (int i = 0; i < L; i++)
		for (int j = 0; j < M; j++)
		{
			b[j][i] = 0.;
		}
	for (int k = 0; k < M; k++)
		b[k][L - 1] = 1.;

	/*RECURSION*/
	for (int i = L - 2; i >= 0; i--)
	{
		for (int k = 0; k < M; k++)
		{
			for (int l = 0; l < M; l++)
				b[k][i] += b[l][i + 1]*transition_matrix[k][l]*emission_matrix[l][x[i]];
		}
	}

	/*TERMINATION*/
	double backward_prob = 0.;
	for (int l = 0; l < M; l++)
		backward_prob += b[l][0] * transition_matrix[0][l] * emission_matrix[l][x[0]];
	cout << backward_prob << endl;

	cout << "Posterior decoding: ";

	double sum;
	for (int i = 0; i < L; i++)
	{
		sum = 0.;
		for (int k = 0; k < M; k++)
		{
			prob[i][k] = f[k][i+1] * b[k][i] / forward_prob;
			sum += prob[i][k];
		}
		/*if (sum > 1.)
			cout << sum << endl;*/
	}
}

int main()
{
	setlocale(LC_ALL, "Russian");

	srand(time(NULL));

	int M = 2;
	int K = 6;
	cout << "M = " << M << endl;
	/*cin >> M;
	cin.get();*/
	cout << "K = " << K << endl;
	/*cin >> K;
	cin.get();*/

	double **tr = new double*[M];
	for (int i = 0; i < M; i++)
		tr[i] = new double[M];
	double **em = new double*[M];
	for (int i = 0; i < M; i++)
		em[i] = new double[K];
	double *st = new double[M];

	define_parameters(tr, em, st, M, K);
	array_output(tr, M, M);
	array_output(em, M, K);
	vector_output(st, M);

	int L;
	cout << "L = ";
	cin >> L;
	cin.get();

	int *rpi = new int[L];
	int *x = new int[L];
	generate_sequence(rpi, x, tr, em, st, M, K, L);
	vector_output_int(rpi, L);
	vector_output_int(x, L);

	int *vpi = new int[L];
	viterbi_pi(vpi, rpi, x, tr, em, st, M, K, L);
	vector_output_int(vpi, L);

	double **p = new double*[L];
	for (int i = 0; i < L; i++)
		p[i] = new double[M];
	for_back_ward_algorithm(p, rpi, x, tr, em, st, M, K, L);
	
	ofstream file_out;
	file_out.open("out.dat");
	for (int i = 0; i < L; i++)
		file_out << rpi[i] << " ";
	for (int i = 0; i < L; i++)
		file_out << vpi[i] << " ";
	for (int i = 0; i < L; i++)
		file_out << p[i][1] << " ";
	file_out.close();

	cin.get();
}