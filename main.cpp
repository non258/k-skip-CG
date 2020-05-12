#include <cmath>
#include <iostream>
#include <vector>

void mvm(std::vector<double> A, std::vector<double> x, std::vector<double> b, int size)
{
  for (int i = 0; i < size; i++)
  {
    b[i] = 0.;
    for (int j = 0; j < size; j++)
    {
      b[i] += A[i * size+j] * x[j];
    }
  }
}

double ip(std::vector<double> a, std::vector<double> b, int size)
{
  double value = 0.;

  for (int i = 0; i < size; i++)
  {
    value += a[i] * b[i];
    return value;
  }
}

double norm(std::vector<double> vec, int size)
{
  double s = 0.0;

  for (int i = 0; i < size; i++)
  {
    s += vec[i] * vec[i];
  }

  return sqrt(s);
}

int main()
{
  int size = 3;
  int k = 1;
  int q = 100000;
  double eps = pow(2, -50);
  std::vector<double> A(size * size, 0);
  std::vector<double> b(size, 0);

  A[0*size+0]=1.;A[0*size+1]=2.;A[0*size+2]=-3.;
  A[1*size+0]=2.;A[1*size+1]=5.;A[1*size+2]=-4.;
  A[2*size+0]=-3.;A[2*size+1]=-4.;A[2*size+2]=8.;

  b[0]=-4;b[1]=-3;b[2]=12;

  double alpha, beta;
  double r0;

  std::vector<double> residual(q, 0);

  std::vector<double> x(size, 0);
  std::vector<double> Ax(size, 0);
  std::vector<std::vector<double>> r(k+2, std::vector<double>(size, 0));
  std::vector<std::vector<double>> p(k+3, std::vector<double>(size, 0));

  std::vector<double> a(2 * k+2, 0);
  std::vector<double> f(2 * k+4, 0);
  std::vector<double> c(2 * k+2, 0);

  mvm(A, x, Ax, size);
  for (int i = 0; i < size; i++)
  {
    r[0][i] = b[i] - Ax[i];
  }

  for (int i = 0; i < size; i++)
  {
    p[0][i] = r[0][i];
  }

//  r0 = sqrt(ip(r, r, size));

  for (int i = 0; i < q; i++)
  {
    std::cout << "i: " << i << std::endl;
    residual[i] = norm(r[0], size) / norm(b, size);
    if (residual[i] < eps)
    {
      std::cout << "break" << std::endl;
      break;
    }

    for (int j = 1; j < k+1; j++)
      for (int k = 0; k < size; k++)
        for (int l = 0; l < size; l++)
          r[j][k] += A[k*size+l]*r[j-1][l];

    for (int j = 1; j < k+1; j++)
      for (int k = 0; k < size; k++)
        for (int l = 0; l < size; l++)
          p[j][k] += A[k*size+l]*p[j-1][l];


     for (int j = 0; j < 2*k+1; j+=2)
     {
       int jj = j / 2;
       for (int k = 0; k < size; k++)
         a[j] = ip(r[jj], r[jj], size);
       for (int k = 0; k < size; k++)
         a[j+1] = ip(r[jj], r[jj+1], size);
     }


    for (int j = 0; j < 2*k+3; j+=2)
    {
      int jj = j / 2;
      for (int k = 0; k < size; k++)
        f[j] = ip(p[jj], p[jj], size);
      for (int k = 0; k < size; k++)
        f[j+1] = ip(p[jj], p[jj+1], size);
    }

    for (int j = 0; j < 2*k+1; j+=2)
    {
      int jj = j / 2;
      for (int k = 0; k < size; k++)
        c[j] = ip(r[jj], p[jj], size);
      for (int k = 0; k < size; k++)
        f[j+1] = ip(r[jj], p[jj+1], size);
    }

    alpha = a[0] / f[1];
    beta = pow(alpha, 2) * f[2] / a[0] -1;
    for (int j = 0; j < p.size(); j++)
      x[j] += alpha * p[0][j];

    for (int j = 0; j < r.size(); j++)
      r[0][j] -= alpha * p[0][j];
    for (int j = 0; j < r.size(); j++)
      p[0][j] = r[0][j] + beta * p[0][j];
    for (int j = 0; j < size; j++)
      for (int k = 0; k < size; k++)
        p[1][j] += A[j*size+k] * p[0][k];

  }

  for (int i = 0; i < size; i++)
    std::cout << x[i] << std::endl;

  return 0;
}