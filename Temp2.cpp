#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

double func(double x, double t)
{
    return (x * x - 2 * t);
}
double leftfunc(double t)
{
    return 0.0;
}

double rightfunc(double t)
{
    return 1200.0 * 1200.0 * t;
}

double zerotimefunc(double x)
{
    return 0.0;
}
void Print(double** data, int string,int collumn, double dt, double dx)
{
  for (int i = 0; i < string; i++)
      for (int j = 0; j < collumn; j++)
      {
         std::cout << std::endl << " (" << i * dx << ", " << dt * j << ")"; 
         std::cout << std::endl << data[i];
      }
}

void FPrint(double** data, int string, int collumn, std::ofstream& out, double dt, double dx)
{
  out.open("..//Temp//Temperature.txt");
  out << "(x; y)";
  for (int j = 0; j <= string; j++)
  for (int i = 0; i <= collumn; i++)
  {
    out << std::endl << " (" << i * dx << "; " << dt * j << ")";
    out << std::endl << data[j][i];
  }
  out.close();
}

double** Temperature(double materialkoef, int xsteps, int tsteps, double maxT, double maxX)
{
    double dx = maxX / xsteps;
    double dt = maxT / tsteps;
    double** TemperatureData = new double* [tsteps+1]; // выделение памяти для указателей 
    for (int i = 0; i <= tsteps; i++) // выделение памяти для элементов
        TemperatureData[i] = new double [xsteps + 1]; 
    for (int i = 0; i <= tsteps; i++) // задание крайних значений по x
        TemperatureData[i][0] = leftfunc(dt * i); 
    for (int i = 0; i <= tsteps; i++) // задание крайних правых значений по x
        TemperatureData[i][xsteps] = rightfunc(dt * i); 
    for (int i = 0; i <= xsteps; i++) // значение в нулевой момент времени.
        TemperatureData[0][i] = zerotimefunc(dx*i);
    int num = 2;
    omp_set_num_threads(num);
    double t1=omp_get_wtime();
#pragma omp parallel for
    for (int i = 1; i < tsteps; i++)
        for(int j = 1; j <= xsteps; j++)
        {
           TemperatureData[i][j] = materialkoef * dt * ((TemperatureData[i - 1][j + 1]) - 2 * (TemperatureData[i - 1] [j]) + (TemperatureData[i - 1][j - 1])) / dx / dx + TemperatureData[i - 1][j] + func((j - 1) * dx, i * dt) * dt;
        }
    double t2=omp_get_wtime();
    std::cout<<"num threads: "<< num <<std::endl;
    std::cout<<"time: "<<t2-t1<<"s"<<std::endl;
    return TemperatureData;
}
double** Error(double** res, double** Temperature, int size_n, int size_m, double maxT, double maxX)
{
    double dx = maxX / size_n;
    double dt = maxT / size_m;
    res = new double* [size_m + 1];
    for (int i = 0; i <= size_m; i++)
        res[i] = new double[size_n];
    for (int i = 0; i <= size_m; i++)
        for(int j = 0; j <= size_n; j++)
        {
            res[i][j] = fabs(Temperature[i][j] - func(dx * j, dt * i));
        }
    return res;
}
int main()
{
    double** TemperatureData = nullptr;
    int n = 20;
    int m = 400;
    double materialkoef = 1;
    double maxT = 5;
    double maxX = 3;
    double dx = maxX / n;
    double dt = maxT / m;
    TemperatureData = Temperature(materialkoef,  n, m,  maxT, maxX);
    std::ofstream out;
    //FPrint(TemperatureData, n, m, out, dt, dx);
    double** ErrorData = nullptr;
    ErrorData = Error(ErrorData, TemperatureData, n, m, maxT, maxX);
    FPrint(ErrorData, n, m, out, dt, dx);
    delete[] TemperatureData;
    delete[] ErrorData;
    return 0;
}
 
