# First-repository
My first repository
```
//#include <iostream>
//#include <cmath>
//int main() {
//    double a, b, c;
//    std::cin >> a >> b >> c;
//    double eps = 1e-10;
//    if (a == 0) {
//        if (b == 0) {
//            if (c == 0) {
//                std::cout << "X any value";
//            }
//            else {
//                std::cout << "No solutions";
//            }
//        }
//        else {
//            double x = -c / b;
//        if (abs(b * x + c) < eps * abs(c)) {
//                std::cout << x;
//        }
//          else {
//                std::cout << "No solutions";
//          }
//        }
//        return 0;
//    }
//    double D = b * b - 4 * a * c;
//    if (D > eps) {
//        double x1 = (-b - sqrt(D)) / (2 * a);
//        double x2 = (c / (a * x1));
//        std::cout << x1 << " " << x2;
//    }
//    else if (abs(D) < eps) {
//        double x = -b / (2 * a);
//    if (abs(a * x * x + b * x + c) < eps * (abs(a * x * x) + abs(b * x) + abs(c))) {
//            std::cout << x;
//        }
//        else {
//            std::cout << "No solutions";
//        }
//    }
//    else {
//        std::cout << "No solutions";
//    }
//     return 0;
//}


















#include <iostream>
#include <cmath>


int main() {
    double a, b, c;
    std::cin >> a >> b >> c;

   
    double eps = 1e-11;
    if (a == 0) {
        if (b == 0) {
            if (c == 0) {
                std::cout << "X any value";
            }
            else {
                std::cout << "No solutions";
            }
        }
        else {
            double x = -c / b;
            
            if (abs(b * x + c) <= eps * std::max(abs(b * x), abs(c))) {
                std::cout << x;
            }
            else {
                std::cout << "No solutions";
            }
        }
        return 0;
    }

    double D = b * b - 4 * a * c;

    if (D > eps * std::max(b * b, 4 * abs(a * c))) {
        double x1, x2;
        x1 = (-b - sqrt(D)) / (2 * a);
        x2 = (-b + sqrt(D)) / (2 * a);
        if (abs(a * x1 * x1 + b * x1 + c) <= eps * (abs(a * x1 * x1) + abs(b * x1) + abs(c))) {
            std::cout << x1 << " " << x2;
        }
        else {
            std::cout << "No solutions";
        }
    }
    else if (abs(D) <= eps * std::max(b * b, 4 * abs(a * c))) {
        double x = -b / (2 * a);
        if (abs(a * x * x + b * x + c) <= eps * (abs(a * x * x) + abs(b * x) + abs(c))) {
            std::cout << x;
        }
        else {
            std::cout << "No solutions";
        }
    }
    else {
        
        std::cout << "No solutions";
    }

    return 0;
}
```
```
#include <iostream>
#include <algorithm>
int main()
{
	int n, m;
	std::cout << "Enter lines: " << std::endl;
	std::cin >> m;
	std::cout << "Enter columns: " << std::endl;
	std::cin >> n;
	int** mas;
	mas = new int* [m];
	for (int i = 0; i < m; ++i) {
		mas[i] = new int[n];
	}
	int max;
	std::cout << "If you want your own massive, please, enter : '1' " << std::endl;
	std::cout << "If you want random massive, please, enter : '2' " << std::endl;
	int choose;
	std::cin >> choose;
	srand(time(NULL));
	if (choose == 1) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				std::cin >> mas[i][j];
			}
		}
	}
	if (choose == 2) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				mas[i][j] = rand();
			}
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	int min = 0;

	for (int j = 0; j < n; ++j) {
		min += mas[0][j];
	}
	/*std::cout << min;*/
	
	//на строчки
	int cnt = 0;
	int cntstring = 0;
	for (int i = 0; i < m; ++i) {
		cnt = 0;
		for (int j = 0; j < n; ++j) {
			cnt += mas[i][j];

		}
		if (cnt < min) {
			min = cnt;
			cntstring = i;
		}
	}

	/*++cntstring;*/
	//std::cout << cntstring << std::endl;
	for (int i = cntstring; i < (m - 1); ++i) {
		int* temp = mas[i];
		mas[i] = mas[i + 1];
		mas[i + 1] = temp;
	}






	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;




	max = mas[0][0];
	// столбцы
	int cnt2 = 0;
	int cntcolumns = 0;
	for (int j = 0; j < n; ++j) {
		cnt2 = 0;
		for (int i = 0; i < m; ++i) {
			cnt2 += mas[i][j];
		}
		if (cnt2 > max) {
			max = cnt2;
			cntcolumns = j;
		}
	}
	/*std::cout << cntcolumns<<std::endl;*/
	for (int j = cntcolumns; j < (n - 1); ++j) {
		for (int i = 0; i < n; ++i) {
			int temp = mas[i][j];
			mas[i][j] = mas[i][j + 1];
			mas[i][j + 1] = temp;
		}

	}





	/*for (int j = cntcolumns; j < n - 1; ++j) {
		for (int i = 0; i < m; ++i) {
			int* ptr_col1 = &mas[i][j];
			int* ptr_col2 = &mas[i][j + 1];
			std::swap(*ptr_col1, *ptr_col2);
		}
	}*/




	/*for (int j = cntcolumns; j < n - 1; ++j) {

			int* ptr_col1 = &mas[i][j];
			int* ptr_col2 = &mas[i][j + 1];
			std::swap(*ptr_col1, *ptr_col2);

	}*/





	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i][j] << " ";
		}
		std::cout << std::endl;
	}
	for (int i = 0; i < m; ++i) {
		delete[] mas[i];
	}
	delete[] mas;
	return 0;
}
```
```
#include <iostream>
#include <cmath>

int main() {
    double a, b, c, x, x1, x2, D;
    std::cin >> a >> b >> c;

    if (a == 0) {
        if (b == 0) {
            if (c == 0) {
                std::cout << "X any value";    
            }
            else {
                std::cout << "No solutions";        
            }
        }
        else {
            x = -c / b;
            std::cout << x;
        }
        return 0;
    }

   
    D = b * b - 4 * a * c;

    if (D > 0) {
        x1 = (-b - sqrt(D)) / (2 * a);
        x2 = (c / (a * x1));
        std::cout   << x1 << " " << x2;
    }
    else if (D == 0) {
        x = -b / (2 * a);
        std::cout  << x;
    }
    else {
        std::cout << "No solutions";
    }

    return 0;
}
```
```
#include <iostream>
#include <algorithm>
int main()
{
	int n, m;
	std::cout << "Enter lines: " << std::endl;
	std::cin >> m;
	std::cout << "Enter columns: " << std::endl;
	std::cin >> n;
	int** mas;
	mas = new int* [m];
	for (int i = 0; i < m; ++i) {
		mas[i] = new int[n];
	}
	int max;
	std::cout << "If you want your own massive, please, enter : '1' " << std::endl;
	std::cout << "If you want random massive, please, enter : '2' " << std::endl;
	int choose;
	std::cin >> choose;
	if (choose == 1) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				std::cin >> mas[i][j];
			}
		}
	}
	if (choose == 2) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				mas[i][j] = rand();
			}
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	int min=0;
	
	for (int j = 0; j < n; ++j) {
		min += mas[0][j];
	}
	/*std::cout << min;*/
	max = mas[0][0];
	//на строчки
	int cnt = 0;
	int cntstring = 0;
	for (int i = 0; i < m; ++i) {
		cnt = 0;
		for (int j = 0; j < n; ++j) {
			cnt += mas[i][j];

		}
		if (cnt < min) {
			min = cnt;
			cntstring = i;
		}
	}
	
	/*++cntstring;*/
	//std::cout << cntstring << std::endl;
	/*for (int i = cntstring; i < (m - 1); ++i) {
		for (int j = 0; j < n; ++j) {
			std::swap(mas[i][j], mas[i + 1][j]);
		}
	}*/

	int* help, *help2;
	for (int i = cntstring; i < (m - 1); ++i) {
	    mas[i] = help;
		mas[i + 1] = help2;
		help2 = mas[i];
		mas[i+1] = help;
		

	}




	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;





	// столбцы
	int cnt2 = 0;
	int cntcolumns = 0;
	for (int j = 0; j < n; ++j) {
		cnt2 = 0;
		for (int i = 0; i < m; ++i) {
			cnt2 += mas[i][j];
		}
		if (cnt2 > max) {
			max = cnt2;
			cntcolumns = j;
		}
	}
	/*std::cout << cntcolumns<<std::endl;*/
	for (int j = cntcolumns; j < (n - 1); ++j) {
		for (int i = 0; i < n; ++i) {
			std::swap(mas[i][j], mas[i][j + 1]);
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i][j] << " ";
		}
		std::cout << std::endl;
	}
	for (int i = 0; i < m; ++i) {
		delete[] mas[i];
	}
	delete[] mas;
	return 0;
}
```
```
#include <iostream>
#include <algorithm>

int main() {
	int m, n;
	std::cout << "Enter lines: " << std::endl;
	std::cin >> m;
	std::cout << "Enter columns: " << std::endl;
	std::cin >> n;

	
	int* mas = new int[ m*n];

	int choose;
	std::cout << "If you want your own massive, please, enter: '1'" << std::endl;
	std::cout << "If you want random massive, please, enter: '2'" << std::endl;
	std::cin >> choose;
	srand(time(NULL));
	if (choose == 1) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				std::cin >> mas[i * n + j];
			}
		}
	}
	else if (choose == 2) {
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				mas[i * n + j] = rand();
			}
		}
	}

	
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i * n + j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	
	int min = 0;
	for (int j = 0; j < n; ++j) {
		min += mas[0 * n + j];
	}
	int max = mas[0];
	
	int cntstring = 0;
	int cnt;
	for (int i = 0; i < m; ++i) {
		cnt = 0;
		for (int j = 0; j < n; ++j) {
			cnt += mas[i * n + j];
		}
		if (cnt < min) {
			min = cnt;
			cntstring = i;
		}
	}

	
	for (int i = cntstring; i < m - 1; ++i) {
		for (int j = 0; j < n; ++j) {
			std::swap(mas[i * n + j], mas[(i + 1) * n + j]);
		}
	}
	std::cout << std::endl;
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i * n + j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	
	
	int cntcolumns = 0;
	int cnt2;
	for (int j = 0; j < n; ++j) {
		cnt2 = 0;
		for (int i = 0; i < m; ++i) {
			cnt2 += mas[i * n + j];
		}
		if (cnt2 > max) {
			max = cnt2;
			cntcolumns = j;
		}
	}

	
	for (int j = cntcolumns; j < n - 1; ++j) {
		for (int i = 0; i < m; ++i) {
			std::swap(mas[i * n + j], mas[i * n + (j + 1)]);
		}
	}

	
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i * n + j] << " ";
		}
		std::cout << std::endl;
	}

	delete[] mas;
	return 0;
}

```
```
#include <iostream>
#include <cmath>
#include <chrono>

void AllocateMemory(double**& mas, int size) {
	mas = new double* [size];
	for (int i = 0; i < size; ++i) {
		mas[i] = new double[size];
	}

}
void FillAB(double**& mas1, double**& mas2, int size) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			mas1[i][j] = std::sin(i - j);
		}
	}
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			mas2[i][j] = 1.0 / (i + j + 1.0);
		}
	}
}
void FillZero(double**& mas, int size) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			mas[i][j] = 0;
		}
	}
}
double MultiplyIJK(double**& mas1, double**& mas2, double**& mas3, int size) {
	auto start = std::chrono::steady_clock::now();
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			for (int k = 0; k < size; ++k) {
				mas3[i][j] += mas1[i][k] * mas2[k][j];
				
			}
		}
	}
	
	
	
	auto finish = std::chrono::steady_clock::now();
	std::chrono::duration<double> time = finish - start;

	return time.count();

}

void ReleaseMemory(double**& mas, int size) {
	for (int i = 0; i < size; ++i) {
		delete[] mas[i];
	}
	delete[] mas;

}


int main()
{
	int  n;
	
	std::cout << "Enter size of matrix: ";
	std::cin >> n;
	
	double** A;
	double** B;
	double** C;
	AllocateMemory(A, n);
	AllocateMemory(B, n);
	AllocateMemory(C, n);


	FillAB(A, B, n);

	FillZero(C, n);

	double assignment = MultiplyIJK(A, B, C, n);
	std::cout << "Time: " << assignment << "seconds" << std::endl;


	


	double middd = C[n / 2][n / 2];
	


	double Numflopspersec;
	std::cout <<"Middle Element: " << middd << std::endl;
	Numflopspersec = (2.0 * n * n * n /  assignment) / (1<<30);
	std::cout <<"Num of Gflops: "<< Numflopspersec << std::endl;
	

	ReleaseMemory(A, n);
	ReleaseMemory(B, n);
	ReleaseMemory(C, n);





	return 0;
}
```

```
#include <iostream>
#include <cmath>
#include <chrono>

void AllocateMemory(double*& mas, int size) {
	
	mas = new double [size];
	/*for (int i = 0; i < size; ++i) {
		mas[i] = new double[size];
	}*/

}
void FillAB(double*& mas1, double*& mas2, int size) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			mas1[i*size + j] = std::sin(i - j);
		}
	}
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			mas2[i*size + j] = 1.0 / (i + j + 1.0);
		}
	}
}
void FillZero(double*& mas, int size) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			mas[i*size + j] = 0;
		}
	}
}
double MultiplyIJK(double*& mas1, double*& mas2, double*& mas3, int size) {
	auto start = std::chrono::steady_clock::now();
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			for (int k = 0; k < size; ++k) {
				mas3[i*size + j] += mas1[i*size + k] * mas2[k*size + j];

			}
		}
	}



	auto finish = std::chrono::steady_clock::now();
	std::chrono::duration<double> time = finish - start;

	return time.count();

}
double MultiplyJKI(double*& mas1, double*& mas2, double*& mas3, int size) {
	auto start = std::chrono::steady_clock::now();
	for (int j = 0; j < size; ++j) {
		for (int k = 0; k < size; ++k) {
			for (int i = 0; i < size; ++i) {
				mas3[i * size + j] += mas1[i * size + k] * mas2[k * size + j];

			}
		}
	}



	auto finish = std::chrono::steady_clock::now();
	std::chrono::duration<double> time = finish - start;

	return time.count();

}
double MultiplyKIJ(double*& mas1, double*& mas2, double*& mas3, int size) {
	auto start = std::chrono::steady_clock::now();
	for (int k = 0; k < size; ++k) {
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				mas3[i * size + j] += mas1[i * size + k] * mas2[k * size + j];

			}
		}
	}



	auto finish = std::chrono::steady_clock::now();
	std::chrono::duration<double> time = finish - start;

	return time.count();

}
double MultiplyIKJ(double*& mas1, double*& mas2, double*& mas3, int size) {
	auto start = std::chrono::steady_clock::now();
	for (int i = 0; i < size; ++i) {
		for (int k = 0; k < size; ++k) {
			for (int j = 0; j < size; ++j) {
				mas3[i * size + j] += mas1[i * size + k] * mas2[k * size + j];

			}
		}
	}



	auto finish = std::chrono::steady_clock::now();
	std::chrono::duration<double> time = finish - start;

	return time.count();

}
double MultiplyJIK(double*& mas1, double*& mas2, double*& mas3, int size) {
	auto start = std::chrono::steady_clock::now();
	for (int j = 0; j < size; ++j) {
		for (int i = 0; i < size; ++i) {
			for (int k = 0; k < size; ++k) {
				mas3[i * size + j] += mas1[i * size + k] * mas2[k * size + j];

			}
		}
	}



	auto finish = std::chrono::steady_clock::now();
	std::chrono::duration<double> time = finish - start;

	return time.count();

}
double MultiplyKJI(double*& mas1, double*& mas2, double*& mas3, int size) {
	auto start = std::chrono::steady_clock::now();
	for (int k = 0; k < size; ++k) {
		for (int j = 0; j < size; ++j) {
			for (int i = 0; i < size; ++i) {
				mas3[i * size + j] += mas1[i * size + k] * mas2[k * size + j];

			}
		}
	}



	auto finish = std::chrono::steady_clock::now();
	std::chrono::duration<double> time = finish - start;

	return time.count();

}

void ReleaseMemory(double*& mas, int size) {
	/*for (int i = 0; i < size; ++i) {
		delete[] mas[i];
	}*/
	delete[] mas;

}


int main()
{
	int  n;

	std::cout << "Enter size of matrix: ";
	std::cin >> n;

	double* A;
	double* B;
	double* C;
	AllocateMemory(A, n);
	AllocateMemory(B, n);
	AllocateMemory(C, n);


	FillAB(A, B, n);

	FillZero(C, n);
	double assignment = MultiplyIJK(A, B, C, n);
	std::cout << "Time IJK: " << assignment << "seconds" << std::endl;
	double middd = C[n / 2][n / 2];
	double Numflopspersec;
	std::cout << "Middle Element: " << middd << std::endl;
	Numflopspersec = (2.0 * n * n * n / assignment) / (1 << 30);
	std::cout << "Num of Gflops: " << Numflopspersec << std::endl;




	FillZero(C, n);
	double assignment = MultiplyJKI(A, B, C, n);
	std::cout << "Time JKI: " << assignment << "seconds" << std::endl;
	double middd = C[n / 2][n / 2];
	double Numflopspersec;
	std::cout << "Middle Element: " << middd << std::endl;
	Numflopspersec = (2.0 * n * n * n / assignment) / (1 << 30);
	std::cout << "Num of Gflops: " << Numflopspersec << std::endl;





	FillZero(C, n);
	double assignment = MultiplyKIJ(A, B, C, n);
	std::cout << "Time KIJ: " << assignment << "seconds" << std::endl;
	double middd = C[n / 2][n / 2];
	double Numflopspersec;
	std::cout << "Middle Element: " << middd << std::endl;
	Numflopspersec = (2.0 * n * n * n / assignment) / (1 << 30);
	std::cout << "Num of Gflops: " << Numflopspersec << std::endl;

	FillZero(C, n);
	double assignment = MultiplyIKJ(A, B, C, n);
	std::cout << "Time IKJ: " << assignment << "seconds" << std::endl;
	double middd = C[n / 2][n / 2];
	double Numflopspersec;
	std::cout << "Middle Element: " << middd << std::endl;
	Numflopspersec = (2.0 * n * n * n / assignment) / (1 << 30);
	std::cout << "Num of Gflops: " << Numflopspersec << std::endl;


	FillZero(C, n);
	double assignment = MultiplyJIK(A, B, C, n);
	std::cout << "Time JIK: " << assignment << "seconds" << std::endl;
	double middd = C[n / 2][n / 2];
	double Numflopspersec;
	std::cout << "Middle Element: " << middd << std::endl;
	Numflopspersec = (2.0 * n * n * n / assignment) / (1 << 30);
	std::cout << "Num of Gflops: " << Numflopspersec << std::endl;




	FillZero(C, n);
	double assignment = MultiplyKJI(A, B, C, n);
	std::cout << "Time KJI: " << assignment << "seconds" << std::endl;
	double middd = C[n / 2][n / 2];
	double Numflopspersec;
	std::cout << "Middle Element: " << middd << std::endl;
	Numflopspersec = (2.0 * n * n * n / assignment) / (1 << 30);
	std::cout << "Num of Gflops: " << Numflopspersec << std::endl;


	ReleaseMemory(A, n);
	ReleaseMemory(B, n);
	ReleaseMemory(C, n);





	return 0;
}
```






```
#include <iostream>
#include <cmath>
int main()
{
    long double  a, b, c, x, D, x1, x2;
    std::cin >> a >> b >> c;
    D = b * b - 4 * a * c;
    if (a == 0 && b == 0 && c == 0) {
        std::cout << "X any value";
    }
    if (a == 0 && b == 0 && c != 0) {
        std::cout << "False";
    }
    if (a == 0 && b != 0 && c == 0) {
        x = 0;
        std::cout << x;
    }
    if (a == 0 && b!=0 && c!=0) {
        std::cout << "Not a quadratic equation";
        x = c / b;
        std::cout << x;
    }
    if (b == 0 && a!=0 && c!=0) {
        if ((a > 0 && c < 0) && (a<0 && c>0)) {
            x1 = sqrt(c / a);
            x2 = -sqrt(c / a);
            std::cout << x1 << " " << x2;
        }
        else {
            std::cout << "No solutions";
        }
    }
    if (a != 0 && b == 0 && c == 0) {
        x = 0;
        std::cout << x;
    }
    if (a != 0 && b != 0 && c == 0) {
        x1 = 0;
        x2 = -(b / a);
        std::cout << x1 << " " << x2;
    }




   
    if (D == 0)
    {
        x = (-b) / (2 * a);
        std::cout << x;
     }
    if (D > 0)
    {
        x1 = (-b + sqrt(D)) / (2 * a);
        x2 = (-b - sqrt(D)) / (2 * a);
        std::cout << x1 << " " << x2 << std::endl;
    }
    if (D < 0) {
        std::cout << "No solutions";
    }
    
   
}
```

```
#include <iostream>
#include <random> 
#include <algorithm>
int main()
{ 
	int n, m;
	std::cout << "Enter lines: " << std::endl;
	std::cin >> m;
	std::cout << "Enter columns: " << std::endl;
	std::cin >> n;
	int** mas;
	mas = new int* [m];
	for (int i = 0; i < m; ++i) {
		mas[i] = new int[n];
	}
	int min, max;
	
	
	std::cout << "If you want your own massive, please, enter : '1' " << std::endl;
	std::cout << "If you want random massive, please, enter : '2' " << std::endl;
	int choose;
	std::cin >> choose;
	if (choose == 1) {
		for (int i = 0; i < m;++i) {
			for (int j = 0; j < n; ++j) {
				std::cin >> mas[i][j];
				
				
			}
		}
	}
	if (choose == 2) {
		for (int i = 0; i < m;++i) {
			for (int j = 0; j < n; ++j) {
				mas[i][j] = rand();
			}
		}
	}

	for (int i = 0; i < m;++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i][j] << " ";
		}
		std::cout << std::endl;
	}
	min = mas[0][0];
	max = mas[0][0];
	//на строчки
	int cnt = 0;
	int cntstring=0;
	
	
	for (int i = 0; i < m; ++i) {
		cnt = 0;
		
		for (int j = 0; j < n; ++j) {
			cnt += mas[i][j];
			
		}
		
		if (cnt < min) {
			min = cnt;
			cntstring = i;
		}
	}
	
	/*++cntstring;*/
	//std::cout << cntstring << std::endl;


	for (int i = cntstring; i < (m-1); ++i) {
		for (int j = 0; j < n; ++j) {
			std::swap(mas[i][j], mas[i + 1][j]);
		}
	}


// столбцы
	int cnt2=0;
	int cntcolumns=0;
	int max = 0;
	for (int j = 0; j < n;++j) {
		cnt2 = 0;
		for (int i = 0; i < m; ++i) {
			cnt2 += mas[i][j];
		}
		if (cnt2 > max) {
			max = cnt2;
			cntcolumns = j;
		}
	}
	/*std::cout << cntcolumns<<std::endl;*/

	for (int j = cntcolumns; j<(n-1) ; ++j) {
		for (int i = 0; i < n; ++i) {
			std::swap(mas[i][j], mas[i][j+1]);
		}
	}





	for (int i = 0; i < m;++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << mas[i][j] << " ";
		}
		std::cout << std::endl;
	}





	/*for (int i = 0; i < m; ++i) {
		delete[] mas[m];
	}*/
	delete[] mas;



	return 0;
}
```
