# First-repository
My first repository
```
//#include <iostream>
//#include <cmath>
//#include <chrono>
//double middleelement;
//void AllocateMemory(double**& mas, int size) {
//	mas = new double* [size];
//	for (int i = 0; i < size; ++i) {
//		mas[i] = new double[size];
//	}
//
//}
//void FillAB(double**& mas1, double**& mas2, int size) {
//	for (int i = 0; i < size; ++i) {
//		for (int j = 0; j < size; ++j) {
//			mas1[i][j] = std::sin(i - j);
//		}
//	}
//	for (int i = 0; i < size; ++i) {
//		for (int j = 0; j < size; ++j) {
//			mas2[i][j] = 1.0 / (i + j + 1.0);
//		}
//	}
//}
//void FillZero(double**& mas, int size) {
//	for (int i = 0; i < size; ++i) {
//		for (int j = 0; j < size; ++j) {
//			mas[i][j] = 0;
//		}
//	}
//}
//double MultiplyIJK(double**& mas1, double**& mas2, double**& mas3, int size) {
//	auto start = std::chrono::steady_clock::now();
//	for (int i = 0; i < size; ++i) {
//		for (int j = 0; j < size; ++j) {
//			for (int k = 0; k < size; ++k) {
//				mas3[i][j] += mas1[i][k] * mas2[k][j];
//				
//			}
//		}
//	}
//	
//	if (size % 2 != 0) {
//		middleelement = mas3[size / 2][size / 2];
//	}
//	else {
//		middleelement = mas3[(size / 2) - 1][(size / 2) - 1];
//	}
//	
//	auto finish = std::chrono::steady_clock::now();
//	std::chrono::duration<double> time = finish - start;
//
//	return time.count();
//
//}
//
//void ReleaseMemory(double**& mas, int size) {
//	for (int i = 0; i < size; ++i) {
//		delete[] mas[i];
//	}
//	delete[] mas;
//
//}
////double MiddleElement(double**& mas3, int size) {
////	double middleelement;
////	if (size % 2 != 0) {
////		middleelement = mas3[size / 2][size / 2];
////	}
////	else {
////		middleelement = mas3[(size / 2) - 1][(size / 2) - 1];
////	}
////	return middleelement;
////}
//
//int main()
//{
//	int  n, helpvariable1;
//	int Numflops, Numflopspersec;
//	std::cout << "Enter sixe of matrix: ";
//	std::cin >> n;
//	Numflops = n * n * n;
//	helpvariable1 = 1ULL << 30;
//	double** A;
//	double** B;
//	double** C;
//	AllocateMemory(A, n);
//	AllocateMemory(B, n);
//	AllocateMemory(C, n);
//	FillAB(A, B, n);
//	FillZero(C, n);
//	double assignment = MultiplyIJK(A, B, C, n);
//	std::cout << "Time: " << assignment << "seconds" << std::endl;
//	ReleaseMemory(A, n);
//	ReleaseMemory(B, n);
//	ReleaseMemory(C, n);
//	/*double middd = MiddleElement(C, n);
//	std::cout << middd<<std::endl;*/
//
//	std::cout <<"Middle Element: " << middleelement << std::endl;
//	Numflopspersec = (Numflops /  assignment); // helpvariable1;
//	std::cout <<"Num of flops: "<< Numflopspersec << std::endl;
//	
//
//
//
//
//
//
//
//	return 0;
//}
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
#include <cstdlib>
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
	//на строчки
	int cnt = 0;
	int cntstring=0;
	int min=2147483647;
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
	
	++cntstring;
	//std::cout << cntstring << std::endl;


	for (int i = (cntstring-1); i < (m-1); ++i) {
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
	++cntcolumns;
	/*std::cout << cntcolumns<<std::endl;*/

	for (int j = (cntcolumns - 1); j<(n-1) ; ++j) {
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





	for (int i = 0; i < m; ++i) {
		delete[] mas[n];
	}
	delete[] mas[m];



	return 0;
}
```
