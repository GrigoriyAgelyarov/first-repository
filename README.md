# First-repository
My first repository



```#include <iostream>

int main()
{
    int a, b, c;
    std::cin >> a >> b >> c;
    if (a + b > c && a + c > b && b + c > a) {
        std::cout << "YES";
    }
    else {
        std::cout << "NO";
    }
    return 0;
}
```

```#include <iostream>

int main() {
    int a, b, c;
    std::cin >> a >> b >> c;
    if ((b < a && a < c) || (c < a && a < b)) {
        std::cout << a;
    }
    if ((a < b && b < c) || (a > b && b > c)) {
        std::cout << b;
    }
    else {
        std::cout << c;
    }
return 0;
}
```
```
#include <iostream>
#include <cmath>
int main()
{
    double a, b, c, x, D, x1, x2;
    std::cin >> a >> b >> c;
    D = b * b - 4 * a * c;
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
        std::cout << "No korney";
    }
    if (a == 0) {
        std::cout << "Ne kvvadratnoye uravnenie";
    }
}
```
```
int main() {
    const double perc1 = 0.0, perc2 = 0.1, perc3 = 0.15, perc4 = 0.2, gate1 = 5000.0, gate2 = 15000.0, gate3 = 35000.0;
    double Sum;
    double cnt = 0, loan1=0, loan2 =0, loan3 =0, loan4=0;

   
    do {
        std::cin >> Sum;
        
        
        
        if (Sum <= gate1 && Sum>=0) {
            loan1 = Sum * perc1;
            std::cout << loan1;
        }
        if (Sum > gate1 && Sum <= gate2) {
            loan2 = (Sum - gate1) * perc2 + perc1*gate1;
            std::cout << loan2;
        }
        if (Sum > gate2 && Sum <= gate3) {
            loan3 = ((Sum - gate2) * perc3) + (gate2-gate1)*perc2 + perc1 * gate1;
            std::cout << loan3;
        }
        
        if (Sum > gate3) {
            loan4 = ((Sum - gate3) * perc4) + (gate3-gate2)*perc3 + (gate2-gate1)*perc2 + gate1*perc1;
            std::cout << loan4;
        }





        if (Sum < 0) {
            break;
        }

    } while (Sum > 0);
    

        

      return 0;
}
```
```
#include <iostream>

int main()
{
	int n, m;
	std::cout << "Enter columns: " << std::endl;
	std::cin >> m;
	std::cout << "Enter lines: "<<std::endl;
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

	}

	if (choose == 2) {

	}













	for (int i = 0; i < m; ++i) {
		delete[] mas[n];
	}
	delete[] mas[m];
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
	std::cout << assignment;
	ReleaseMemory(A, n);
	ReleaseMemory(B, n);
	ReleaseMemory(C, n);
	return 0;
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
