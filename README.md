# First-repository
My first repository
```
#include <iostream>
#include <string>
#include <vector>
#include <list>
class LiverpoolFC{
	int goals;
	std::string name;
};
class HashTable {
private:
	struct Entry {
		std::string key;
		LiverpoolFC value;
		bool isDeleted;
	};
	std::vector<Entry*> table;
	size_t capacity;
	size_t size_;
	size_t hashString(const std::string& s) const {
		const long long p = 31;
		const long long M = 1e9 + 7;
		long long h = 0;
		for (size_t i = 0; i < s.length(); ++i) {
			h = (h * p + s[i]) & M;
		}
		return static_cast<size_t>(h);

	}
	size_t hashInt(int g) const {
		return ((g << 5) ^ (g >> 3)) & 0x7FFFFFFF;
	}
	size_t getIndex(const std::string& key, int goals) const {
		size_t strHash = hashString(key);
		size_t intHash = hashString(key);
		return(strHash ^ intHash) % capacity;
	}
};


int main()
{
   
}


```
```
#include <iostream>
#include <vector>


int main() {
	int n;
	int tmp = 1000000007;
	std::cin >> n;
	int help;
	std::vector<int> a;
	for (int i = 0; i < n; ++i) {
		std::cin >> help;
		a.push_back(help);
	}
	int l = 0;

	long long sum = 0;
	

	while ((l + 2) < n) {
    sum += (a[l] * a[l+1]*a[l+2])%tmp;
    
		
		++l;
		


	}
	
	std::cout << sum;
	return 0;
}
```
```
#include <iostream>
#include <vector>
#include <string>

struct Discipline {
    std::string name;
    int mark;
};

struct Student {
    std::string fio;
    std::vector<Discipline> disciplines;
};

int main() {
    std::cout << "Enter number of students: ";
    int n;
    std::cin >> n;
    std::string help;
    std::getline(std::cin, help);
    std::vector<Student> group;

    for (int i = 0; i < n; ++i) {
        Student student;
        

        std::cout << "Enter FIO of the student: " << std::endl;
        std::getline(std::cin, student.fio); 
        int n2;
        std::cout << "Enter number of disciplines: ";
        std::cin >> n2;
        std::getline(std::cin, help); 

        if (n2 <= 3) {
            std::cout << "Less than 4 disciplines :(\n";
            continue;
        }

        for (int j = 0; j < n2; ++j) {
            Discipline discipline;
            std::cout << "Enter the name of discipline: " << std::endl;
            std::getline(std::cin, discipline.name); 

            std::cout << "Enter mark: ";
            std::cin >> discipline.mark;
            std::getline(std::cin, help); 

            student.disciplines.push_back(discipline);
        }

        group.push_back(student);
    }

    int cnt = 0;

    for (int i = 0; i < group.size(); ++i) {
        Student student = group[i];
        for (int j = 0; j < student.disciplines.size(); ++j) {
            Discipline discipline = student.disciplines[j];

            if (discipline.mark == 2) {
                cnt++;
                break;
            }
        }
    }

    std::cout << "Number of students with '2': " << cnt << std::endl;
    return 0;
}
```
```
#include <iostream>
#include <vector>
int main() {
    int n, dist;
    std::cin >> n >> dist;
    std::vector<int> a;
    int hlp;
    for (int i = 0; i < n;++i) {
        std::cin >> hlp;
        a.push_back(hlp);
    }
    int l = 0;
    int r = 1;
    long long cnt = 0;
    while ((a[r] - a[l]) <= dist) {
        r++;
        while (r < n-1) {
            if ((a[r] - a[l]) > dist) {
                l++;
            }
            else {
                cnt += n - r;
            }
        }
    }
    std::cout << cnt;
}
```
```
//#include <iostream>
//#include <cmath>
//int main() {
//    double a, b, c;
//    std::cin >> a >> b >> c;
//    double eps = 2.2e-16;
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
//            std::cout << x;
//        }
//        return 0;
//    }
//    double D = b * b - 4 * a * c;
//    double help1 = b * b;
//    double help2 = 4 * fabs(a * c);
//    double maxD = help1;
//    if (help2 > maxD) {
//        maxD = help2;
//    }
//    if (D > eps * maxD) {
//        double x1, x2;
//        x1 = (-b - sqrt(D)) / (2 * a);
//        x2 = c / (x1 * a);
//        std::cout << x1 << " " << x2;
//    }
//    else if (fabs(D) <= eps * maxD) {
//        double x = -b / (2 * a);
//        std::cout << x;
//    }
//    else {
//        std::cout << "No solutions";
//    }
//    return 0;
//}

```
```
#include <iostream>
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
	for (int i = 0; i < m; ++i) {
		for (int j = cntcolumns; j < (n - 1); ++j){
			int* ptr1 = &mas[i][j];
			int* ptr2 = &mas[i][j + 1];
			int temp = *ptr1;
			*ptr1 = *ptr2;
			*ptr2 = temp;
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
