# First-repository
My first repository
```
#include <iostream>

class vector {
	double x;
	double y;
	double z;
public:
	vector operator+(const vector& b);

	vector operator-(const vector& b);

	void set(double x_, double y_, double z_) { x = x_; y = y_; z = z_; };

	friend std::ostream& operator<<(std::ostream& str, const vector& c); //мы типо не меняем значения ну и типо доступ к приватным x y z

	vector operator++(int); // int если постфиксный декремент 
	vector operator++();
	vector operator--(int);
	vector operator--();

	double& operator[](int i);
	const double& operator[](int i) const;
	vector operator+();
	vector operator-();
	vector operator*(const vector& b);
	vector operator/(const vector& b);
	vector operator*(double i);
	friend vector operator*(double i, const vector& v); //???
	vector operator/(double i);
	vector operator*=(double i);
	vector operator/=(double i);
	double operator&(const vector& b);
	double length() { return sqrt(x * x + y * y + z * z); };
	vector operator^(const vector& b);
	vector normalize();
	vector normalize(double i);
	vector unit();
	vector unit(double i);
	vector proj(vector& b);
	vector rotate(vector& b, double fi);
};
```
```
#include "vec.h"
#include <cmath>
vector vector::operator+ (const vector& b) { // слева объект справа то, что принимаем
	vector c;
	c.x = x + b.x;
	c.y = y + b.y;
	c.z = z + b.z;
	return c;
}
std::ostream& operator<<(std::ostream& str, const vector& c) { // тут типо слева когда выводим идет строка, а справа вектор как раз, поэтому такая запись 
	str << c.x << " " << c.y << " " << c.z;
	return str;
}
vector vector::operator-(const vector& b) {
	vector c;
	c.x = x - b.x;
	c.y = y - b.y;
	c.z = z - b.z;
	return c;
}
vector vector::operator++(int) {
	vector c= *this;
	this->x += 1;
	this->y += 1;
	this->z += 1;
	return c;
}
vector vector::operator++() {
	x += 1;
	y += 1;
	z += 1;
	return *this;
	/*this->x += 1;
	this->y += 1;
	this->z += 1;
	return *this;*/ // разница способов
}
vector vector::operator--(int) {
	vector c = *this;
	this->x -= 1;
	this->y -= 1;
	this->z -= 1;
	return c;
}
vector vector::operator--() {
	this->x -= 1;
	this->y -= 1;
	this->z -= 1;
	return *this;

}
double& vector::operator[](int i) { //(возможно изменение)
	switch (i) {
	case 0:
		return x;
	case 1:
		return y;
	case 2:
		return z;
	}
}
const double& vector::operator[](int i) const { // 2 конст - мы не меняем какие-то значения в функции, а 1 конст нельзя изменять элементы(чтение)
	switch (i) {
	case 0:
		return x;
	case 1:
		return y;
	case 2:
		return z;
	}
}
vector vector::operator+() {
	return *this;
}
vector vector::operator-() {
	x *= (-1);
	y *= (-1);
	z *= (-1);
	return *this;
}
vector vector::operator*(const vector& b) {
	vector c;
	c.x = x * b.x;
	c.y = y * b.y;
	c.z = z * b.z;
	return c;
}
vector vector::operator/(const vector& b) {
	vector c;
	c.x = x / b.x;
	c.y = y / b.y;
	c.z = z / b.z;
	return c;
}
vector vector::operator*(double i) {
	vector c;
	c.x = x * i;
	c.y = y * i;
	c.z = z * i;
	return c;
}
vector operator*(double i, const vector& b) {
	vector c;
	c.x = i * b.x;
	c.y = i * b.y;
	c.z = i * b.z;
	return c;
}
vector vector::operator/(double i) {
	vector c;
	c.x = x / i;
	c.y = y / i;
	c.z = z / i;
	return c;
}
vector vector::operator*=(double i) {
	x *= i;
	y *= i;
	z *= i;
	return *this;
}
vector vector::operator/=(double i) {
	x /= i;
	y /= i;
	z /= i;
	return *this;
}
double vector::operator&(const vector& b) { //без ссылки, так как c временная переменная
	double c;
	c = x * b.x + y * b.y + z * b.z;
	return c;
}
vector vector::operator^(const vector& b) { // ?????? вернуть вектор или длину вектора? "используем eps" теперь во всех задачах с double?
	vector c;
	double eps = 2.2e-16;
	c.x = y * b.z - z * b.y;
	c.y = - (x * b.z - z * b.x);
	c.z = x * b.y - y * b.x;
	if (fabs(c.x) < eps) {
		c.x = 0;
	}
	if (fabs(c.y) < eps) {
		c.y = 0;
	}
	if (fabs(c.z) < eps) {
		c.z = 0;
	}
	return c;
}
vector vector::normalize() {
	double w = (*this).length();
	this->x = x / w;
	this->y = y / w;
	this->z = z / w;
	return *this;
}
vector vector::normalize(double i) {
	double w = (*this).length();
	this->x = (x / w) * i;
	this->y = (y / w) * i;
	this->z = (z / w) * i;
	return *this;
	
}
vector vector::unit() {
	vector c = *this;
	c.normalize();
	return c;
}
vector vector::unit(double i) {
	vector c = *this;
	c.normalize(i);
	return c;
}
vector vector::proj(vector& b) {
	vector c = *this;
	c = ((*this & b) / (b & b)) * b;
	return c;
}
vector vector::rotate(vector& b, double fi) {

}
```
```
#include <iostream>
#include "vec.h"

int main()
{
    vector a, b, c, d;
    double w;
    a.set(1.0, 2.0, 3.0);
    b.set(0.1, 0.2, 0.3);

    c = a + b;
    std::cout << c << std::endl;

    c = a - b;
    std::cout << c << std::endl;

    c = a++;
    std::cout << a << " " << c << std::endl;

    a.set(1.0, 2.0, 3.0);
    c = ++a;
    std::cout << a << " " << c << std::endl;

    a.set(1.0, 2.0, 3.0);
    c = a--;
    std::cout << a << " " << c << std::endl;

    a.set(1.0, 2.0, 3.0);
    //a[0] = 3.2;
    c = --a;
    std::cout << a << " " << c << std::endl;

    a.set(1.0, -2.0, 3.0);
    /*c = +a;
    d = -a;
    std::cout << c << " " << d << std::endl;*/
    std::cout << (+a) << " " << (-a) << std::endl;

    a.set(1.0, 2.0, 3.0);
    b.set(0.1, 0.2, 0.3);

    c = a * b;
    std::cout << c << std::endl;

    c = (a / b);
    std::cout << c << std::endl;

    
    c = a * 2.0;
    std::cout << c << std::endl;
 
    c = 2.0 * a;
    std::cout << c << std::endl;

    c = a / 2.0;
    std::cout << c << std::endl;

    a *= 2.0;
    std::cout << a << std::endl;

    a /= 2.0;
    std::cout << a << std::endl;

    w = a & b;
    std::cout << w << std::endl;

    std::cout << a.length() << std::endl;

    a.set(1.0, 2.0, 3.0);
    b.set(0.5, 1.25, 1.5);
    c = a ^ b;
    std::cout << c << std::endl;
    b.set(0.1, 0.2, 0.3);

    std::cout << (a.normalize()) << std::endl;

    std::cout << (a.normalize(2)) << std::endl;

    a.set(1.0, 2.0, 3.0);
    std::cout << (a.unit()) << " " << a << std::endl;

    std::cout << (a.unit(2)) << " " << a << std::endl;

    a.set(1.0, 2.0, 3.0);
    b.set(3.0, 1.5, 2.0);
    c = a.proj(b);
    std::cout << a.proj(b) << " " << a << " " << b << " " << std::endl;


}
```
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
