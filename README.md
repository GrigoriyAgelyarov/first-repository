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
