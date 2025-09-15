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
