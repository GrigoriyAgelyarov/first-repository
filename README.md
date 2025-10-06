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
/* 1)/#include <iostream>

int main()
{
    const double aperc = 0.1, bperc = 0.05;
    double abal, bbal, peep;
    int cnt = 0;
    std::cin >> abal >> bbal;
    peep = abal;
    if (bbal >= abal) {
        std::cout << cnt;
    }
    while (abal > bbal) {
        cnt++;
        abal += aperc * peep;
        bbal += bperc * bbal;
        std::cout << "Year: " << cnt << " Anna's balance: " << abal << " Boris's balance: " << bbal << std::endl;
    }





}
*/
/* 2)#include <iostream>

int main(){
    const double perc1 = 0.0, perc2 = 0.1, perc3 = 0.15, perc4 = 0.2, gate1 = 5000, gate2 = 10000, gate3=20000, gate4=35000;
    double Sum;
    std::cin >> Sum;
    double cnt = 0;
    if (Sum < 0) {

    }
    while (Sum > 0) {
        if (Sum <= gate1) {
            Sum = Sum - gate1;
            std::cout << cnt;
        }
        if (Sum <= gate2 && Sum > gate1) {
            cnt = (gate1 * perc1) + ((gate2 - gate1) * perc2);
            Sum = Sum - gate1 - gate2;
            std::cout << cnt;
        }
        if (Sum <= gate3 && Sum > gate2) {
            cnt = (gate1 * perc1) + ((gate2) * perc2) + ((Sum - gate1 - gate2)*perc3);
            Sum = Sum - gate1 - gate2 - gate3;
            std::cout << cnt;
        }
        if (Sum <= gate4 && Sum > gate3) {
            cnt = (gate1 * perc1) + ((gate2) * perc2) + ((gate3) * perc3) + ((Sum - gate1 - gate2 - gate3) * perc4);
            Sum = Sum - gate1 - gate2 - gate3 - gate4;
            std::cout << cnt;
        }
        if (Sum > gate4) {
            cnt = (gate1 * perc1) + ((gate2) * perc2) + ((gate3) * perc3) + ((gate4) * perc4) + ((Sum - gate1 - gate2 - gate3 - gate4) * perc4);
            Sum = Sum - gate1 - gate2 - gate3 - gate4;
            std::cout << cnt;
        }

    }
}
*/

#include <iostream>

int main() {
    int n;
    double sum = 0, middle, g;
    std::cin >> n;

    for (int i = 1; i <= n; ++i) {
        
        std::cin >> g;
        sum += g;
        middle = sum / i;  
        std::cout << "Arithmetic mean: " << middle << std::endl;
    }

    return 0;
}



```
