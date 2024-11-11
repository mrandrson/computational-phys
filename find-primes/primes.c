#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int gcd(int n1, int n2) {
    int a, b, r;
    if (n1 > n2) {
        a = n1;
        b = n2;
    } else {
        a = n2;
        b = n1;
    }
   while (b != 0) {
      r = a % b;
      a = b;
      b = r;
   }
   return a;
}
