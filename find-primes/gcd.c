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

int main() {
    int n1, n2;
    printf("Enter two integers: ");
    scanf("%d %d", &n1, &n2);
    printf("GCD of %d and %d is %d\n", n1, n2, gcd(n1, n2));
    return 0;
}
