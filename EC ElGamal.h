// Elliptic Curve ElGamal.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
using namespace std;
#define ui unsigned
#define llint long long int



// mod operator
llint mod(llint a, llint b) {
	llint res = round(remainder((double)a, b));
	while (res <= 0)
		res += b;
	return res;
}

// Utility function to do modular exponentiation.
// It returns (x^y) % p
llint power(llint x, llint y, llint p) {
	llint res = 1; // Initialize result
	x = mod(x, p); // Update x if it is more than or equal to p

	y = mod(y, p - 1); // we will exponentiate in mod phi(p) = p-1

	while (y > 0)
	{
		// If y is odd, multiply x with result
		if (y & 1)
			res = mod(res * x, p);

		// y must be even now
		y = y >> 1; // y = y/2
		x = mod(x * x, p);
	}
	return res;
}

// This function is called for all k trials. 
// It returns false if n is composite and returns true if n is probably prime.
// d is an odd number such that d*2<sup>r</sup> = n-1 for some r >= 1
bool MR_Test(llint d, llint n)
{
	// Pick a random number in [2..n-2]
	// Corner cases make sure that n > 4
	llint a = 2 + rand() % (n - 4);

	// Compute a^d % n
	llint x = power(a, d, n);

	if (x == 1 || x == n - 1)
		return true;

	// Keep squaring x while one of the following doesn't
	// happen
	// (i) d does not reach n-1
	// (ii) (x^2) % n is not 1
	// (iii) (x^2) % n is not n-1
	while (d != n - 1)
	{
		x = (x * x) % n;
		d *= 2;

		if (x == 1)	 return false;
		if (x == n - 1) return true;
	}

	// Return composite
	return false;
}

// It returns false if n is composite and returns true if n
// is probably prime. k is an input parameter that determines
// accuracy level. Higher value of k indicates more accuracy.
bool isPrime(llint n, llint k = 5)
{
	// Corner cases
	if (n <= 1 || n == 4) return false;
	if (n <= 3) return true;

	// Find r such that n = 2^d * r + 1 for some r >= 1
	llint d = n - 1;
	while (d % 2 == 0)
		d /= 2;

	// Iterate given number of 'k' times
	for (llint i = 0; i < k; i++)
		if (!MR_Test(d, n))
			return false;

	return true;
}

// Function to generate a size-bit prime ( 2^(size-1) <p < 2^size), 20-bit by default
// uses the Miller-Rabin test with k passes, 5 passes by default
bool genPrime(llint& p, llint size = 20, llint k = 5) {
	srand(time(NULL));
	size -= 1;
	p = pow(2, size);

	llint SIZE = pow(2, size);
	llint MAX = pow(2, size + 2);

	llint r = rand() % (SIZE / 2);
	if (r % 2 == 0)
		r++;
	p = SIZE + r;//p is in the range [2^size, 2^(size+1)
	//cout << SIZE << " " << r << " " << p << " " << SIZE * 2 << endl;

	while (!isPrime(p, k) && p < MAX) {
		p += 2;
	}
	if (isPrime(p, 5)) {
		return true;
	}
	return false;
}

// Function for extended Euclidean Algorithm
llint gcdExt(llint a, llint b, llint & x, llint & y){	
	// Base Case
	if (a == 0)
	{
		x = 0;
		y = 1;
		return b;
	}

	llint  x_, y_; // To store results of recursive call
	llint gcd = gcdExt(b % a, a, x_, y_);

	// Update x and y using results of the recursive call
	x = y_ - (b / a) * x_;
	y = x_;
	return gcd;
}

// Function to invert an element in mod p
llint inv(llint a, llint p) {
	//llint a = x;
	//llint b = p;
	llint x, y;
	if (a < 0)
		a = a + 10 * p;
	if (p < 0)
		p = -p;
	gcdExt(a, p, x, y);
	return mod(x,p);
}// verified!


llint list_EC_points(llint a, llint b, llint p) {
	llint num_el = 0;
	for (llint i = 0; i < p; i++) {
		for (llint j = 0; j < p; j++) {
			if ( mod(power(j,2,p) , p) == mod( power(i, 3, p) + mod(a * i, p) + b , p) ) {
				//cout << "(" << i << "," << j << ")" << endl;
				num_el++;
			}
		}
	}
	cout << "The number of elements is " << num_el + 1 << endl;
	return num_el + 1;
}

// This function tests if (x,y) is on the elliptic curve y^2 = x^3 + ax +b (mod p)
llint is_on_EC(llint x, llint y, llint a,llint b, llint p){
	//return ((y * y) % p == (x * x * x + a * x + b) % p);
	return mod(power(y, 2, p),p) == mod( ( power(x, 3, p) + (a * x ) %p + b ) + 10* p , p);
}


llint verify_parameters(llint a, llint b, llint p) {
	llint D = mod( 4 * power(a,3,p)  + (27 * power(b,2,p) % p) , p);
	//cout << 4 * a * a * a + 27 * b * b<<endl;
	return (D != 0);
}


llint order(llint x, llint y, llint a, llint b, llint p);
// Create p, a, b, and a (non-) generator G=(gx,gy) of order ord_gen
// p = 3 (mod 4) required for easy encryption/decryption, i.e. for the encode/decode functions
void generate_parameters(llint& a, llint& b, llint& p, llint& gx, llint& gy, llint& ord_gen) {
	
	bool found = false;

	do {
		genPrime(p, 20, 5);
	} while (p % 4 != 3);
	
	// genPrime(p);
	//p = 1979;

	llint Hasse_min_SIZE, bit_size_of_gen, min_gen_order;

	while (!found) {
		srand(time(NULL));
		a = mod(rand(), p);

		srand(time(NULL) + 1);
		gx = 1 + rand() % 5;

		srand(time(NULL) + 2);
		gy = 1 + rand() % 10;
		//a = 3, gx = 4, gy = 11, p = 8831;
		b = mod(power(gy, 2, p) - power(gx, 3, p) - a * gx, p);

		Hasse_min_SIZE = p + 1 - 2 * sqrt(p);
		bit_size_of_gen = log2(Hasse_min_SIZE) - 2;
		min_gen_order = pow(2, bit_size_of_gen);
		ord_gen = order(gx, gy, a, b, p);
		//cout << "Trying\n" << Hasse_min_SIZE << " " << bit_size_of_gen << " " << min_gen_order << " " << ord_gen << endl;

		if (ord_gen >= min_gen_order)
			found = true;
	}
	
	printf("p=%lld\na=%lld\nb=%lld\nG=(%lld,%lld) with order %lld\n",p,a,b,gx,gy,ord_gen);
	//printf("E = EllipticCurve(GF(%lld),[0,0,0,%lld,%lld])\nG = E(%lld , %lld)\n", p, a, b, gx, gy, ord_gen);
}

 
void dbl(llint x, llint y, llint& rx, llint& ry, llint a, llint b, llint p) {
	/*Takes the x and y and the parameters a,b,p 
	doubles the point (x,y) in y^2 = x^3 + ax + b (mod p), and returns the output in x and y*/
	llint idy = inv(2 * y, p);
	//cout << idy << endl;
	llint temp;
	temp = (3 * power(x,2,p) + a) * idy;
	rx = mod( power(temp, 2, p) - 2 * x , p) ;
	ry = ( (temp * (x- rx) - y ) % p + p) % p;
}// verified!

void add(llint x1, llint y1, llint x2, llint y2 , llint& rx, llint& ry, llint a, llint b, llint p) {
	//printf("Adding (%d,%d) and (%d,%d)\n", x1, y1, x2, y2);
	if ((x1 %p == 0) && (y1 % p == 0) ){
		// the case where (x1,y1) = inf = zero
		rx = x2; ry = y2;
		return;
	}// verified!

	if ((x2 % p == 0) && (y2 % p == 0)) {
		// the case where (x2,y2) = inf = zero
		rx = x1; ry = y1;
		return;
	}// verified!

	if ((mod(x1 , p) == mod(x2 , p) )&& (mod(y1, p) == mod(y2, p))) {
		// the case of doubling
		dbl(x1, y1, rx, ry, a, b, p);
		return;
	}// verified!
	
	if ((mod(x1, p) == mod(x2, p)) && ( (y1+y2) % p==0) ) {
		// the case of adding (x1,y1) and -(x1,y1) = (x1, -y1) 
		// we should return inf = zero
		rx = 0; ry = 0;
		return;
	}// verified!

	//point addition case when P != +/- Q
	//cout << inv(x2 - x1, p) << endl;
	llint lambda;
	//if (x2 > x1) {
	//	lambda = (((y2 - y1) * inv(x2 - x1, p)) % p + 10 * p) % p;
	//}
	//else {
	//	lambda = (((y1 - y2) * inv(x1 - x2, p)) % p + 10 * p) % p;
	//}
	lambda = (((y2 - y1) * inv(x2 - x1, p)) % p + 10 * p) % p;
	rx = ( (lambda * lambda - x1 - x2 ) % p + 10* p ) %p ;
	ry = ( (lambda * (x1 - rx) - y1 ) %p + 10*p ) %p;
}


void k_times(llint k, llint x, llint y, llint& rx, llint& ry, llint a, llint b, llint p) {
	llint sign = 1;
	if (k == 0) {
		rx = 0; ry = 0;
		//cout << "k=0" << endl;
		return;
	}

	if (k < 0) {
		sign = -1;
		k = -k;
		//cout << "sign(k)=-1 and k=" << k<< endl;
	}
	//cout << "---"<< k << "---" << endl;
	// at this point k=abs(k)>0 and sign holds the sign of k
	// get binary representation of k
	llint size=0;
	llint bin_k[64] = { 0 };
	while (k > 0) {
		//cout << k <<" ";
		bin_k[size] = k % 2;
		size++;
		k= k >> 1;
	}
	//cout << endl;
	//for (llint i=size-1; i>= 0; i--)
	//	cout << bin_k[i] << " ";
	//cout << endl;
	// at this point, bin_k holds "size" elements which are the binary representation of k

	llint x_ = x, y_=y;
	//cout << "The size is " << size << endl;
	for (llint i = size - 2; i >= 0; i--) {
		//cout << "i=" << i << " and doubling. ";
		dbl(x_, y_, x_, y_, a, b, p);
		if (bin_k[i] == 1) {
			//cout << "bin_k[" << i << "]=" << bin_k[i] << " hence Q=Q+P";
			add(x, y, x_, y_, x_, y_, a, b, p);
		}
		//cout << endl;
	}
	rx = x_; 
	ry = y_;
	if (sign == -1)
		ry = (( - y_) % p +p )%p;
	
	return;
}

llint order(llint x, llint y, llint a, llint b, llint p) {
	if (x == 0 && y == 0)
		return 1;
	//else {
	//	llint x2, y2;
	//	llint sqrt_p = sqrt(p);
	//	cout << sqrt_p << endl;
	//	for (llint k = 2; k<p+3+2*sqrt_p; k++) {
	//		k_times(k, x, y, x2, y2, a, b, p);
	//		//printf("%d(%d,%d) = (%d,%d)\n", k, x, y, x2, y2);
	//		if (x2 == 0 && y2 == 0)
	//			return k;
	//	}
	//}
	llint k = 1;
	llint x2=x, y2=y;
	while (!(x2 == 0 && y2 == 0)) {
		//printf("%d(%d,%d) = (%d,%d)\n", k, x, y, x2, y2);
		add(x, y, x2, y2, x2, y2, a, b, p);
		k++;
	}
	return k;
}

// This encodes the message, but the message number has to satisfy (m+1)K<p
llint encode(llint m, llint &mx, llint &my, llint a, llint b, llint p, llint K=10) {
	// m: message represented as a number 
	// K: failure rate is smaller than 1/2^K 

	if ((m + 1) * K >= p ) {
		//This is the case NOT-explained in the book
		cout << "Problem encoding\n";
		exit(4130);
		return 0; 
	}


	llint x, rhs, sr;

	for (llint j = 0; j < K; j++) {
		x = mod(m * K + j ,p);
		//rhs = ( ((x * x * x) % p + a * x + b) + 10 * p ) % p;
		rhs = mod (power(x, 3, p) + mod(a * x, p) + b, p);
		sr = power(rhs, (p + 1) / 4, p);
		if (power(sr,2,p) == rhs) {
			mx = x;
			my = sr;
			return 1;
		}
	}
	return 0;
}

llint decode(llint& m, llint mx, llint my, llint a, llint b, llint p, llint K = 10) {
	m = mx / K;
	return m;
}

void text_to_number(llint& res, char* text, unsigned len) {
	//// Encode text into mpz_t so that the leftmost character is the most significant
	res = 0;

	for (unsigned i = 0; i < len; i++) {
		res = res + (unsigned) text[i];
		if (i < len - 1) {
			res = res * 256;
		}
	}
}

void number_to_text(char* text, llint xc, unsigned & len) {
	//// Decode mpz_t xc into text so that the leftmost is the most significant
	if (xc < 0) {
		cerr << "\nError!: Message number negative!\n" << endl;
		exit(2);
	}

	llint r = 0;
	llint tmp = xc;
	char rev_text[64];
	unsigned i = 0;
	unsigned ch = r;
	//the following decodes from int to string, but in reverse order, also calculates the length of the outcome string.
	while (tmp > 0) {
		r = tmp % 256;
		ch = (unsigned)(r);
		rev_text[i] = ch;
		tmp = tmp / 256; // mpz_tdiv_q_ui(tmp, tmp, 256);
		i++;
	}
	rev_text[i] = 0; // finish up string
	len = i; // set len for stringlength for later use and this is returned in the function

	for (i = 0; i < len; i++) {
		text[i] = rev_text[len - i - 1];
	}
	text[len] = 0;// last char of string 
}







