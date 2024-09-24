// EC ElGamal.cpp : This file contains the 'main' function. Program execution begins and ends there.
//



#include <iostream>
#include "EC ElGamal.h"
using namespace std;

#define ui unsigned

int main()
{
    cout << "Hello World!\n" << endl;
	cout << "The following demonstrates an \nElliptic Curve ElGamal Cryptosystem.\n";
	
	llint a, b, p;
	llint mx, my; // point coordinates of the original message
	llint gx, gy, ord_gen; //parameters related to the generator
	
	cout << "\nGenerating parameters...\n\n";
	generate_parameters(a, b, p, gx, gy, ord_gen);
	cout << "Our prime p is between 2^19 and 2^20 with p = 3 (mod 4).\n";
	cout << "a and G (with small coordinates) are chosen at random,\n then b is chosen so that G is on the elliptic curve \ny^2=x^3+ax+b (mod p).\n";


	// Now we have E: y^2 = x^3+ax+b (mod p)
	// and a generator G=(gx,gy) with order = ord_gen
	// alpha will be made public later by Bob
	llint alpha_x = gx;
	llint alpha_y = gy;
	// beta will be made public later by Bob
	llint beta_x, beta_y;

	//Bob's secret is the number prk_B
	srand(time(NULL) + rand()%1907 );
	llint prk_B = rand() % p;
	
	// Bob calculates 
	k_times(prk_B, alpha_x, alpha_y, beta_x, beta_y, a, b, p);

	printf("\nBob's private key=prk_B=%lld\n", prk_B);
	printf("Bob publishes \nalpha=(%lld,%lld) \nbeta=prk_B*alpha=(%lld,%lld)\n", alpha_x, alpha_y, beta_x, beta_y);
	
	//Bob makes p,a, n, alpha = (alpha_x, alpha_y) , beta = (beta_x, beta_y) public

	unsigned prec = 4; // the precision in encoding/decoding: Has to be the same in both encode and decode functions

	unsigned max_mess_size = log2(p / prec);
	cout << "\nWe can safely encode/decode numbers up to \n"<< pow(2, max_mess_size) << "=2^" << max_mess_size << "-bits=" << max_mess_size / 8 << " bytes/characters at once." << endl;
	//llint m = rand() % (llint) pow(2, max_mess_size); //message of Alice (randomized for now)

	char message_text[] = "AB";
	llint m;
	unsigned len = strlen(message_text);
	text_to_number(m, message_text, len);

	//Alice's secret is the number prk_A
	srand(time(NULL) + rand() % m);
	llint prk_A = rand() % p;

	// m is represented as (mx,my) on the elliptic curve
	encode(m, mx, my, a, b, p, prec);
	printf("\nAlice's message \"%s\" was translated to the number %lld \nand encoded as P=(%lld,%lld).\n",message_text, m, mx, my);

	
	// Encryption
	
	cout << "\nEncrypting...\n";
	// Alice calculates
	// y1 = (y1_x , y1_y) will be prk_A * alpha
	// y2 = (y2_x , y2_y) will be m + prk_A * beta
	llint y1_x, y1_y, y2_x, y2_y;
	
	k_times(prk_A, alpha_x, alpha_y, y1_x, y1_y, a, b, p);
	k_times(prk_A,  beta_x,  beta_y, y2_x, y2_y, a, b, p);
	add(mx, my, y2_x, y2_y, y2_x, y2_y, a, b, p);


	printf("Alice's private key=prk_A=%lld\n", prk_A);
	printf("Alice publishes \ny1=prk_A*alpha=(%lld, %lld) \ny2=P+prk_A*beta=(%lld, %lld)\n", y1_x, y1_y, y2_x, y2_y);

	
	// Decryption

	cout << "\nDecrypting...\n";
	// Bob receives y1 and y2
	// He calculates the message from Alice as m = y2 - prk_B *y1

	llint rm, rx, ry; // recovered message and its coordinates

	k_times(prk_B, y1_x, y1_y,rx,ry ,a, b, p);
	add(y2_x, y2_y, rx, -ry, rx,ry ,a, b, p);
	// Bob has now calculated the message from Alice as (rx,ry)

	
	decode(rm, rx, ry, a, b, p, prec);
	printf("Bob decrypts the message as \ny2-prk_B*y1=(%lld,%lld)\nand decodes it to receive Alice's original message as the number %lld.\n", rx, ry,rm);

	char recovered_text[3] ;
	number_to_text(recovered_text, rm, len);
	printf("Bob translates %lld to recover the text \"%s\"\n\n",rm, recovered_text);

	cout << "Have a nice day!" << endl;
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

