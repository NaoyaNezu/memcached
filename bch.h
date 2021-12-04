/*
 * m = order of the Galois field GF(2**m) 
 * n = 2**m - 1 = size of the multiplicative group of GF(2**m)
 * length = length of the BCH code
 * t = error correcting capability (max. no. of errors the code corrects)
 * d = 2*t + 1 = designed min. distance = no. of consecutive roots of g(x) + 1
 * k = n - deg(g(x)) = dimension (no. of information bits/codeword) of the code
*/
#define BB_SIZE 548576

struct ecc_prms{
	int m;
	int n;
	int length;
	int k;
	int t;
	int d;
	int bb[BB_SIZE];
	int p[21];
};

// should not be refferd by external file.
// void generate_gf(struct ecc_prms* prmtr);
// void gen_poly(struct ecc_prms* prmtr);
// void encode_bch(struct ecc_prms* prmtr,int *data,int *bb);

//for external call
void call_encode_bch(struct ecc_prms* new, int* data);
int call_decode_bch(struct ecc_prms*prmtr,int *data);
//void read_p(struct ecc_prms* prmtr);
