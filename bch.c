
/*
 * File:    bch3.c
 * Title:   Encoder/decoder for binary BCH codes in C (Version 3.1)
 * Author:  Robert Morelos-Zaragoza
 * Date:    August 1994
 * Revised: June 13, 1997
 *
 * ===============  Encoder/Decoder for binary BCH codes in C =================
 *
 * Version 1:   Original program. The user provides the generator polynomial
 *              of the code (cumbersome!).
 * Version 2:   Computes the generator polynomial of the code.
 * Version 3:   No need to input the coefficients of a primitive polynomial of
 *              degree m, used to construct the Galois Field GF(2**m). The
 *              program now works for any binary BCH code of length such that:
 *              2**(m-1) - 1 < length <= 2**m - 1
 *
 * Note:        You may have to change the size of the arrays to make it work.
 *
 * The encoding and decoding methods used in this program are based on the
 * book "Error Control Coding: Fundamentals and Applications", by Lin and
 * Costello, Prentice Hall, 1983.
 *
 * Thanks to Patrick Boyle (pboyle@era.com) for his observation that 'bch2.c'
 * did not work for lengths other than 2**m-1 which led to this new version.
 * Portions of this program are from 'rs.c', a Reed-Solomon encoder/decoder
 * in C, written by Simon Rockliff (simon@augean.ua.oz.au) on 21/9/89. The
 * previous version of the BCH encoder/decoder in C, 'bch2.c', was written by
 * Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu) on 5/19/92.
 *
 * NOTE:    
 *          The author is not responsible for any malfunctioning of
 *          this program, nor for any damage caused by it. Please include the
 *          original program along with these comments in any redistribution.
 *
 *  For more information, suggestions, or other ideas on implementing error
 *  correcting codes, please contact me at:
 *
 *                           Robert Morelos-Zaragoza
 *                           5120 Woodway, Suite 7036
 *                           Houston, Texas 77056
 *
 *                    email: r.morelos-zaragoza@ieee.org
 *
 * COPYRIGHT NOTICE: This computer program is free for non-commercial purposes.
 * You may implement this program for any non-commercial application. You may 
 * also implement this program for commercial purposes, provided that you
 * obtain my written permission. Any modification of this program is covered
 * by this copyright.
 *
 * == Copyright (c) 1994-7,  Robert Morelos-Zaragoza. All rights reserved.  ==
 *
 * m = order of the Galois field GF(2**m) 
 * n = 2**m - 1 = size of the multiplicative group of GF(2**m)
 * length = length of the BCH code
 * t = error correcting capability (max. no. of errors the code corrects)
 * d = 2*t + 1 = designed min. distance = no. of consecutive roots of g(x) + 1
 * k = n - deg(g(x)) = dimension (no. of information bits/codeword) of the code
 * p[] = coefficients of a primitive polynomial used to generate GF(2**m)
 * g[] = coefficients of the generator polynomial, g(x)
 * alpha_to [] = log table of GF(2**m) 
 * index_of[] = antilog table of GF(2**m)
 * data[] = information bits = coefficients of data polynomial, i(x)
 * bb[] = coefficients of redundancy polynomial x^(length-k) i(x) modulo g(x)
 * numerr = number of errors 
 * errpos[] = error positions 
 * recd[] = coefficients of the received polynomial 
 * decerror = number of decoding errors (in _message_ positions) 
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bch.h"

#define RECD_SIZE		1048576
#define ALPHA_TO_SIZE	1048576
#define INDEX_OF_SIZE	1048576
#define G_SIZE			548576

//int             p[21];
//int             alpha_to[1048576], index_of[1048576], g[548576];
int             recd[1048576]; //data[1048576];//, bb[548576];
//int             seed;
//int             numerr, errpos[1024], decerror = 0;



static void read_p(struct ecc_prms* ecc_p)
/*
 *	Read m, the degree of a primitive polynomial p(x) used to compute the
 *	Galois field GF(2**m). Get precomputed coefficients p[] of p(x). Read
 *	the code length.
 */
{
    int			i;

	//nezu add----------------------------------------------------------------------
	int		m = ecc_p->m;
	int 	*p = ecc_p->p;
	int n;
	//---------------------------------------------------------------------------------
	
	for (i=1; i<m; i++)
		p[i] = 0;
	p[0] = p[m] = 1;
	if (m == 2)			p[1] = 1;
	else if (m == 3)	p[1] = 1;
	else if (m == 4)	p[1] = 1;
	else if (m == 5)	p[2] = 1;
	else if (m == 6)	p[1] = 1;
	else if (m == 7)	p[1] = 1;
	else if (m == 8)	p[4] = p[5] = p[6] = 1;
	else if (m == 9)	p[4] = 1;
	else if (m == 10)	p[3] = 1;
	else if (m == 11)	p[2] = 1;
	else if (m == 12)	p[3] = p[4] = p[7] = 1;
	else if (m == 13)	p[1] = p[3] = p[4] = 1;
	else if (m == 14)	p[1] = p[11] = p[12] = 1;
	else if (m == 15)	p[1] = 1;
	else if (m == 16)	p[2] = p[3] = p[5] = 1;
	else if (m == 17)	p[3] = 1;
	else if (m == 18)	p[7] = 1;
	else if (m == 19)	p[1] = p[5] = p[6] = 1;
	else if (m == 20)	p[3] = 1;
	//printf("p(x) = ");
    n = 1;
	for (i = 0; i <= m; i++) {
        n *= 2;
		//printf("%1d", p[i]);
    }
	//printf("\n");
	n = n / 2 - 1;

	//nezu add--------------------------------------------------------------------------
	ecc_p->n = n;
	ecc_p->length = n;
	//---------------------------------------------------------------------------------------
}


static void generate_gf(struct ecc_prms* ecc_p,int* alpha_to,int* index_of)
/*
 * Generate field GF(2**m) from the irreducible polynomial p(X) with
 * coefficients in p[0]..p[m].
 *
 * Lookup tables:
 *   index->polynomial form: alpha_to[] contains j=alpha^i;
 *   polynomial form -> index form:	index_of[j=alpha^i] = i
 *
 * alpha=2 is the primitive element of GF(2**m) 
 */
{
	register int    i, mask;

	//nezu add--------------------------------------------------------------------------
	int m = ecc_p->m;
	int n = ecc_p->n;
	int *p = ecc_p->p;
	//-----------------------------------------------------------------------------------

	mask = 1;
	alpha_to[m] = 0;
	for (i = 0; i < m; i++) {
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;
		if (p[i] != 0)
			alpha_to[m] ^= mask;
		mask <<= 1;
	}
	index_of[alpha_to[m]] = m;
	mask >>= 1;
	for (i = m + 1; i < n; i++) {
		if (alpha_to[i - 1] >= mask)
		  alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
		else
		  alpha_to[i] = alpha_to[i - 1] << 1;
		index_of[alpha_to[i]] = i;
	}
	index_of[0] = -1;
}

static void gen_poly(struct ecc_prms* ecc_p,int*alpha_to,int*index_of,int*g)
/*
 * Compute the generator polynomial of a binary BCH code. Fist generate the
 * cycle sets modulo 2**m - 1, cycle[][] =  (i, 2*i, 4*i, ..., 2^l*i). Then
 * determine those cycle sets that contain integers in the set of (d-1)
 * consecutive integers {1..(d-1)}. The generator polynomial is calculated
 * as the product of linear factors of the form (x+alpha^i), for every i in
 * the above cycle sets.
 */
{
	register int	ii, jj, ll, kaux;
	register int	test, aux, nocycles, root, noterms, rdncy;
	int             cycle[1024][21], size[1024], min[1024], zeros[1024];

	//int m = ecc_p->m;
	int n = ecc_p->n;
	int t = ecc_p->t;
	int d,k;
	int length = ecc_p->length;

	/* Generate cycle sets modulo n, n = 2**m - 1 */
	cycle[0][0] = 0;
	size[0] = 1;
	cycle[1][0] = 1;
	size[1] = 1;
	jj = 1;			/* cycle set index */
	// if (m > 9)  {
	// 	printf("Computing cycle sets modulo %d\n", n);
	// 	printf("(This may take some time)...\n");
	// }
	do {
		/* Generate the jj-th cycle set */
		ii = 0;
		do {
			ii++;
			cycle[jj][ii] = (cycle[jj][ii - 1] * 2) % n;
			size[jj]++;
			aux = (cycle[jj][ii] * 2) % n;
		} while (aux != cycle[jj][0]);
		/* Next cycle set representative */
		ll = 0;
		do {
			ll++;
			test = 0;
			for (ii = 1; ((ii <= jj) && (!test)); ii++)	
			/* Examine previous cycle sets */
			  for (kaux = 0; ((kaux < size[ii]) && (!test)); kaux++)
			     if (ll == cycle[ii][kaux])
			        test = 1;
		} while ((test) && (ll < (n - 1)));
		if (!(test)) {
			jj++;	/* next cycle set index */
			cycle[jj][0] = ll;
			size[jj] = 1;
		}
	} while (ll < (n - 1));
	nocycles = jj;		/* number of cycle sets modulo n */

	// printf("Enter the error correcting capability, t: ");
	// scanf("%d", &t);

	d = 2 * t + 1;
	ecc_p->d = d;

	/* Search for roots 1, 2, ..., d-1 in cycle sets */
	kaux = 0;
	rdncy = 0;
	for (ii = 1; ii <= nocycles; ii++) {
		min[kaux] = 0;
		test = 0;
		for (jj = 0; ((jj < size[ii]) && (!test)); jj++)
			for (root = 1; ((root < d) && (!test)); root++)
				if (root == cycle[ii][jj])  {
					test = 1;
					min[kaux] = ii;
				}
		if (min[kaux]) {
			rdncy += size[min[kaux]];
			kaux++;
		}
	}
	noterms = kaux;
	kaux = 1;
	for (ii = 0; ii < noterms; ii++)
		for (jj = 0; jj < size[min[ii]]; jj++) {
			zeros[kaux] = cycle[min[ii]][jj];
			kaux++;
		}

	k = length - rdncy;

    if (k<0)
      {
         printf("\tParameters invalid!\n");
         exit(0);
      }
	ecc_p->k = k;

	//printf("This is a (%d, %d, %d) binary BCH code\n", length, k, d);

	/* Compute the generator polynomial */
	g[0] = alpha_to[zeros[1]];
	g[1] = 1;		/* g(x) = (X + zeros[1]) initially */
	for (ii = 2; ii <= rdncy; ii++) {
	  g[ii] = 1;
	  for (jj = ii - 1; jj > 0; jj--)
	    if (g[jj] != 0)
	      g[jj] = g[jj - 1] ^ alpha_to[(index_of[g[jj]] + zeros[ii]) % n];
	    else
	      g[jj] = g[jj - 1];
	  g[0] = alpha_to[(index_of[g[0]] + zeros[ii]) % n];
	}
	// printf("Generator polynomial:\ng(x) = ");
	// for (ii = 0; ii <= rdncy; ii++) {
	//   printf("%d", g[ii]);
	//   if (ii && ((ii % 50) == 0))
	//     printf("\n");
	// }
	// printf("\n");
}


static void encode_bch(struct ecc_prms* ecc_p,int *data,int *bb,int*g)
/*
 * Compute redundacy bb[], the coefficients of b(x). The redundancy
 * polynomial b(x) is the remainder after dividing x^(length-k)*data(x)
 * by the generator polynomial g(x).
 */
{
	register int    i, j;
	register int    feedback;
	int length = ecc_p->length;
	int k = ecc_p->k;

	for (i = 0; i < length - k; i++)
		bb[i] = 0;
	for (i = k - 1; i >= 0; i--) {
		feedback = data[i] ^ bb[length - k - 1];
		if (feedback != 0) {
			for (j = length - k - 1; j > 0; j--)
				if (g[j] != 0)
					bb[j] = bb[j - 1] ^ feedback;
				else
					bb[j] = bb[j - 1];
			bb[0] = g[0] && feedback;
		} else {
			for (j = length - k - 1; j > 0; j--)
				bb[j] = bb[j - 1];
			bb[0] = 0;
		}
	}

	// nezu add for debug----------------------------------------------------------
	// fprintf(stderr,"original data  = ");
	// for (int t = 0; t < k; t++) {
	// 	fprintf(stderr,"%1d", data[t]);
	// }
	// fprintf(stderr,"\n");
	// fprintf(stderr,"encode\n");
	// for(int t = 0; t < length - k ;t++){
	// 	fprintf(stderr,"%d",bb[t]);
	// }
	// fprintf(stderr,"\n");
}


static int decode_bch(struct ecc_prms* ecc_p, int* recd,int *alpha_to,int *index_of)
/*
 * Simon Rockliff's implementation of Berlekamp's algorithm.
 *
 * Assume we have received bits in recd[i], i=0..(n-1).
 *
 * Compute the 2*t syndromes by substituting alpha^i into rec(X) and
 * evaluating, storing the syndromes in s[i], i=1..2t (leave s[0] zero) .
 * Then we use the Berlekamp algorithm to find the error location polynomial
 * elp[i].
 *
 * If the degree of the elp is >t, then we cannot correct all the errors, and
 * we have detected an uncorrectable error pattern. We output the information
 * bits uncorrected.
 *
 * If the degree of elp is <=t, we substitute alpha^i , i=1..n into the elp
 * to get the roots, hence the inverse roots, the error location numbers.
 * This step is usually called "Chien's search".
 *
 * If the number of errors located is not equal the degree of the elp, then
 * the decoder assumes that there are more than t errors and cannot correct
 * them, only detect them. We output the information bits uncorrected.
 */
{
	register int    i, j, u, q, t2, count = 0, syn_error = 0;
	int             elp[1026][1024], d[1026], l[1026], u_lu[1026], s[1025];
	//int             root[200], loc[200], err[1024], reg[201];
	int				loc[200], reg[201];

	//nezu add---------------------------------------------------------------
	int t = ecc_p->t;
	int n = ecc_p->n;
	int length = ecc_p->length;

	t2 = 2 * t;

	/* first form the syndromes */
	//printf("S(x) = ");
	for (i = 1; i <= t2; i++) {
		s[i] = 0;
		for (j = 0; j < length; j++)
			if (recd[j] != 0)
				s[i] ^= alpha_to[(i * j) % n];
		if (s[i] != 0)
			syn_error = 1; /* set error flag if non-zero syndrome */
/*
 * Note:    If the code is used only for ERROR DETECTION, then
 *          exit program here indicating the presence of errors.
 */
		/* convert syndrome from polynomial form to index form  */
		s[i] = index_of[s[i]];
		//printf("%3d ", s[i]);
	}
	//printf("\n");

	if (syn_error) {	/* if there are errors, try to correct them */
		/*
		 * Compute the error location polynomial via the Berlekamp
		 * iterative algorithm. Following the terminology of Lin and
		 * Costello's book :   d[u] is the 'mu'th discrepancy, where
		 * u='mu'+1 and 'mu' (the Greek letter!) is the step number
		 * ranging from -1 to 2*t (see L&C),  l[u] is the degree of
		 * the elp at that step, and u_l[u] is the difference between
		 * the step number and the degree of the elp. 
		 */
		/* initialise table entries */
		d[0] = 0;			/* index form */
		d[1] = s[1];		/* index form */
		elp[0][0] = 0;		/* index form */
		elp[1][0] = 1;		/* polynomial form */
		for (i = 1; i < t2; i++) {
			elp[0][i] = -1;	/* index form */
			elp[1][i] = 0;	/* polynomial form */
		}
		l[0] = 0;
		l[1] = 0;
		u_lu[0] = -1;
		u_lu[1] = 0;
		u = 0;
 
		do {
			u++;
			if (d[u] == -1) {
				l[u + 1] = l[u];
				for (i = 0; i <= l[u]; i++) {
					elp[u + 1][i] = elp[u][i];
					elp[u][i] = index_of[elp[u][i]];
				}
			} else
				/*
				 * search for words with greatest u_lu[q] for
				 * which d[q]!=0 
				 */
			{
				q = u - 1;
				while ((d[q] == -1) && (q > 0))
					q--;
				/* have found first non-zero d[q]  */
				if (q > 0) {
				  j = q;
				  do {
				    j--;
				    if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
				      q = j;
				  } while (j > 0);
				}
 
				/*
				 * have now found q such that d[u]!=0 and
				 * u_lu[q] is maximum 
				 */
				/* store degree of new elp polynomial */
				if (l[u] > l[q] + u - q)
					l[u + 1] = l[u];
				else
					l[u + 1] = l[q] + u - q;
 
				/* form new elp(x) */
				for (i = 0; i < t2; i++)
					elp[u + 1][i] = 0;
				for (i = 0; i <= l[q]; i++)
					if (elp[q][i] != -1)
						elp[u + 1][i + u - q] = 
                                   alpha_to[(d[u] + n - d[q] + elp[q][i]) % n];
				for (i = 0; i <= l[u]; i++) {
					elp[u + 1][i] ^= elp[u][i];
					elp[u][i] = index_of[elp[u][i]];
				}
			}
			u_lu[u + 1] = u - l[u + 1];
 
			/* form (u+1)th discrepancy */
			if (u < t2) {	
			/* no discrepancy computed on last iteration */
			  	if (s[u + 1] != -1){
			    	d[u + 1] = alpha_to[s[u + 1]];
				  }
			  	else{
			    	d[u + 1] = 0;
				}
				for (i = 1; i <= l[u + 1]; i++){
			    	if ((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0)){
			        d[u + 1] ^= alpha_to[(s[u + 1 - i] 
			                      + index_of[elp[u + 1][i]]) % n];
					}
				}
			  /* put d[u+1] into index form */
				d[u + 1] = index_of[d[u + 1]];	
			}
		} while ((u < t2) && (l[u + 1] <= t));
 
		u++;
		if (l[u] <= t) {/* Can correct errors */
			/* put elp into index form */
			for (i = 0; i <= l[u]; i++)
				elp[u][i] = index_of[elp[u][i]];

			// printf("sigma(x) = ");
			// for (i = 0; i <= l[u]; i++)
			// 	printf("%3d ", elp[u][i]);
			// printf("\n");
			// printf("Roots: ");

			/* Chien search: find roots of the error location polynomial */
			for (i = 1; i <= l[u]; i++)
				reg[i] = elp[u][i];
			count = 0;
			for (i = 1; i <= n; i++) {
				q = 1;
				for (j = 1; j <= l[u]; j++)
					if (reg[j] != -1) {
						reg[j] = (reg[j] + j) % n;
						q ^= alpha_to[reg[j]];
					}
				if (!q) {	/* store root and error
						 * location number indices */
					//root[count] = i;  //nezu erased
					loc[count] = n - i; //nezu erased
					count++;
					// printf("%3d ", n - i);
				}
			}
			//printf("\n");
			if (count == l[u])	
			/* no. roots = degree of elp hence <= t errors */
				for (i = 0; i < l[u]; i++)
					recd[loc[i]] ^= 1;
			else{	/* elp has degree >t hence cannot solve */
				printf("\tIncomplete decoding: errors detected\n");
				return 2;
			}
			return 1;
		}
		return 2;
	}
	return 0;
}


void call_encode_bch(struct ecc_prms* ecc_p, int* data){

	if(ecc_p == NULL){
		fprintf(stderr,"ecc_p is NULL\n");
		return;
	}

	int *alpha_to,*index_of,*g;
	alpha_to = malloc(ALPHA_TO_SIZE);
	index_of = malloc(INDEX_OF_SIZE);
	g 		 = malloc(G_SIZE);

	read_p(ecc_p);
	generate_gf(ecc_p,alpha_to,index_of);          /* Construct the Galois Field GF(2**m) */
	gen_poly(ecc_p,alpha_to,index_of,g);             /* Compute the generator polynomial of BCH code */

	encode_bch(ecc_p,data,ecc_p->bb,g);           /* encode data */

	// fprintf(stderr,"bbbbbb\n");
	// for(int s=0;s<ecc_p->length-ecc_p->k;s++){
	// 	fprintf(stderr,"%d",ecc_p->bb[s]);
	// }
	// fprintf(stderr,"\n");
	free(alpha_to);
	free(index_of);
	free(g);
}

static int check_errors(int length,int k,int*bb,int*cmp_bb){
	for (int i = 0; i < length -k ;i ++){
		//fprintf(stderr,"bb:%d,cbb:%d\n",bb[i],cmp_bb[i]);
		if(bb[i] != cmp_bb[i]){
			return 1;
		}
	}
	return 0;
}

int call_decode_bch(struct ecc_prms*ecc_p,int *data){
	int cmp_bb[BB_SIZE];
	int length = ecc_p->length;
	int *bb = ecc_p->bb;
	int k = ecc_p->k;
	int *alpha_to,*index_of,*g,*recd;
	int is_error = 0;
	int ret = 0;
	alpha_to = malloc(ALPHA_TO_SIZE);
	index_of = malloc(INDEX_OF_SIZE);
	g 		 = malloc(G_SIZE);

	generate_gf(ecc_p,alpha_to,index_of);          /* Construct the Galois Field GF(2**m) */
	gen_poly(ecc_p,alpha_to,index_of,g);             /* Compute the generator polynomial of BCH code */

	encode_bch(ecc_p,data,cmp_bb,g);

	if((is_error = check_errors(length,k,bb,cmp_bb))){
		fprintf(stderr,"\tThere is some error and try to correct error.\n");
		//prepare recd array
		recd = malloc(RECD_SIZE);
		for (int i = 0; i < length - k; i++)
				recd[i] = bb[i];
		for (int i = 0; i < k; i++)
			recd[i + length - k] = data[i];

		ret = decode_bch(ecc_p,recd,alpha_to,index_of);
		if(ret == 1 || ret == 2){
			if(ret == 1) // error is detected and it is corrected
			{
				for (int i = 0; i < k; i++) {
					data[i] = recd[i + length -k];
				}
				fprintf(stderr,"\t\terror is modified\n");
			}
			else{
				fprintf(stderr,"\t\tcouldn't correct error.\n");
			}
		}
		free(recd);
	}
	else{
		//fprintf(stderr,"\tThere is no error\n");
	}

	free(alpha_to);
	free(index_of);
	free(g);
	return ret;
}