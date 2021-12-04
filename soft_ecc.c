#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "memcached.h"

#ifdef SOFT_ECC
#define ITEM_META_SIZE (sizeof(item) + (settings.use_cas ? 8 : 0))
#define DATA_ARRAY_SIZE 4000
#define CORRECT_ERROR_NUM	2

struct item_ecc* make_item_ecc(item *it){
	int nkey = it->nkey + 1;
	int nbytes = it->nbytes;
	struct item_ecc * new;

	new = malloc(sizeof(struct item_ecc));
	new->m_size = ITEM_META_SIZE;
	new->k_size = nkey;
	new->v_size = nbytes;
	new->m_bch = new->k_bch = new->v_bch = NULL; //zeroing data array field

	//fprintf(stderr,"make_item_ecc: msize:%d,ksize:%d,vsize:%d\n",new->m_size,new->k_size,new->v_size);
	return new;
}

static void init_ecc_struct(struct item_ecc* it_ec){
	struct ecc_prms *m,*k,*v;

	m = malloc(sizeof(struct ecc_prms));
	k = malloc(sizeof(struct ecc_prms));
	v = malloc(sizeof(struct ecc_prms));
	it_ec->m_bch = m;
	it_ec->k_bch = k;
	it_ec->v_bch = v;
}

static int calc_prmtr(int d_size){
	int threshold = 32;
	int m;
	for(m = 5;m < 20;m++){
		if(d_size*8 < threshold - m*CORRECT_ERROR_NUM - 1){
			return m;
		}
		threshold *= 2;
	}
	return -1;
}

static void set_prmtr(int d_size, struct ecc_prms* ecc_p){
	ecc_p->m = calc_prmtr(d_size);
	ecc_p->t = CORRECT_ERROR_NUM;
	return;
}

static void calc_num_to_bin(int n, int* bin){
	//fprintf(stderr,"n:%d\n",n);
	for(int i = 7;i >= 0;i--){
		bin[i] = n % 2;
		n = n / 2;
	}
}

static int calc_bin_to_num(int*bin){
	int ret = 0;
	int mul = 1;
	for(int i = 7; i >= 0;i--){
		ret += bin[i]*mul;
		mul = mul * 2;
	}
	return ret;
}

static void generate_data_array(unsigned char* p, int size,int*bin){
	unsigned char* data;
	int b_index = 0;
	//fprintf(stderr,"printf data:\n");
	for(int i = 0;i < size; i++){
		data = &p[i];
		calc_num_to_bin((int)(*data),&bin[b_index]);
		//fprintf(stderr,"%d,",*data);
		//fprintf(stderr,"%d,%d%d%d%d%d%d%d%d\n",*data,bin[b_index + 0],bin[b_index +1],bin[b_index +2],bin[b_index +3],bin[b_index +4],bin[b_index +5],bin[b_index +6],bin[b_index +7]);
		b_index += 8;
	}
	//fprintf(stderr,"\n");
}


static void correct_data(unsigned char*p,int size,int *bin){
	int b_index = 0;
	for(int i = 0;i < size; i++){
		p[i] = calc_bin_to_num(&bin[b_index]);
		b_index += 8;
	}
}


static void set_zero(int* data, int size){
	memset(data,0,size*sizeof(int));
	// for(int i = 0;i < size; i++){
	// 	data[i] = 0;
	// }
}

struct item_ecc* set_ecc(item * it){
	struct item_ecc *it_ec = it->ecc_addr;
	int *m_data,*k_data,*v_data;

	init_ecc_struct(it_ec);
	// for meta data field
	set_prmtr(it_ec->m_size,it_ec->m_bch);
	m_data = malloc(sizeof(int)*DATA_ARRAY_SIZE);
	set_zero(m_data,DATA_ARRAY_SIZE);
	generate_data_array((unsigned char*)it,it_ec->m_size,m_data);

	//for key data field
	set_prmtr(it_ec->k_size,it_ec->k_bch);
	k_data = malloc(sizeof(int)*DATA_ARRAY_SIZE);
	set_zero(k_data,DATA_ARRAY_SIZE);
	generate_data_array((unsigned char*)ITEM_key(it),it_ec->k_size,k_data);

	//for valu data field
	set_prmtr(it_ec->v_size,it_ec->v_bch);
	v_data = malloc(sizeof(int)*DATA_ARRAY_SIZE);
	set_zero(v_data,DATA_ARRAY_SIZE);
	generate_data_array((unsigned char*)ITEM_data(it),it_ec->v_size,v_data);

	call_encode_bch(it_ec->m_bch,m_data);
	call_encode_bch(it_ec->k_bch,k_data);
	call_encode_bch(it_ec->v_bch,v_data);

	free(m_data);
	free(k_data);
	free(v_data);
	fprintf(stderr,"\tencode is done\n");
	return it_ec;
}

void check_ecc(item* it){
	if(it != NULL){
	struct item_ecc *it_ec = it->ecc_addr;
	int *m_data,*k_data,*v_data;
 
	// for meta data field
	m_data = malloc(sizeof(int)*DATA_ARRAY_SIZE);
	set_zero(m_data,DATA_ARRAY_SIZE);
	generate_data_array((unsigned char*)it,it_ec->m_size,m_data);

	//for key data field
	k_data = malloc(sizeof(int)*DATA_ARRAY_SIZE);
	set_zero(k_data,DATA_ARRAY_SIZE);
	generate_data_array((unsigned char*)ITEM_key(it),it_ec->k_size,k_data);


	//for valu data field
	v_data = malloc(sizeof(int)*DATA_ARRAY_SIZE);
	set_zero(v_data,DATA_ARRAY_SIZE);
	generate_data_array((unsigned char*)ITEM_data(it),it_ec->v_size,v_data);

	if(call_decode_bch(it_ec->m_bch,m_data) == 1){
		correct_data((unsigned char*)it,it_ec->m_size,m_data);
		fprintf(stderr,"meta data is corrected\n");
	}
	else{
		//fprintf(stderr,"there is no error\n");
	}

	// //fprintf(stderr,"before key:%s\n",ITEM_key(it));
	// if(call_decode_bch(it_ec->k_bch,k_data) == 1){
	// 	correct_data((unsigned char*)ITEM_key(it),it_ec->k_size,k_data);
	// }
	// //fprintf(stderr,"after key:%s\n",ITEM_key(it));

	// //fprintf(stderr,"before value:%s\n",ITEM_data(it));
	// if(call_decode_bch(it_ec->v_bch,v_data) == 1){
	// 	correct_data((unsigned char*)ITEM_data(it),it_ec->v_size,v_data);
	// }
	// //fprintf(stderr,"after value:%s\n",ITEM_data(it));
	
	free(m_data);
	free(k_data);
	free(v_data);
	}
	
}

void update_ecc(item* it,bool do_m,bool do_k, bool do_v){
	struct item_ecc *it_ec = it->ecc_addr;
	int *m_data,*k_data,*v_data;
	if(do_m){
		//encode metadata
		m_data = malloc(sizeof(int)*DATA_ARRAY_SIZE);
		generate_data_array((unsigned char*)it,it_ec->m_size,m_data);
		call_encode_bch(it_ec->m_bch,m_data);
		free(m_data);
	}
	if(do_k){
		//encode key data
		k_data = malloc(sizeof(int)*DATA_ARRAY_SIZE);
		generate_data_array((unsigned char*)ITEM_key(it),it_ec->k_size,k_data);
		call_encode_bch(it_ec->k_bch,k_data);
		free(k_data);
	}
	if(do_v){
		//encode value data
		v_data = malloc(sizeof(int)*DATA_ARRAY_SIZE);
		generate_data_array((unsigned char*)ITEM_data(it),it_ec->v_size,v_data);
		call_encode_bch(it_ec->v_bch,v_data);
		free(v_data);
	}

	return;
}

int incr_refcount(item* it){
	//fprintf(stderr,"\tin_bf:%p,%p,%p,%d,%d,%d,%d,%d,%d\n",(void*)it->next,(void*)it->prev,(void*)it->h_next,it->time,it->exptime,it->nbytes,it->refcount,it->it_flags,it->slabs_clsid);
	
	if(it && it->ecc_addr->m_bch){
		check_ecc(it);
	}

	int ret = ++(it->refcount);
	if(it->ecc_addr->m_bch != NULL){
		update_ecc(it,true,false,false);
		//fprintf(stderr,"\tin_af:%p,%p,%p,%d,%d,%d,%d,%d,%d\n",(void*)it->next,(void*)it->prev,(void*)it->h_next,it->time,it->exptime,it->nbytes,it->refcount,it->it_flags,it->slabs_clsid);
	}

	return ret;
}

int decr_refcount(item* it){
	//fprintf(stderr,"\tde_bf:%p,%p,%p,%d,%d,%d,%d,%d,%d\n",(void*)it->next,(void*)it->prev,(void*)it->h_next,it->time,it->exptime,it->nbytes,it->refcount,it->it_flags,it->slabs_clsid);
	if(it && it->ecc_addr->m_bch){
		check_ecc(it);
	}
	
	int ret = --(it->refcount);
	if(it->ecc_addr->m_bch != NULL){
		update_ecc(it,true,false,false);
		//fprintf(stderr,"\tde_af:%p,%p,%p,%d,%d,%d,%d,%d,%d\n",(void*)it->next,(void*)it->prev,(void*)it->h_next,it->time,it->exptime,it->nbytes,it->refcount,it->it_flags,it->slabs_clsid);
           
	}

	return ret;
}

void free_ecc(item *it){
	struct item_ecc* it_p;
	struct ecc_prms* pr_p;
	it_p = it->ecc_addr;

	//free meta data field
	pr_p  = it_p->m_bch;
	free(pr_p);

	//free key data field
	pr_p = it_p->k_bch;
	free(pr_p);

	//free value data field
	pr_p = it_p->v_bch;
	free(pr_p);

	//free item_ecc
	free(it_p);

	return;
}

#endif


// void soft_ecc(void){
//     fprintf(stderr,"hello\n");
// 	struct ecc_prms *new;
// 	new = malloc(sizeof(ecc_prms));
// 	new->m = 5;
// 	new->t = 2;
// 	new->length = 30;
// 	int data[1048576];
// 	data[0] = 0;
// 	data[1] = 1;
// 	data[2] = 0;
// 	data[3] = 1;
// 	data[4] = 0;
// 	data[5] = 1;
// 	data[6] = 0;
// 	data[7] = 1;
// 	data[8] = 0;
// 	data[9] = 1;

// 	call_encode_bch(new,data);
//     return;
// }


// #define MAX_LENGTH_OF_CODEWORD 64
// #define FAILUER_PROBABILITY 1e-10

// uint calcurate_k(uint s,uint l){
//     uint ret = (l + s - 1) / s;
//     return ret;
// }

// long double calcuratit_ec(uint s){
//     long double ret = FAILUER_PROBABILITY/(3 * s);
//     return ret;
// }

// void set_prmtr(uint length,uint* k,uint* n,uint* t){
//     uint s;
//     long double prob;
//     s = (length + MAX_LENGTH_OF_CODEWORD -1) / MAX_LENGTH_OF_CODEWORD;
//     prob = calcuratit_ec(s);
    
//     *k = calcurate_k(s,length);

// }







