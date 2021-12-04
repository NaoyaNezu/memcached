#include "memcached.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef FI
void print_address(char *obj_name,unsigned long head,unsigned long tail,size_t size,int elements){
    fprintf(stderr,"%s : \n\tAddress : 0x%lx ~ 0x%lx , element size: %3ld Byte , number of elements : %5d\n",obj_name,head,tail,size,elements);
}

void fault_inject(void* p, int offset){
    char *addr = p;
    *addr ^= (1 << (offset-1));
    return;
}

#endif