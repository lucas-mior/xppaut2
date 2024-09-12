#ifndef _storage_h_
#define _storage_h_
#include "integers.h"

void init_alloc_info(void);
void alloc_meth(void);
void init_stor(int32 nrow, int32 ncol);
void free_storage(int32 ncol);
int32 reallocstor(int32 ncol, int32 nrow);

#endif
