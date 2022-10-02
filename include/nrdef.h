/* --------------- */
/* --- nrdef.h --- */
/* --------------- */

/*
 * Copyright (c) 2000 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique 
 */

#ifndef __NRDEF_H__
#define __NRDEF_H__

#define TRUE 1
#define FALSE 0

#define IMAGE_EXPORT(X) X

#define load1(X,i) X[i]
#define load2(X,i,j) X[i][j]
#define store1(X,i,x) X[i] = x
#define store2(X,i,j,x) X[i][j] = x


#define min2(e1,e2) e1&e2
#define min3(e1,e2,e3) e1&e2&e3

#define max2(e1,e2) e1|e2
#define max3(e1,e2,e3) e1|e2|e3



#endif // __NRDEF_H__
