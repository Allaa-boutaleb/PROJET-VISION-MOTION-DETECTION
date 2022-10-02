/* -------------------- */
/* --- sigmadelta.h --- */
/* -------------------- */

/*
 * Copyright (c) 2004 - 2021, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique
 * Sorbonne University - Laboratoire d'Informatique de Paris 6
 */

#ifndef __SIGMADELTA_H__
#define __SIGMADELTA_H__

#ifdef VERBOSE_PRAGMA
#pragma message ("- *** include sigmadelta.h ***")
#endif

#ifdef __cplusplus
#ifdef PRAGMA_VERBOSE
#pragma message ("C++")
#endif
extern "C" {
#endif

#define SD_VMAX 253
#define SD_VMIN 2


// Macro pour renvoyer le max ou le min entre deux entiers
#define MAX(x,y) ( ((x) > (y)) ? (x) : (y)  )  
#define MIN(x,y) ( ((x) < (y)) ? (x) : (y)  )  


void SigmaDelta_Step0_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E,        int j0, int j1);
void SigmaDelta_1Step_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int k, int j0, int j1);

void SigmaDelta_Step0(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E,        int i0, int i1, int j0, int j1);
void SigmaDelta_1Step(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E, int k, int i0, int i1, int j0, int j1);

#ifdef __cplusplus
}
#endif

#endif // __SIGMADELTA_H__
