/* -------------------- */
/* --- sigmadelta.c --- */
/* -------------------- */

/*
* Copyright (c) 2004 - 2021, Lionel Lacassagne, All rights reserved
* University of Paris Sud, Laboratoire de Recherche en Informatique
* Sorbonne University - Laboratoire d'Informatique de Paris 6
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"

#include "sigmadelta.h"

// -----------------------------------------------------------------------------------------
void SigmaDelta_Step0_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int j0, int j1)
// -----------------------------------------------------------------------------------------
{
    // Initialisation
    for (int j = j0 ; j < j1 ; j++) {
        M[j] = I[j];
        V[j] = SD_VMIN;
    }
 
}
// ------------------------------------------------------------------------------------------------
void SigmaDelta_1Step_line(uint8 *I, uint8 *M, uint8 *O, uint8 *V, uint8 *E, int k, int j0, int j1)
// ------------------------------------------------------------------------------------------------
{

    // Mt estimation
    for (int j = j0 ; j < j1 ; j++) {
        if (M[j] < I[j]) M[j] += 1;
        else if (M[j] > I[j]) M[j] -= 1;
    }  

    // Ot computation
    for (int j = j0 ; j < j1 ; j++) {
        O[j] = abs(M[j] - I[j]);
    }

    // Vt updating and clamping
    for (int j = j0 ; j < j1 ; j++) {
        if (V[j] < k * O[j]) V[j] += 1;
        else if (V[j] > k * O[j]) V[j] -= 1;
        else V[j] = MAX(MIN(V[j], SD_VMAX),SD_VMIN);
    }

    // Et estimation
    for (int j = j0 ; j < j1 ; j++) {
        if (O[j] < V[j]) E[j] = 0;
        else E[j] = 1;
    }

}
// ---------------------------------------------------------------------------------------------------------
void SigmaDelta_Step0(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E, int i0, int i1, int j0, int j1)
// ---------------------------------------------------------------------------------------------------------
{
    for (int i = i0 ; i < i1 ; i++) {
        SigmaDelta_Step0_line(I[i], M[i], O[i], V[i], E[i], j0, j1);
    }
}
// ----------------------------------------------------------------------------------------------------------------
void SigmaDelta_1Step(uint8 **I, uint8 **M, uint8 **O, uint8 **V, uint8 **E, int k, int i0, int i1, int j0, int j1)
// ----------------------------------------------------------------------------------------------------------------
{
    for (int i = i0 ; i < i1 ; i++) {
        SigmaDelta_1Step_line(I[i], M[i], O[i], V[i], E[i], k, j0, j1);
    }
}
