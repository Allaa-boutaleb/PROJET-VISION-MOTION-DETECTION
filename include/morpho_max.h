/* -------------------- */
/* --- morpho_max.h --- */
/* -------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique 
 */

#ifndef __MORPHO_MAX_H__
#define __MORPHO_MAX_H__

#ifdef __cplusplus
#ifdef PRAGMA_VERBOSE
#pragma message ("C++")
#endif
extern "C" {
#endif

// line functions
void line_max3_ui8matrix_basic               (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_max3_ui8matrix_reg                 (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_max3_ui8matrix_rot                 (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_max3_ui8matrix_red                 (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_max3_ui8matrix_ilu3                (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_max3_ui8matrix_ilu3_red            (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_max3_ui8matrix_elu2_red            (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_max3_ui8matrix_elu2_red_factor     (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_max3_ui8matrix_ilu3_elu2_red       (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_max3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y);
    
// full-matrix functions
void max3_ui8matrix_basic               (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void max3_ui8matrix_reg                 (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void max3_ui8matrix_rot                 (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void max3_ui8matrix_red                 (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void max3_ui8matrix_ilu3                (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void max3_ui8matrix_ilu3_red            (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void max3_ui8matrix_elu2_red            (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void max3_ui8matrix_elu2_red_factor     (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void max3_ui8matrix_ilu3_elu2_red       (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void max3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
    
#ifdef __cplusplus
}
#endif

#endif // __MORPHO_MAX_H__
