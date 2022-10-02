/* -------------------------- */
/* --- morpho_ouverture.h --- */
/* -------------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique 
 */

#ifndef __MORPHO_OUVERTURE_H__
#define __MORPHO_OUVERTURE_H__

#ifdef __cplusplus
#ifdef PRAGMA_VERBOSE
#pragma message ("C++")
#endif
extern "C" {
#endif

// =====================
// == fonction lignes ==
// =====================

// fusion d'operateurs
void line_ouverture3_ui8matrix_fusion                     (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_ouverture3_ui8matrix_fusion_ilu5_red            (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_ouverture3_ui8matrix_fusion_ilu5_elu2_red       (uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_ouverture3_ui8matrix_fusion_ilu5_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y);
void line_ouverture3_ui8matrix_fusion_ilu15_red           (uint8 **X, int i, int j0, int j1, uint8 **Y);

// pipeline d'operateurs
void line_ouverture3_ui8matrix_pipeline_bassic              (uint8 **X, int i, int j0, int j1, uint8 **T, uint8 **Y);
void line_ouverture3_ui8matrix_pipeline_red                 (uint8 **X, int i, int j0, int j1, uint8 **T, uint8 **Y);
void line_ouverture3_ui8matrix_pipeline_ilu3_red            (uint8 **X, int i, int j0, int j1, uint8 **T, uint8 **Y);
void line_ouverture3_ui8matrix_pipeline_elu2_red            (uint8 **X, int i, int j0, int j1, uint8 **T, uint8 **Y);
void line_ouverture3_ui8matrix_pipeline_elu2_red_factor     (uint8 **X, int i, int j0, int j1, uint8 **T, uint8 **Y);
void line_ouverture3_ui8matrix_pipeline_ilu3_elu2_red       (uint8 **X, int i, int j0, int j1, uint8 **T, uint8 **Y);
void line_ouverture3_ui8matrix_pipeline_ilu3_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **T, uint8 **Y);

// ========================
// == fonction completes ==
// ========================

// composition d'operateurs
void ouverture3_ui8matrix_basic                      (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y, uint8 **Z);

// fusion d'operateurs
void ouverture3_ui8matrix_fusion                     (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void ouverture3_ui8matrix_fusion_ilu5_red            (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void ouverture3_ui8matrix_fusion_ilu5_elu2_red       (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void ouverture3_ui8matrix_fusion_ilu5_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);
void ouverture3_ui8matrix_fusion_ilu15_red           (uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y);

// pipeline d'operateurs
void ouverture3_ui8matrix_pipeline_basic               (uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y);
void ouverture3_ui8matrix_pipeline_red                 (uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y);
void ouverture3_ui8matrix_pipeline_ilu3_red            (uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y);
void ouverture3_ui8matrix_pipeline_elu2_red            (uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y);
void ouverture3_ui8matrix_pipeline_elu2_red_factor     (uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y);
void ouverture3_ui8matrix_pipeline_ilu3_elu2_red       (uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y);
void ouverture3_ui8matrix_pipeline_ilu3_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y);

#ifdef __cplusplus
}
#endif



#endif // __MORPHO_OUVERTURE_H__
