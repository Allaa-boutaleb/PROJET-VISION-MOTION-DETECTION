/* --------------------- */
/* ---- swp_test.h ---- */
/* --------------------- */

#ifndef __SWP_TEST_H__
#define __SWP_TEST_H__

#ifdef __cplusplus
#ifdef PRAGMA_VERBOSE
#pragma message ("C++")
#endif
extern "C" {
#endif

int test_swp(int argc, char* argv[]);
void line_max3_ui8matrix_basic_swp(uint8 **X_Pack, int i, int j0, int j1, uint8 **Y_Pack);
void line_min3_ui8matrix_basic_swp(uint8 **X_Pack, int i, int j0, int j1, uint8 **Y_Pack);
void line_min3_ui8matrix_basic_swp_32(uint32 **X_Pack, int i, int j0, int j1, uint32 **Y_Pack);
void line_max3_ui8matrix_basic_swp_32(uint32 **X_Pack, int i, int j0, int j1, uint32 **Y_Pack);
void line_min3_ui8matrix_elu2_red_swp_32(uint32 **X_Pack, int i, int j0, int j1, uint32 **Y_Pack);
void line_max3_ui8matrix_elu2_red_swp_32(uint32 **X_Pack, int i, int j0, int j1, uint32 **Y_Pack);

void max3_ui8matrix_basic_swp(uint8 **X_Pack, int i0, int i1, int j0, int j1, uint8 **Y_Pack);
void min3_ui8matrix_basic_swp(uint8 **X_Pack, int i0, int i1, int j0, int j1, uint8 **Y_Pack);
void min3_ui8matrix_basic_swp_32(uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **Y_Pack);
void max3_ui8matrix_basic_swp_32(uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **Y_Pack);
void min3_ui8matrix_elu2_red_swp_32(uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **Y_Pack);
void max3_ui8matrix_elu2_red_swp_32(uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **Y_Pack);

void ouverture3_ui8matrix_pipeline_swp_8(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y);
void ouverture3_ui8matrix_pipeline_swp_32(uint32 **X, int i0, int i1, int j0, int j1, uint32 **T, uint32 **Y);
void ouverture3_ui8matrix_pipeline_elu2_red_swp_32(uint32 **X, int i0, int i1, int j0, int j1, uint32 **T, uint32 **Y);




#ifdef __cplusplus
}
#endif

#endif // __MOTION_TEST_H__

