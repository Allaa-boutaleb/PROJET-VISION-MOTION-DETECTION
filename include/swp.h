/* ------------- */
/* --- swp.h --- */
/* ------------- */

/*
 * Copyright (c) 2021 - 2021, Lionel Lacassagne, All rights reserved
 * Laboratoire de Recherche de Paris6 (LIP6), Sorbonne University, CNRS
 */

#ifndef __SWP_TEST_H__
#define __SWP_TEST_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

//Macro left/right ici
#define ui8left1(a, b) (((a)>>7) | (b << 1)) 
#define ui8right1(b, c) ((b >> 1) | (c << 7))


uint8   ui8left( uint8 a,  uint8 b, int k);
uint16 ui16left(uint16 a, uint16 b, int k);
uint32 ui32left(uint32 a, uint32 b, int k);

uint8   ui8right( uint8 b,  uint8 c, int k);
uint16 ui16right(uint16 b, uint16 c, int k);
uint32 ui32right(uint32 b, uint32 c, int k);

void printfR8(uint8_t x8r, char* name);
void printfM8(uint8_t x8m, char* name);

void printfR16(uint16_t x16r, char* name);
void printfM16(uint16_t x16m, char* name);

void printfR(uint32_t x, int n, char* name);
void printfM(uint32_t x, int n, char* name);

uint32_t str2Reg(char* str);
uint32_t str2Mem(char* str);

void set_ui8vector_str(uint8 *v, int i0, int i1, char *str);

int ishexa(char x);
// extension of isdigit to hexa number

int atob(char x); // atoi for bit
int atox(char x); // atoi for hexa (return 0 otherwise)

void set_ui8vector_str1 (uint8* v, int i0, int i1, char *str);
void set_ui8vector_str16(uint8* v, int i0, int i1, char *str);

void pack_ui4vector (uint8 *X1, int src_width, uint8  *Y4 );
void pack_ui8vector (uint8 *X1, int src_width, uint8  *Y8 );
void pack_ui16vector(uint8 *X1, int src_width, uint16 *Y16);
void pack_ui32vector(uint8 *X1, int src_width, uint32 *Y32);

void unpack_ui4vector (uint8  *X4 , int src_width, uint8* Y1);
void unpack_ui8vector (uint8  *X8 , int src_width, uint8* Y1);
void unpack_ui16vector(uint16 *X16, int src_width, uint8* Y1);
void unpack_ui32vector(uint32 *X32, int src_width, uint8* Y1);


void pack_ui4matrix (uint8  **X1, int src_height, int src_width, uint8  **Y4 );
void pack_ui8matrix (uint8  **X1, int src_height, int src_width, uint8  **Y8 );
void pack_ui16matrix(uint8  **X1, int src_height, int src_width, uint16 **Y16);
void pack_ui32matrix(uint8  **X1, int src_height, int src_width, uint32 **Y32);

void unpack_ui4matrix (uint8  **X4 , int src_height, int src_width, uint8** Y1);
void unpack_ui8matrix (uint8  **X8 , int src_height, int src_width, uint8** Y1);
void unpack_ui16matrix(uint16 **X16, int src_height, int src_width, uint8** Y1);
void unpack_ui32matrix(uint32 **X32, int src_height, int src_width, uint8** Y1);

void displayR_ui4vector (uint8  *v, int i0, int i1, char *name);
void displayR_ui8vector (uint8  *v, int i0, int i1, char *name);

void displayM_ui4vector (uint8  *v, int i0, int i1, char *name);
void displayM_ui8vector (uint8  *v, int i0, int i1, char *name);
void displayM_ui16vector(uint16 *v, int i0, int i1, char *name);
void displayM_ui32vector(uint32 *v, int i0, int i1, char *name);

void displayR_ui4matrix (uint8  **m, int i0, int i1, int j0, int j1, char *name);
void displayR_ui8matrix (uint8  **m, int i0, int i1, int j0, int j1, char *name);

void displayM_ui4matrix (uint8  **m, int i0, int i1, int j0, int j1, char *name);
void displayM_ui8matrix (uint8  **m, int i0, int i1, int j0, int j1, char *name);
void displayM_ui16matrix(uint16 **m, int i0, int i1, int j0, int j1, char *name);
void displayM_ui32matrix(uint32 **m, int i0, int i1, int j0, int j1, char *name);

#ifdef __cplusplus
}
#endif

#endif // __SWP_H__
