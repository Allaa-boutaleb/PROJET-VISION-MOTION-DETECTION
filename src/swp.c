/* ------------- */
/* --- swp.c --- */
/* ------------- */

/*
 * Copyright (c) 2021 - 2021, Lionel Lacassagne, All rights reserved
 * LIP6 - Sorbonne University
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <ctype.h> // isdigit

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"




// fonctions left/right = versions lentes

// -----------------------------------
uint8 ui8left(uint8 a, uint8 b, int k)
// -----------------------------------
{
    return (a >> (8-k)) | (b << k);
}
// ---------------------------------------
uint16 ui16left(uint16 a, uint16 b, int k)
// ---------------------------------------
{
    return (a >> (16-k)) | (b << k);
}
// ---------------------------------------
uint32 ui32left(uint32 a, uint32 b, int k)
// ---------------------------------------
{
    return (a >> (32-k)) | (b << k);
}
// ------------------------------------
uint8 ui8right(uint8 b, uint8 c, int k)
// ------------------------------------
{
    return (b >> k) | (c << (8-k));
}
// ----------------------------------------
uint16 ui16right(uint16 b, uint16 c, int k)
// ----------------------------------------
{
    return (b >> k) | (c << (16-k));
}
// ----------------------------------------
uint32 ui32right(uint32 b, uint32 c, int k)
// ----------------------------------------
{
    return (b >> k) | (c << (32-k));
}
// -------------------
bit_t char2bit(char c)
// -------------------
{
    bit_t b;
    
    if(isdigit(c)) {
        uint8_t d = atoi((const char*) &c);
        if(d>1) d = 1; // saturation
        b = (bit_t) d;
    } else {
        b = 0;
    }
    return b;
}
// ------------------------
uint32_t str2Reg(char* str)
// ------------------------
{
    int len = strlen(str);
    int n = (len < 32) ? len : 32; // clamping au min
    
    uint32_t x = 0;
    
    for(int i=0; i<n; i++) {
    
        char c = str[i];
        
        bit_t b = char2bit(c);
        
        x = (x << 1) | b;
    }
    return x;
}
// ------------------------
uint32_t str2Mem(char* str)
// ------------------------
{
    int len = strlen(str);
    int n = (len < 32) ? len : 32; // clamping au min
    
    uint32_t x = 0;
    
    for(int i=0; i<n; i++) {
        
        char c = str[i];
        
        bit_t b = char2bit(c);
        
        x = x | (b << i);
    }
    return x;
}
// -----------------------------------------------------------
void set_ui8vector_str_v0(uint8 *v, int i0, int i1, char *str)
// -----------------------------------------------------------
{
    // unsafe
    char c[2];
    
    for(int i=i0; i<=i1; i++) {
        *c = *str++;
        if(isdigit(*c))
            v[i] = (uint8) atoi(c);
        else
            v[i] = (uint8) 0;
    } // i
}
// -----------------------------------------------------------
void set_ui8vector_str_v2(uint8 *v, int i0, int i1, char *str)
// -----------------------------------------------------------
{
    int len = strlen(str);
    int ncol = i1-i0+1;
    int m;
    char c[2];
    
    m = (len < ncol) ? len : ncol;
    
    for(int i=i0; i<=i0+m-1; i++) {
        *c = *str++;
        if(isdigit(*c))
            v[i] = (uint8) atoi(c);
        else
            v[i] = (uint8) 0;
    }
    
    for(int i=i0+m; i<=i1; i++) {
        v[i] = (uint8) 0;
    }
}
// --------------------------------------------------------
void set_ui8vector_str(uint8 *v, int i0, int i1, char *str)
// --------------------------------------------------------
{
    //set_ui8vector_str_v0(v, i0, i1, str);
    //set_ui8vector_str_v1(v, i0, i1, str);
    set_ui8vector_str_v2(v, i0, i1, str);
}
// ---------------
int ishexa(char x)
// ---------------
{
    if( ((x >= '0') && (x <= '9')) ||
       ((x >= 'a') && (x <= 'f')) ||
       ((x >= 'A') && (x <= 'F')) ) {
        return 1;
    } else {
        return 0;
    }
}
// -------------
int atob(char x)
// -------------
{
    // restriction of atoi for bit
    if( ((x >= '0') && (x <= '1')) ) return (x - '0');
    return 0; // default value
}
// -------------
int atox(char x)
// -------------
{
    // extension of atoi for hexa number
    if( ((x >= '0') && (x <= '9')) ) return (x - '0');
    if( ((x >= 'a') && (x <= 'f')) ) return (x - 'a' + 10);
    if( ((x >= 'A') && (x <= 'F')) ) return (x - 'A' + 10);
    return 0; // default value
}
// ---------------------------------------------------------
void set_ui8vector_str1(uint8 *v, int i0, int i1, char *str)
// ---------------------------------------------------------
{
    int len = strlen(str);
    int ncol = i1-i0+1;
    int m;
    char c;
    
    m = (len < ncol) ? len : ncol;
    
    for(int i=i0; i<=i0+m-1; i++) {
        c = *str++;
        v[i] = (uint8) atob(c);
    }
    
    for(int i=i0+m; i<=i1; i++) {
        v[i] = (uint8) 0;
    }
}
// ---------------------------------------------------------
void set_ui8vector_str16(uint8 *v, int i0, int i1, char *str)
// ---------------------------------------------------------
{
    int len = strlen(str);
    int ncol = i1-i0+1;
    int m;
    char c;
    
    m = (len < ncol) ? len : ncol;
    
    for(int i=i0; i<=i0+m-1; i++) {
        c = *str++;
        v[i] = (uint8) atox(c);
    }
    
    for(int i=i0+m; i<=i1; i++) {
        v[i] = (uint8) 0;
    }
}
// -------------------------------------------------------
void unpack_ui4vector(uint8* X4, int src_width, uint8* Y1)
// -------------------------------------------------------
{
    for(int i=0; i<src_width; i++) {
    
        uint8 x = X4[i];
    
        for(int k=0; k<4; k++) {
            Y1[4*i + k] = x & 1;
            x = x >> 1;
        }
    }
}
// -------------------------------------------------------
void unpack_ui8vector(uint8* X8, int src_width, uint8* Y1)
// -------------------------------------------------------
{
    for(int i=0; i<src_width; i++) {
        uint8 x = X8[i];
        for(int k=0; k<8; k++) {
            Y1[8*i + k] = x & 1;
            x = x >> 1;
        }
    }
}
// ----------------------------------------------------------
void unpack_ui16vector(uint16* X16, int src_width, uint8* Y1)
// ----------------------------------------------------------
{
    for(int i=0; i<src_width; i++) {
        
        uint8 x = X16[i];
        
        for(int k=0; k<16; k++) {
            Y1[16*i + k] = x & 1;
            x = x >> 1;
        }
    }
}
// ----------------------------------------------------------
void unpack_ui32vector(uint32* X32, int src_width, uint8* Y1)
// ----------------------------------------------------------
{
    for(int i=0; i<src_width; i++) {
        
        uint32 x = X32[i];
        
        for(int k=0; k<32; k++) {
            Y1[32*i + k] = x & 1;
            x = x >> 1;
        }
    }
}
// -----------------------------------------------------
void pack_ui4vector(uint8* X1, int src_width, uint8* Y4)
// -----------------------------------------------------
{
    int dst_width = src_width/4;
    
    for(int i=0; i<dst_width; i++) {
    
        uint8 x = 0;
        
        for(int k=0; k<4; k++) {
            
            uint8 b = X1[4*i + k];
            x = x | (b << k);
        }
        Y4[i] = x;
    }
}
// -----------------------------------------------------
void pack_ui8vector(uint8* X1, int src_width, uint8* Y8)
// -----------------------------------------------------
{
    int dst_width = src_width/8;
    
    for(int i=0; i<dst_width; i++) {
        
        uint8 x = 0;
        
        for(int k=0; k<8; k++) {
            
            uint8 b = X1[8*i + k];
            x = x | (b << k);
        }
        Y8[i] = x;
    }
}
// --------------------------------------------------------
void pack_ui16vector(uint8* X1, int src_width, uint16* Y16)
// --------------------------------------------------------
{
    int dst_width = src_width/16;
    
    for(int i=0; i<dst_width; i++) {
        
        uint16 x = 0;
        
        for(int k=0; k<16; k++) {
            
            uint8 b = X1[16*i + k];
            x = x | (b << k);
        }
        Y16[i] = x;
    }
}
// --------------------------------------------------------
void pack_ui32vector(uint8* X1, int src_width, uint32* Y32)
// --------------------------------------------------------
{
    int dst_width = src_width/32;
    
    for(int i=0; i<dst_width; i++) {
        
        uint32 x = 0;
        
        for(int k=0; k<32; k++) {
            
            uint8 b = X1[32*i + k];
            x = x | (b << k);
        }
        Y32[i] = x;
    }
}
// -----------------------------------------------------------------------
void pack_ui4matrix(uint8** X1, int src_height, int src_width, uint8** Y8)
// -----------------------------------------------------------------------
{
    for(int i=0; i<src_height; i++) {
        pack_ui4vector(X1[i], src_width, Y8[i]);
    }
}
// -----------------------------------------------------------------------
void pack_ui8matrix(uint8** X1, int src_height, int src_width, uint8** Y8)
// -----------------------------------------------------------------------
{
    for(int i=0; i<src_height; i++) {
        pack_ui8vector(X1[i], src_width, Y8[i]);
    }
}
// --------------------------------------------------------------------------
void pack_ui16matrix(uint8** X1, int src_height, int src_width, uint16** Y16)
// --------------------------------------------------------------------------
{
    for(int i=0; i<src_height; i++) {
        pack_ui16vector(X1[i], src_width, Y16[i]);
    }
}
// --------------------------------------------------------------------------
void pack_ui32matrix(uint8** X1, int src_height, int src_width, uint32** Y32)
// --------------------------------------------------------------------------
{
    for(int i=0; i<src_height; i++) {
        pack_ui32vector(X1[i], src_width, Y32[i]);
    }
}
// -------------------------------------------------------------------------
void unpack_ui4matrix(uint8** X8, int src_height, int src_width, uint8** Y1)
// -------------------------------------------------------------------------
{
    for(int i=0; i<src_height; i++) {
        unpack_ui4vector(X8[i], src_width, Y1[i]);
    }
}
// -------------------------------------------------------------------------
void unpack_ui8matrix(uint8** X8, int src_height, int src_width, uint8** Y1)
// -------------------------------------------------------------------------
{
    for(int i=0; i<src_height; i++) {
        unpack_ui8vector(X8[i], src_width, Y1[i]);
    }
}
// ----------------------------------------------------------------------------
void unpack_ui16matrix(uint16** X16, int src_height, int src_width, uint8** Y1)
// ----------------------------------------------------------------------------
{
    for(int i=0; i<src_height; i++) {
        unpack_ui16vector(X16[i], src_width, Y1[i]);
    }
}
// ----------------------------------------------------------------------------
void unpack_ui32matrix(uint32** X32, int src_height, int src_width, uint8** Y1)
// ----------------------------------------------------------------------------
{
    for(int i=0; i<src_height; i++) {
        unpack_ui32vector(X32[i], src_width, Y1[i]);
    }
}
// ----------------------------------------
void printfR(uint32_t x, int n, char* name)
// ----------------------------------------
{
    // affichage en binaire "Register-like"
    if(name) printf("%s = (%x)h = (", name, x);
    
    for(int i=n; i>0; i--) {
        if( !(i&7) && (i != n) ) putchar(' ');
        bit_t b = (x >> (i-1)) & 1;
        printf("%d", b);
        
    }
    if(name) puts(")b");
}
// ----------------------------------------
void printfM(uint32_t x, int n, char* name)
// ----------------------------------------
{
    // affichage en binaire "Memory-like"
    if(name) printf("%s = (%x)h = (", name, x);
    
    for(int i=0; i<n; i++) {
        if( !(i&7) && i ) putchar(' ');
        bit_t b = (x >> i) & 1;
        printf("%d", b);
    }
    if(name) puts(")b");
}

// ----------------------------------------------------------
void displayR_ui4vector(uint8 *v, int i0, int i1, char *name)
// ----------------------------------------------------------
{
    if(name != NULL) printf("%s", name);
    
    for(int i=i0; i<=i1; i++) {
        printfR(v[i], 4, NULL); // 8-bit display
        putchar(' ');
    }
    if(name) putchar('\n');
}
// ----------------------------------------------------------
void displayR_ui8vector(uint8 *v, int i0, int i1, char *name)
// ----------------------------------------------------------
{
    if(name != NULL) printf("%s", name);
    
    for(int i=i0; i<=i1; i++) {
        printfR(v[i], 8, NULL); // 8-bit display
        putchar(' ');
    }
    if(name) putchar('\n');
}
// ----------------------------------------------------------
void displayM_ui4vector(uint8 *v, int i0, int i1, char *name)
// ----------------------------------------------------------
{
    if(name != NULL) printf("%s", name);
    
    for(int i=i0; i<=i1; i++) {
        printfM(v[i], 4, NULL); // 4-bit display
        putchar(' ');
    }
    if(name) putchar('\n');
}
// ----------------------------------------------------------
void displayM_ui8vector(uint8 *v, int i0, int i1, char *name)
// ----------------------------------------------------------
{
    if(name != NULL) printf("%s", name);
    
    for(int i=i0; i<=i1; i++) {
        printfM(v[i], 8, NULL); // 8-bit display
        putchar(' ');
    }
    if(name) putchar('\n');
}
// ------------------------------------------------------------
void displayM_ui16vector(uint16 *v, int i0, int i1, char *name)
// ------------------------------------------------------------
{
    if(name != NULL) printf("%s", name);
    
    for(int i=i0; i<=i1; i++) {
        printfM(v[i], 16, NULL); // 16-bit display
        putchar(' ');
    }
    if(name) putchar('\n');
}
// ------------------------------------------------------------
void displayM_ui32vector(uint32 *v, int i0, int i1, char *name)
// ------------------------------------------------------------
{
    if(name != NULL) printf("%s", name);
    
    for(int i=i0; i<=i1; i++) {
        printfM(v[i], 32, NULL); // 32-bit display
        putchar(' ');
    }
    if(name) putchar('\n');
}
// ---------------------------------------------------------------------------
void displayR_ui4matrix(uint8 **m, int i0, int i1, int j0, int j1, char *name)
// ---------------------------------------------------------------------------
{
    if(name != NULL) printf("%s\n", name);
    
    for(int i=i0; i<=i1; i++) {
        displayR_ui4vector(m[i], j0, j1, NULL);
        putchar('\n');
    }
    if(name) putchar('\n');
}
// ---------------------------------------------------------------------------
void displayR_ui8matrix(uint8 **m, int i0, int i1, int j0, int j1, char *name)
// ---------------------------------------------------------------------------
{
    if(name != NULL) printf("%s\n", name);
    
    for(int i=i0; i<=i1; i++) {
        displayR_ui8vector(m[i], j0, j1, NULL);
        putchar('\n');
    }
    if(name) putchar('\n');
}
// ---------------------------------------------------------------------------
void displayM_ui4matrix(uint8 **m, int i0, int i1, int j0, int j1, char *name)
// ---------------------------------------------------------------------------
{
    if(name != NULL) printf("%s\n", name);
    
    for(int i=i0; i<=i1; i++) {
        displayM_ui4vector(m[i], j0, j1, NULL);
        putchar('\n');
    }
    if(name) putchar('\n');
}
// ---------------------------------------------------------------------------
void displayM_ui8matrix(uint8 **m, int i0, int i1, int j0, int j1, char *name)
// ---------------------------------------------------------------------------
{
    if(name != NULL) printf("%s\n", name);
    
    for(int i=i0; i<=i1; i++) {
        displayM_ui8vector(m[i], j0, j1, NULL);
        putchar('\n');
    }
    if(name) putchar('\n');
}
// -----------------------------------------------------------------------------
void displayM_ui16matrix(uint16 **m, int i0, int i1, int j0, int j1, char *name)
// -----------------------------------------------------------------------------
{
    if(name != NULL) printf("%s\n", name);
    
    for(int i=i0; i<=i1; i++) {
        displayM_ui16vector(m[i], j0, j1, NULL);
        putchar('\n');
    }
    if(name) putchar('\n');
}
// -----------------------------------------------------------------------------
void displayM_ui32matrix(uint32 **m, int i0, int i1, int j0, int j1, char *name)
// -----------------------------------------------------------------------------
{
    if(name != NULL) printf("%s\n", name);
    
    for(int i=i0; i<=i1; i++) {
        displayM_ui32vector(m[i], j0, j1, NULL);
        putchar('\n');
    }
    if(name) putchar('\n');
}