/* --------------------- */
/* ---- swp_test.c ---- */
/* --------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"

#include "swp.h"
#include "swp_test.h"

#include "macro_bench.h"
#include "x86intrin.h" // _rdtsc()

void random_ui8matrix(uint8 **X, int i0, int i1, int j0, int j1, int density)
// -------------------------------------------------------------------------
{
    // binary generator

    for (int i = i0; i <= i1; i++)
    {
        for (int j = j0; j <= j1; j++)
        {

            int r = rand() % 100;
            if (r <= density)
                X[i][j] = 1;
            else
                X[i][j] = 0;
        }
    }
}

// ---------------------------------------
void traitement_bordure_8(uint8 **X, int i0, int i1, int j0, int j1, int b)
// ---------------------------------------
{
    for (int k = 1; k < b; k++)
    {
        store2(X, i0 - k, j0 - k,  load2(X, i0, j0));       
        store2(X, i0 - k, j1 + k,  load2(X, i0, j1));
        store2(X, i1 + k, j0 - k,  load2(X, i1, j0));
        store2(X, i1 + k, j1 + k,  load2(X, i1, j1));

        for (int j = j0 - k + 1; j <= j1 + k - 1; j++)
        {
            store2(X, i0 - k, j, load2(X, i0 - k + 1, j));
            store2(X, i1 + k, j, load2(X, i1 + k - 1, j));
        }

        for (int i = i0 - k + 1; i <= i1 + k - 1; i ++)
        {
            store2(X, i, j0 - k, load2(X, i, j0 - k + 1));
            store2(X, i, j1 + k, load2(X, i, j1 + k - 1));
        }
    }
}

// ---------------------------------------
void traitement_bordure_16(uint16 **X, int i0, int i1, int j0, int j1, int b)
// ---------------------------------------
{
    for (int k = 1; k < b; k++)
    {
        store2(X, i0 - k, j0 - k,  load2(X, i0, j0));       
        store2(X, i0 - k, j1 + k,  load2(X, i0, j1));
        store2(X, i1 + k, j0 - k,  load2(X, i1, j0));
        store2(X, i1 + k, j1 + k,  load2(X, i1, j1));

        for (int j = j0 - k + 1; j <= j1 + k - 1; j++)
        {
            store2(X, i0 - k, j, load2(X, i0 - k + 1, j));
            store2(X, i1 + k, j, load2(X, i1 + k - 1, j));
        }

        for (int i = i0 - k + 1; i <= i1 + k - 1; i ++)
        {
            store2(X, i, j0 - k, load2(X, i, j0 - k + 1));
            store2(X, i, j1 + k, load2(X, i, j1 + k - 1));
        }
    }
}

// ---------------------------------------
void traitement_bordure_32(uint32 **X, int i0, int i1, int j0, int j1, int b)
// ---------------------------------------
{
    for (int k = 1; k < b; k++)
    {
        store2(X, i0 - k, j0 - k,  load2(X, i0, j0));       
        store2(X, i0 - k, j1 + k,  load2(X, i0, j1));
        store2(X, i1 + k, j0 - k,  load2(X, i1, j0));
        store2(X, i1 + k, j1 + k,  load2(X, i1, j1));

        for (int j = j0 - k + 1; j <= j1 + k - 1; j++)
        {
            store2(X, i0 - k, j, load2(X, i0 - k + 1, j));
            store2(X, i1 + k, j, load2(X, i1 + k - 1, j));
        }

        for (int i = i0 - k + 1; i <= i1 + k - 1; i ++)
        {
            store2(X, i, j0 - k, load2(X, i, j0 - k + 1));
            store2(X, i, j1 + k, load2(X, i, j1 + k - 1));
        }
    }
}

// -----------------------------------------------------------------------------
void line_max3_ui8matrix_basic_swp(uint8 **X_Pack, int i, int j0, int j1, uint8 **Y_Pack)
// -----------------------------------------------------------------------------
//  t1  a1 a2 a3 a4 a5 a6 a7 a8   t4
//  t2  x1 x2 x3 x4 x5 x6 x7 x8   t5
//  t3  b1 b2 b3 b4 b5 b6 b7 b8   t6
{
    uint8 a, x, b;
    uint8 t1, t2, t3, t4, t5, t6;
    uint8 r1, r2, r3, resu;

    for (int j = j0; j <= j1; j += 1)
    {
        a = load2(X_Pack, i + 1, j);
        x = load2(X_Pack, i, j);
        b = load2(X_Pack, i - 1, j);

        t1 = load2(X_Pack, i - 1, j - 1);
        t2 = load2(X_Pack, i, j - 1);
        t3 = load2(X_Pack, i - 1, j - 1);

        t4 = load2(X_Pack, i - 1, j + 1);
        t5 = load2(X_Pack, i, j + 1);
        t6 = load2(X_Pack, i - 1, j + 1);

        r1 = max3(t1, t2, t3);
        r2 = max3(a, x, b);
        r3 = max3(t4, t5, t6);

        r3 = ui8left1(r3, r2);
        r1 = ui8right1(r2, r1);

        resu = max3(r1, r2, r3);

        store2(Y_Pack, i, j, resu);
    }
}

// ----------------------------------------------------------------------------
void max3_ui8matrix_basic_swp(uint8 **X_Pack, int i0, int i1, int j0, int j1, uint8 **Y_Pack)
// ----------------------------------------------------------------------------
{

    for (int i = i0; i <= i1; i++)
    {
        line_max3_ui8matrix_basic_swp(X_Pack, i, j0, j1, Y_Pack);
    }
}

// -----------------------------------------------------------------------------
void line_min3_ui8matrix_basic_swp(uint8 **X_Pack, int i, int j0, int j1, uint8 **Y_Pack)
// -----------------------------------------------------------------------------
//  t1  b1 b2 b3 b4 b5 b6 b7 b8   t4
//  t2  x1 x2 x3 x4 x5 x6 x7 x8   t5
//  t3  a1 a2 a3 a4 a5 a6 a7 a8   t6
{
    uint8 a, x, b;
    uint8 t1, t2, t3, t4, t5, t6;
    uint8 r1, r2, r3, resu;

    for (int j = j0; j <= j1; j += 1)
    {
        a = load2(X_Pack, i + 1, j);
        x = load2(X_Pack, i, j);
        b = load2(X_Pack, i - 1, j);

        t1 = load2(X_Pack, i - 1, j - 1);
        t2 = load2(X_Pack, i, j - 1);
        t3 = load2(X_Pack, i - 1, j - 1);

        t4 = load2(X_Pack, i - 1, j + 1);
        t5 = load2(X_Pack, i, j + 1);
        t6 = load2(X_Pack, i - 1, j + 1);

        r1 = min3(t1, t2, t3);
        r2 = min3(a, x, b);
        r3 = min3(t4, t5, t6);

        r3 = ui8left1(r3, r2);
        r1 = ui8right1(r2, r1);

        resu = min3(r1, r2, r3);

        store2(Y_Pack, i, j, resu);
    }
}

// ----------------------------------------------------------------------------
void min3_ui8matrix_basic_swp(uint8 **X_Pack, int i0, int i1, int j0, int j1, uint8 **Y_Pack)
// ----------------------------------------------------------------------------
{

    for (int i = i0; i <= i1; i++)
    {
        line_min3_ui8matrix_basic_swp(X_Pack, i, j0, j1, Y_Pack);
    }
}

// -----------------------------------------------------------------------------
void line_max3_ui8matrix_basic_swp_16(uint16 **X_Pack, int i, int j0, int j1, uint16 **Y_Pack)
// -----------------------------------------------------------------------------
//  t1  a1 a2 a3 a4 a5 a6 a7 a8   t4
//  t2  x1 x2 x3 x4 x5 x6 x7 x8   t5
//  t3  b1 b2 b3 b4 b5 b6 b7 b8   t6
{
    uint16 a, x, b;
    uint16 t1, t2, t3, t4, t5, t6;
    uint16 r1, r2, r3, resu;

    for (int j = j0; j <= j1; j += 1)
    {
        a = load2(X_Pack, i + 1, j);
        x = load2(X_Pack, i, j);
        b = load2(X_Pack, i - 1, j);

        t1 = load2(X_Pack, i - 1, j - 1);
        t2 = load2(X_Pack, i, j - 1);
        t3 = load2(X_Pack, i - 1, j - 1);

        t4 = load2(X_Pack, i - 1, j + 1);
        t5 = load2(X_Pack, i, j + 1);
        t6 = load2(X_Pack, i - 1, j + 1);

        r1 = max3(t1, t2, t3);
        r2 = max3(a, x, b);
        r3 = max3(t4, t5, t6);

        r3 = ui8left1(r3, r2);
        r1 = ui8right1(r2, r1);

        resu = max3(r1, r2, r3);

        store2(Y_Pack, i, j, resu);
    }
}

// ----------------------------------------------------------------------------
void max3_ui8matrix_basic_swp_16(uint16 **X_Pack, int i0, int i1, int j0, int j1, uint16 **Y_Pack)
// ----------------------------------------------------------------------------
{

    for (int i = i0; i <= i1; i++)
    {
        line_max3_ui8matrix_basic_swp_16(X_Pack, i, j0, j1, Y_Pack);
    }
}

// -----------------------------------------------------------------------------
void line_min3_ui8matrix_basic_swp_16(uint16 **X_Pack, int i, int j0, int j1, uint16 **Y_Pack)
// -----------------------------------------------------------------------------
//  t1  b1 b2 b3 b4 b5 b6 b7 b8   t4
//  t2  x1 x2 x3 x4 x5 x6 x7 x8   t5
//  t3  a1 a2 a3 a4 a5 a6 a7 a8   t6
{
    uint16 a, x, b;
    uint16 t1, t2, t3, t4, t5, t6;
    uint16 r1, r2, r3, resu;

    for (int j = j0; j <= j1; j += 1)
    {
        a = load2(X_Pack, i + 1, j);
        x = load2(X_Pack, i, j);
        b = load2(X_Pack, i - 1, j);

        t1 = load2(X_Pack, i - 1, j - 1);
        t2 = load2(X_Pack, i, j - 1);
        t3 = load2(X_Pack, i - 1, j - 1);

        t4 = load2(X_Pack, i - 1, j + 1);
        t5 = load2(X_Pack, i, j + 1);
        t6 = load2(X_Pack, i - 1, j + 1);

        r1 = min3(t1, t2, t3);
        r2 = min3(a, x, b);
        r3 = min3(t4, t5, t6);

        r3 = ui8left1(r3, r2);
        r1 = ui8right1(r2, r1);

        resu = min3(r1, r2, r3);

        store2(Y_Pack, i, j, resu);
    }
}

// ----------------------------------------------------------------------------
void min3_ui8matrix_basic_swp_16(uint16 **X_Pack, int i0, int i1, int j0, int j1, uint16 **Y_Pack)
// ----------------------------------------------------------------------------
{

    for (int i = i0; i <= i1; i++)
    {
        line_min3_ui8matrix_basic_swp_16(X_Pack, i, j0, j1, Y_Pack);
    }
}

// -----------------------------------------------------------------------------
void line_min3_ui8matrix_basic_swp_32(uint32 **X_Pack, int i, int j0, int j1, uint32 **Y_Pack)
// -----------------------------------------------------------------------------
//  t1  a   t4
//  t2  x   t5
//  t3  b   t6
{
    uint32 a, x, b;
    uint32 t1, t2, t3, t4, t5, t6;
    uint32 r1, r2, r3, resu;

    for (int j = j0; j <= j1; j += 1)
    {
        a = load2(X_Pack, i + 1, j);
        x = load2(X_Pack, i, j);
        b = load2(X_Pack, i - 1, j);

        t1 = load2(X_Pack, i - 1, j - 1);
        t2 = load2(X_Pack, i, j - 1);
        t3 = load2(X_Pack, i - 1, j - 1);

        t4 = load2(X_Pack, i - 1, j + 1);
        t5 = load2(X_Pack, i, j + 1);
        t6 = load2(X_Pack, i - 1, j + 1);

        r1 = min3(t1, t2, t3);
        r2 = min3(a, x, b);
        r3 = min3(t4, t5, t6);

        r3 = ui8left1(r3, r2);
        r1 = ui8right1(r2, r1);

        resu = min3(r1, r2, r3);

        store2(Y_Pack, i, j, resu);
    }
}

// ----------------------------------------------------------------------------
void min3_ui8matrix_basic_swp_32(uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **Y_Pack)
// ----------------------------------------------------------------------------
{

    for (int i = i0; i <= i1; i++)
    {
        line_min3_ui8matrix_basic_swp_32(X_Pack, i, j0, j1, Y_Pack);
    }
}

// -----------------------------------------------------------------------------
void line_max3_ui8matrix_basic_swp_32(uint32 **X_Pack, int i, int j0, int j1, uint32 **Y_Pack)
// -----------------------------------------------------------------------------
//  t1  a   t4
//  t2  x   t5
//  t3  b   t6
{
    uint32 a, x, b;
    uint32 t1, t2, t3, t4, t5, t6;
    uint32 r1, r2, r3, resu;

    for (int j = j0; j <= j1; j += 1)
    {
        a = load2(X_Pack, i + 1, j);
        x = load2(X_Pack, i, j);
        b = load2(X_Pack, i - 1, j);

        t1 = load2(X_Pack, i - 1, j - 1);
        t2 = load2(X_Pack, i, j - 1);
        t3 = load2(X_Pack, i - 1, j - 1);

        t4 = load2(X_Pack, i - 1, j + 1);
        t5 = load2(X_Pack, i, j + 1);
        t6 = load2(X_Pack, i - 1, j + 1);

        r1 = max3(t1, t2, t3);
        r2 = max3(a, x, b);
        r3 = max3(t4, t5, t6);

        r3 = ui8left1(r3, r2);
        r1 = ui8right1(r2, r1);

        resu = max3(r1, r2, r3);

        store2(Y_Pack, i, j, resu);
    }
}

// ----------------------------------------------------------------------------
void max3_ui8matrix_basic_swp_32(uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **Y_Pack)
// ----------------------------------------------------------------------------
{

    for (int i = i0; i <= i1; i++)
    {
        line_max3_ui8matrix_basic_swp_32(X_Pack, i, j0, j1, Y_Pack);
    }
}

// ----------------------------------------------------------------------------
void line_min3_ui8matrix_elu2_red_swp_32(uint32 **X_Pack, int i, int j0, int j1, uint32 **Y_Pack)
// ----------------------------------------------------------------------------
//  t1  t2   t3
//  t4  x1   t5
//  t6  x2   t7
//  t8  t9   t10
{
    uint32 x1, x2;
    uint32 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    uint32 r1, r2, r3, r4, r5, r6;
    uint32 tmp1, tmp2, resu;

    t1 = load2(X_Pack, i - 1, j0 - 1);
    t2 = load2(X_Pack, i - 1, j0);
    t4 = load2(X_Pack, i, j0 - 1);
    x1 = load2(X_Pack, i, j0);
    t6 = load2(X_Pack, i + 1, j0 - 1);
    x2 = load2(X_Pack, i + 1, j0);
    t8 = load2(X_Pack, i + 2, j0 - 1);
    t9 = load2(X_Pack, i + 2, j0);

    r1 = min3(t1, t4, t6);
    r2 = min3(t2, x1, x2);
    r4 = min3(t4, t6, t8);
    r5 = min3(x1, x2, t9);

    for (int j = j0; j <= j1; j += 1)
    {
        t3 = load2(X_Pack, i - 1, j0 + 1);
        t5 = load2(X_Pack, i, j0 + 1);
        t7 = load2(X_Pack, i + 1, j0 + 1);
        t10 = load2(X_Pack, i + 2, j0 + 1);

        r3 = min3(t3, t5, t7);
        r6 = min3(t5, t7, t10);

        tmp1 = ui8left1(r3, r2);
        tmp2 = ui8right1(r2, r1);
        resu = min3(tmp1, r2, tmp2);
        store2(Y_Pack, i, j, resu);

        tmp1 = ui8left1(r6, r5);
        tmp2 = ui8right1(r5, r4);
        resu = min3(tmp1, r5, tmp2);
        store2(Y_Pack, i + 1, j, resu);

        r1 = r2;
        r2 = r3;
        r4 = r5;
        r5 = r6;
    }
}

// ----------------------------------------------------------------------------
void min3_ui8matrix_elu2_red_swp_32(uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **Y_Pack)
// ----------------------------------------------------------------------------
{
    int r = (i1 - i0 + 1) % 2;

    for (int i = i0; i <= i1 - r; i += 2)
    {
        line_min3_ui8matrix_elu2_red_swp_32(X_Pack, i, j0, j1, Y_Pack);
    }

    if (r)
        line_max3_ui8matrix_basic_swp_32(X_Pack, i1, j0, j1, Y_Pack);
}

// ----------------------------------------------------------------------------
void line_max3_ui8matrix_elu2_red_swp_32(uint32 **X_Pack, int i, int j0, int j1, uint32 **Y_Pack)
// ----------------------------------------------------------------------------
//  t1  t2   t3
//  t4  x1   t5
//  t6  x2   t7
//  t8  t9   t10
{
    uint32 x1, x2;
    uint32 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    uint32 r1, r2, r3, r4, r5, r6;
    uint32 tmp1, tmp2, resu;

    t1 = load2(X_Pack, i - 1, j0 - 1);
    t2 = load2(X_Pack, i - 1, j0);
    t4 = load2(X_Pack, i, j0 - 1);
    x1 = load2(X_Pack, i, j0);
    t6 = load2(X_Pack, i + 1, j0 - 1);
    x2 = load2(X_Pack, i + 1, j0);
    t8 = load2(X_Pack, i + 2, j0 - 1);
    t9 = load2(X_Pack, i + 2, j0);

    r1 = max3(t1, t4, t6);
    r2 = max3(t2, x1, x2);
    r4 = max3(t4, t6, t8);
    r5 = max3(x1, x2, t9);

    for (int j = j0; j <= j1; j += 1)
    {
        t3 = load2(X_Pack, i - 1, j0 + 1);
        t5 = load2(X_Pack, i, j0 + 1);
        t7 = load2(X_Pack, i + 1, j0 + 1);
        t10 = load2(X_Pack, i + 2, j0 + 1);

        r3 = max3(t3, t5, t7);
        r6 = max3(t5, t7, t10);

        tmp1 = ui8left1(r3, r2);
        tmp2 = ui8right1(r2, r1);
        resu = max3(tmp1, r2, tmp2);
        store2(Y_Pack, i, j, resu);

        tmp1 = ui8left1(r6, r5);
        tmp2 = ui8right1(r5, r4);
        resu = max3(tmp1, r5, tmp2);
        store2(Y_Pack, i + 1, j, resu);

        r1 = r2;
        r2 = r3;
        r4 = r5;
        r5 = r6;
    }
}

// ----------------------------------------------------------------------------
void max3_ui8matrix_elu2_red_swp_32(uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **Y_Pack)
// ----------------------------------------------------------------------------
{
    int r = (i1 - i0 + 1) % 2;

    for (int i = i0; i <= i1 - r; i += 2)
    {
        line_max3_ui8matrix_elu2_red_swp_32(X_Pack, i, j0, j1, Y_Pack);
    }

    if (r)
        line_max3_ui8matrix_basic_swp_32(X_Pack, i1, j0, j1, Y_Pack);
}

// -----------------------------------------------------------------------------
void test_swp_routine_16(int h, int w)
// -----------------------------------------------------------------------------
{
    int dh = 4;
    int dw = 3;

    int b = 1; // bord

    int w16 = w / 16;
    if (w % 16)
        w16 = w16 + 1;
    int w1 = 16 * w16; // w1 >= w

    puts("--------------------------------------------------");
    printf("test_morpho_min_routine h = %d w0 = %d w8 = %d, w1 = %d\n", h, w, w16, w1);
    if (w1 > w)
        puts("w1 > w0");

    uint8 **Y, **X;
    uint16 **X_Packed, **Y_Packed;

    int c; // error

    // puts("malloc");
    X = ui8matrix(0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    Y = ui8matrix(0, h - 1, 0, w1 - 1);
    X_Packed = ui16matrix(0 - b, h - 1 + b, 0 - b, w16 - 1 + b);
    Y_Packed = ui16matrix(0, h - 1, 0, w16 - 1);

    zero_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    zero_ui16matrix(X_Packed, 0, h - 1, 0, w16 - 1);
    zero_ui16matrix(Y_Packed, 0, h - 1, 0, w16 - 1);
    zero_ui8matrix(Y, 0, h - 1, 0, w1 - 1);

    random_ui8matrix(X, 0, h - 1, 0, w1 - 1, 80); // binaire [0,1]
    traitement_bordure_8(X, 0, h - 1, 0, w1 - 1, b);
    display_ui8matrix(X, 0, h - 1, 0, w1 - 1, "%5d", "X0");

    pack_ui16matrix(X, h, w1, X_Packed);
    traitement_bordure_16(X_Packed, 0, h - 1, 0, w16 - 1, b);
    displayM_ui16matrix(X_Packed, 0, h - 1, 0, w16 - 1, "X_PACKED");
    // unpack_ui32matrix(X_Packed, h, w1, X);
    // display_ui8matrix(X, 0, h - 1, 0, w1 - 1, "%5d", "X_UNPACKED");

    min3_ui8matrix_basic_swp_16(X_Packed, 0, h - 1, 0, w16 - 1, Y_Packed);
    displayM_ui16matrix(Y_Packed, 0, h - 1, 0, w16 - 1, "MIN_PACKED");
    unpack_ui16matrix(Y_Packed, h, w16, Y);
    display_ui8matrix(Y, 0, h - 1, 0, w1 - 1, "%5d", "MIN");

    max3_ui8matrix_basic_swp_16(X_Packed, 0, h - 1, 0, w16 - 1, Y_Packed);
    displayM_ui16matrix(Y_Packed, 0, h - 1, 0, w16 - 1, "MAX_PACKED");
    unpack_ui16matrix(Y_Packed, h, w16, Y);
    display_ui8matrix(Y, 0, h - 1, 0, w1 - 1, "%5d", "MAX");

    free_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    free_ui16matrix(X_Packed, 0 - b, h - 1 + b, 0 - b, w16 - 1 + b);
    free_ui8matrix(Y, 0, h - 1, 0, w1 - 1);
    free_ui16matrix(Y_Packed, 0, h - 1, 0, w16 - 1);
}
// -----------------------------------------------------------------------------
void test_swp_routine_32(int h, int w)
// -----------------------------------------------------------------------------
{
    int dh = 4;
    int dw = 3;

    int b = 1; // bord

    int w32 = w / 32;
    if (w % 32)
        w32 = w32 + 1;
    int w1 = 32 * w32; // w1 >= w

    puts("--------------------------------------------------");
    printf("test_morpho_min_routine h = %d w0 = %d w8 = %d, w1 = %d\n", h, w, w32, w1);
    if (w1 > w)
        puts("w1 > w0");

    uint8 **Y, **X;
    uint32 **X_Packed, **Y_Packed;

    int c; // error

    // puts("malloc");
    X = ui8matrix(0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    Y = ui8matrix(0, h - 1, 0, w1 - 1);
    X_Packed = ui32matrix(0 - b, h - 1 + b, 0 - b, w32 - 1 + b);
    Y_Packed = ui32matrix(0, h - 1, 0, w32 - 1);

    zero_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    zero_ui32matrix(X_Packed, 0, h - 1, 0, w32 - 1);
    zero_ui32matrix(Y_Packed, 0, h - 1, 0, w32 - 1);
    zero_ui8matrix(Y, 0, h - 1, 0, w1 - 1);

    random_ui8matrix(X, 0, h - 1, 0, w1 - 1, 80); // binaire [0,1]
    traitement_bordure_8(X, 0, h - 1, 0, w1 - 1, b);
    display_ui8matrix(X, 0, h - 1, 0, w1 - 1, "%5d", "X0");

    pack_ui32matrix(X, h, w1, X_Packed);
    traitement_bordure_32(X_Packed, 0, h - 1, 0, w32 - 1, b);
    displayM_ui32matrix(X_Packed, 0, h - 1, 0, w32 - 1, "X_PACKED");
    // unpack_ui32matrix(X_Packed, h, w1, X);
    // display_ui8matrix(X, 0, h - 1, 0, w1 - 1, "%5d", "X_UNPACKED");

    // min3_ui8matrix_basic_swp_32(X_Packed, 0, h - 1, 0, w32 - 1, Y_Packed);
    min3_ui8matrix_elu2_red_swp_32(X_Packed, 0, h - 1, 0, w32 - 1, Y_Packed);
    displayM_ui32matrix(Y_Packed, 0, h - 1, 0, w32 - 1, "MIN_PACKED");
    unpack_ui32matrix(Y_Packed, h, w32, Y);
    display_ui8matrix(Y, 0, h - 1, 0, w1 - 1, "%5d", "MIN");

    // max3_ui8matrix_basic_swp_32(X_Packed, 0, h - 1, 0, w32 - 1, Y_Packed);
    max3_ui8matrix_elu2_red_swp_32(X_Packed, 0, h - 1, 0, w32 - 1, Y_Packed);
    displayM_ui32matrix(Y_Packed, 0, h - 1, 0, w32 - 1, "MAX_PACKED");
    unpack_ui32matrix(Y_Packed, h, w32, Y);
    display_ui8matrix(Y, 0, h - 1, 0, w1 - 1, "%5d", "MAX");

    free_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    free_ui32matrix(X_Packed, 0 - b, h - 1 + b, 0 - b, w32 - 1 + b);
    free_ui8matrix(Y, 0, h - 1, 0, w1 - 1);
    free_ui32matrix(Y_Packed, 0, h - 1, 0, w32 - 1);
}

// -----------------------------------------------------------------------------
void test_swp_routine_8(int h, int w)
// -----------------------------------------------------------------------------
{
    int dh = 4;
    int dw = 3;

    int b = 1; // bord

    int w8 = w / 8;
    if (w % 8)
        w8 = w8 + 1;
    int w1 = 8 * w8; // w1 >= w
    // int wp = (w1) / 8;

    puts("--------------------------------------------------");
    printf("test_morpho_min_routine h = %d w0 = %d w8 = %d, w1 = %d\n", h, w, w8, w1);
    if (w1 > w)
        puts("w1 > w0");
    uint8 **X, **X_Packed;
    uint8 **Y, **Y_Packed;

    int c; // error

    // puts("malloc");
    X = ui8matrix(0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    Y = ui8matrix(0, h - 1, 0, w1 - 1);
    X_Packed = ui8matrix(0 - b, h - 1 + b, 0 - b, w8 - 1 + b);
    Y_Packed = ui8matrix(0, h - 1, 0, w8 - 1);

    zero_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    zero_ui8matrix(X_Packed, 0, h - 1, 0, w8 - 1);
    zero_ui8matrix(Y_Packed, 0, h - 1, 0, w8 - 1);
    zero_ui8matrix(Y, 0, h - 1, 0, w1 - 1);

    random_ui8matrix(X, 0, h - 1, 0, w1 - 1, 80); // binaire [0,1]
    traitement_bordure_8(X, 0, h - 1, 0, w1 - 1, b);
    display_ui8matrix(X, 0, h - 1, 0, w1 - 1, "%5d", "X0");

    pack_ui8matrix(X, h, w1, X_Packed);
    traitement_bordure_8(X_Packed, 0, h - 1, 0, w1 - 1, b);
    displayM_ui8matrix(X_Packed, 0, h - 1, 0, w8 - 1, "X_PACKED");

    min3_ui8matrix_basic_swp(X_Packed, 0, h - 1, 0, w8 - 1, Y_Packed);
    displayM_ui8matrix(Y_Packed, 0, h - 1, 0, w8 - 1, "MIN_PACKED");
    unpack_ui8matrix(Y_Packed, h, w8, Y);
    display_ui8matrix(Y, 0, h - 1, 0, w1 - 1, "%5d", "MIN");

    max3_ui8matrix_basic_swp(X_Packed, 0, h - 1, 0, w8 - 1, Y_Packed);
    displayM_ui8matrix(Y_Packed, 0, h - 1, 0, w8 - 1, "MAX_PACKED");
    unpack_ui8matrix(Y_Packed, h, w8, Y);
    display_ui8matrix(Y, 0, h - 1, 0, w1 - 1, "%5d", "MAX");

    free_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    free_ui8matrix(X_Packed, 0 - b, h - 1 + b, 0 - b, w8 - 1 + b);
    free_ui8matrix(Y, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_Packed, 0, h - 1, 0, w8 - 1);
}

// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_swp_8(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ------------------------------------------------------------------------------------------------------
{
    // Prologue
    line_min3_ui8matrix_basic_swp(X, i0 - 1, j0, j1, T);
    line_min3_ui8matrix_basic_swp(X, i0, j0, j1, T);

    for (int i = i0; i <= i1; i++)
    {
        line_min3_ui8matrix_basic_swp(X, i + 1, j0, j1, T);
        line_max3_ui8matrix_basic_swp(T, i, j0, j1, Y);
    }
}

// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_swp_16(uint16 **X, int i0, int i1, int j0, int j1, uint16 **T, uint16 **Y)
// ------------------------------------------------------------------------------------------------------
{
    // Prologue
    line_min3_ui8matrix_basic_swp_16(X, i0 - 1, j0, j1, T);
    line_min3_ui8matrix_basic_swp_16(X, i0, j0, j1, T);

    for (int i = i0; i <= i1; i++)
    {
        line_min3_ui8matrix_basic_swp_16(X, i + 1, j0, j1, T);
        line_max3_ui8matrix_basic_swp_16(T, i, j0, j1, Y);
    }
}

// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_swp_32(uint32 **X, int i0, int i1, int j0, int j1, uint32 **T, uint32 **Y)
// ------------------------------------------------------------------------------------------------------
{

    // Prologue
    line_min3_ui8matrix_basic_swp_32(X, i0 - 1, j0, j1, T);
    line_min3_ui8matrix_basic_swp_32(X, i0, j0, j1, T);

    for (int i = i0; i <= i1; i++)
    {
        line_min3_ui8matrix_basic_swp_32(X, i + 1, j0, j1, T);
        line_max3_ui8matrix_basic_swp_32(T, i, j0, j1, Y);
    }
}

// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_elu2_red_swp_32(uint32 **X, int i0, int i1, int j0, int j1, uint32 **T, uint32 **Y)
// ------------------------------------------------------------------------------------------------------
{
    // Prologue
    line_min3_ui8matrix_elu2_red_swp_32(X, i0 - 1, j0, j1, T);

    for (int i = i0; i <= i1; i += 2)
    {
        line_min3_ui8matrix_elu2_red_swp_32(X, i + 1, j0, j1, T);
        line_max3_ui8matrix_elu2_red_swp_32(T, i, j0, j1, Y);
    }
}

// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_swp_8_P(uint8 **X, uint8 **X_Pack, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y_Pack, uint8 **Y, int h, int w1, int w8)
// ------------------------------------------------------------------------------------------------------
{
    pack_ui8matrix(X, h, w1, X_Pack);

    // Prologue
    line_min3_ui8matrix_basic_swp(X_Pack, i0 - 1, j0, j1, T);
    line_min3_ui8matrix_basic_swp(X_Pack, i0, j0, j1, T);

    for (int i = i0; i <= i1; i++)
    {
        line_min3_ui8matrix_basic_swp(X_Pack, i + 1, j0, j1, T);
        line_max3_ui8matrix_basic_swp(T, i, j0, j1, Y_Pack);
    }

    unpack_ui8matrix(Y_Pack, h, w8, Y);
}

// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_swp_16_P(uint8 **X, uint16 **X_Pack, int i0, int i1, int j0, int j1, uint16 **T, uint16 **Y_Pack, uint8 **Y, int h, int w1, int w16)
// ------------------------------------------------------------------------------------------------------
{
    pack_ui16matrix(X, h, w1, X_Pack);

    // Prologue
    line_min3_ui8matrix_basic_swp_16(X_Pack, i0 - 1, j0, j1, T);
    line_min3_ui8matrix_basic_swp_16(X_Pack, i0, j0, j1, T);

    for (int i = i0; i <= i1; i++)
    {
        line_min3_ui8matrix_basic_swp_16(X_Pack, i + 1, j0, j1, T);
        line_max3_ui8matrix_basic_swp_16(T, i, j0, j1, Y_Pack);
    }
    unpack_ui16matrix(Y_Pack, h, w16, Y);
}

// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_swp_32_P(uint8 **X, uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **T, uint32 **Y_Pack, uint8 **Y, int h, int w1, int w32)
// ------------------------------------------------------------------------------------------------------
{
    pack_ui32matrix(X, h, w1, X_Pack);

    // Prologue
    line_min3_ui8matrix_basic_swp_32(X_Pack, i0 - 1, j0, j1, T);
    line_min3_ui8matrix_basic_swp_32(X_Pack, i0, j0, j1, T);

    for (int i = i0; i <= i1; i++)
    {
        line_min3_ui8matrix_basic_swp_32(X_Pack, i + 1, j0, j1, T);
        line_max3_ui8matrix_basic_swp_32(T, i, j0, j1, Y_Pack);
    }
    unpack_ui32matrix(Y_Pack, h, w32, Y);
}

// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_elu2_red_swp_32_P(uint8 **X, uint32 **X_Pack, int i0, int i1, int j0, int j1, uint32 **T, uint32 **Y_Pack, uint8 **Y, int h, int w1, int w32)
// ------------------------------------------------------------------------------------------------------
{
    pack_ui32matrix(X, h, w1, X_Pack);
    // Prologue
    line_min3_ui8matrix_elu2_red_swp_32(X_Pack, i0 - 1, j0, j1, T);

    for (int i = i0; i <= i1; i += 2)
    {
        line_min3_ui8matrix_elu2_red_swp_32(X_Pack, i + 1, j0, j1, T);
        line_max3_ui8matrix_elu2_red_swp_32(T, i, j0, j1, Y_Pack);
    }
    unpack_ui32matrix(Y_Pack, h, w32, Y);
}

// ---------------------------------------------------
void bench_swp(int n0, int n1, int nstep)
// ---------------------------------------------------
{
    int r = 1;
    int h = n1; // max size
    int w = n1; // max size
    int w0 = w;
    int w8 = w0 / 8;
    if (w0 % 8)
        w8 = w8 + 1;
    int w1 = 8 * w8; // w1 >= w

    int w16 = w0 / 16;
    if (w0 % 16)
        w16 = w16 + 1;
    w1 = 16 * w16; // w1 >= w

    int w32 = w0 / 32;
    if (w0 % 32)
        w32 = w32 + 1;
    w1 = 32 * w32; // w1 >= w

    uint8 **X;
    uint8 **X8;
    uint16 **X16;
    uint32 **X32;

    uint8 **T8;
    uint16 **T16;
    uint32 **T32;

    uint8 **Y, **Y8;
    uint16 **Y16;
    uint32 **Y32, **Y32_elu2;

    double cpp8, cpp8p;
    double cpp16, cpp16p;
    double cpp32, cpp32p;
    double cpp32_red, cpp32_redp;

    char *format = "%8.2f";

    format = "%5.0f";
    format = "%6.1f";

    puts("=============================");
    puts("== bench_morpho_ouverture ===");
    puts("=============================");

    // puts("malloc");

    // X 2r-border
    X = ui8matrix(0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);
    X8 = ui8matrix(0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w8 - 1 + 2 * r);
    X16 = ui16matrix(0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w16 - 1 + 2 * r);
    X32 = ui32matrix(0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w32 - 1 + 2 * r);

    // T 1r-border
    T8 = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w8 - 1 + 1 * r);
    T16 = ui16matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w16 - 1 + 1 * r);
    T32 = ui32matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w32 - 1 + 1 * r);

    // Y 0r-border
    Y = ui8matrix(0, h - 1, 0, w1 - 1);

    Y8 = ui8matrix(0, h - 1, 0, w8 - 1);
    Y16 = ui16matrix(0, h - 1, 0, w16 - 1);
    Y32 = ui32matrix(0, h - 1, 0, w32 - 1);
    Y32_elu2 = ui32matrix(0, h - 1, 0, w32 - 1);

    zero_ui8matrix(X, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);

    zero_ui8matrix(T8, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w8 - 1 + 1 * r);
    zero_ui16matrix(T16, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w16 - 1 + 1 * r);
    zero_ui32matrix(T32, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w32 - 1 + 1 * r);

    zero_ui8matrix(Y, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y8, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w8 - 1 + 0 * r);
    zero_ui16matrix(Y16, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w16 - 1 + 0 * r);
    zero_ui32matrix(Y32, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w32 - 1 + 0 * r);
    zero_ui32matrix(Y32_elu2, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w32 - 1 + 0 * r);

    puts("temps de calcul en ccp (cycle/point)");

    for (int n = n0; n <= n1; n += nstep)
    {

        h = n;
        w8 = n / 8;
        w16 = n / 16;
        w32 = n / 32;
        w1 = n / 1;

        // printf("i = %3d\n", n);

        resize_ui8matrix(X, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);
        resize_ui8matrix(X8, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w8 - 1 + 2 * r);
        resize_ui16matrix(X16, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w16 - 1 + 2 * r);
        resize_ui32matrix(X32, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w32 - 1 + 2 * r);

        // T 1r-border
        resize_ui8matrix(T8, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w8 - 1 + 1 * r);
        resize_ui16matrix(T16, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w16 - 1 + 1 * r);
        resize_ui32matrix(T32, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w32 - 1 + 1 * r);

        // Y 0r-border
        resize_ui8matrix(Y, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y8, 0, h - 1, 0, w8 - 1);
        resize_ui16matrix(Y16, 0, h - 1, 0, w16 - 1);
        resize_ui32matrix(Y32, 0, h - 1, 0, w32 - 1);

        pack_ui8matrix(X, h, w1, X8);
        pack_ui16matrix(X, h, w1, X16);
        pack_ui32matrix(X, h, w1, X32);

        printf("n: %d\t", n);

        BENCH(ouverture3_ui8matrix_pipeline_swp_8(X8, 0, h - 1, 0, w8 - 1, T8, Y8), n, cpp8);
        BENCH(ouverture3_ui8matrix_pipeline_swp_32(X32, 0, h - 1, 0, w32 - 1, T32, Y32), n, cpp32);
        BENCH(ouverture3_ui8matrix_pipeline_swp_16(X16, 0, h - 1, 0, w16 - 1, T16, Y16), n, cpp16);
        BENCH(ouverture3_ui8matrix_pipeline_elu2_red_swp_32(X32, 0, h - 1, 0, w32 - 1, T32, Y32_elu2), n, cpp32_red);

        BENCH(ouverture3_ui8matrix_pipeline_swp_8_P(X, X8, 0, h - 1, 0, w8 - 1, T8, Y8, Y, h, w1, w8), n, cpp8p);
        BENCH(ouverture3_ui8matrix_pipeline_swp_16_P(X, X16, 0, h - 1, 0, w16 - 1, T16, Y16, Y, h, w1, w16), n, cpp16p);
        BENCH(ouverture3_ui8matrix_pipeline_swp_32_P(X, X32, 0, h - 1, 0, w32 - 1, T32, Y32, Y, h, w1, w32), n, cpp32p);
        BENCH(ouverture3_ui8matrix_pipeline_elu2_red_swp_32_P(X, X32, 0, h - 1, 0, w32 - 1, T32, Y32_elu2, Y, h, w1, w32), n, cpp32_redp);

        printf(format, cpp8);
        printf(format, cpp16);
        printf(format, cpp32);
        printf(format, cpp32_red);
        printf("\t");
        printf(format, cpp8p);
        printf(format, cpp16p);
        printf(format, cpp32p);
        printf(format, cpp32_redp);

        putchar('\n');
    }
}

void test_morpho_swp()
{
    puts("=== test_swp ===");
    int h0 = 8;
    int w0 = 16;

    int dh = 4;
    int dw = 3;
    for (int h = h0; h <= h0 + dh; h++)
    { // pour tester elu2
        for (int w = w0; w <= w0 + dw; w++)
        { // pour tester ilu3
            test_swp_routine_8(h, w);
            test_swp_routine_16(h, w);
            test_swp_routine_32(h, w);
        }
    }
}

void start_bench_swp()
{
    puts("=== bench_swp ===");
    bench_swp(128, 512, 8);
    // bench_swp(128, 1024, 8);
}

// ====================================
int test_swp(int argc, char *argv[])
// ====================================
{
    test_morpho_swp();
    // start_bench_swp();
    return 0;
}