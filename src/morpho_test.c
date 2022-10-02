/* --------------------- */
/* --- morpho_test.c --- */
/* --------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique
 */

/* -------------- */
/* --- main.c --- */
/* -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"

#include "macro_bench.h"
#include "x86intrin.h" // _rdtsc()

#include "swp.h"
#include "morpho_max.h"
#include "morpho_min.h"
#include "morpho_ouverture.h"

#include "morpho_test.h"

// -------------------------------------------------------------------------------
void rand_ui8matrix(uint8 **X, int i0, int i1, int j0, int j1, uint8 x0, uint8 x1)
// -------------------------------------------------------------------------------
{
    int d = x1 - x0 + 1; // range [x0, x1]

    for (int i = i0; i <= i1; i++)
    {
        for (int j = j0; j <= j1; j++)
        {

            int x = x0 + rand() % d;

            X[i][j] = (uint8)x;
        }
    }
}
// -------------------------------------------------------------------------
void rand1_ui8matrix(uint8 **X, int i0, int i1, int j0, int j1, int density)
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
// ------------------------------------------------------------------------------------
int compare_ui8matrix(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y, char *name)
// ------------------------------------------------------------------------------------
{
    int c = 0;

    // display_ui8matrix(Y, i0, i1, j0, j1, "%5d", name);

    for (int i = i0; i <= i1; i++)
    {
        for (int j = j0; j <= j1; j++)
        {
            if (X[i][j] != Y[i][j])
            {
                printf("%s [%3d][%3d] = %3d NEQ %3d\n", name, i, j, X[i][i], Y[i][i]);
                c++;
            }
        }
    }
    if (c)
    {
        printf("%d errors !\n", c);
        exit(-1);
    }
    else
    {
        printf("[compare]: %s is OK\n", name);
        return c;
    }
    puts("----------------------");
}

// ---------------------------------------
void traitement_bordure(uint8 **X, int i0, int i1, int j0, int j1, int b)
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
void test_morpho_max_routine(int h, int w0)
// ---------------------------------------
{
    int b = 1; // bord

    int w8 = w0 / 8;
    if (w0 % 8)
        w8 = w8 + 1;
    int w1 = 8 * w8; // w1 >= w

    puts("--------------------------------------------------");
    printf("test_morpho_max_routine h = %d w0 = %d w8 = %d, w1 = %d\n", h, w0, w8, w1);
    if (w1 > w0)
        puts("w1 > w0");
    uint8 **X;
    uint8 **Y_bas, **Y_reg, **Y_rot, **Y_red, **Y_ilu3, **Y_ilu3r, **Y_elu2r, **Y_elu2rf, **Y_ilu3_elu2r, **Y_ilu3_elu2rf;

    int c; // error

    // puts("malloc");
    X = ui8matrix(0 - b, h - 1 + b, 0 - b, w1 - 1 + b);

    Y_bas = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_reg = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_rot = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_ilu3 = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_ilu3r = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_elu2r = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_elu2rf = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_ilu3_elu2r = ui8matrix(0, h - 1, 0, w1 - 1);  /**/
    Y_ilu3_elu2rf = ui8matrix(0, h - 1, 0, w1 - 1); /**/

    // puts("zero");
    zero_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);

    zero_ui8matrix(Y_bas, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_reg, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_rot, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_red, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_ilu3, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_ilu3r, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_elu2r, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_elu2rf, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_ilu3_elu2r, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_ilu3_elu2rf, 0, h - 1, 0, w1 - 1);

    // puts("rand");

    // deux generateurs aleatoire au choix
    // rand_ui8matrix(X, 0, h-1, 0, w0-1, 1, 255); // niveau de gris [1,255]
    rand1_ui8matrix(X, 0, h - 1, 0, w0 - 1, 20); // binaire [0,1]: pourcentage de point a 1
    traitement_bordure(X, 0, h - 1,0, w0 - 1, b);

    // puts("display");
    // display_ui8matrix (X,  0, h-1, 0, w0-1, "%5d", "X0");

    // puts("\n-- max3 --");
    max3_ui8matrix_basic(X, 0, h - 1, 0, w0 - 1, Y_bas);                        // puts("bas");
    max3_ui8matrix_reg(X, 0, h - 1, 0, w0 - 1, Y_reg);                          // puts("reg");
    max3_ui8matrix_rot(X, 0, h - 1, 0, w0 - 1, Y_rot);                          // puts("rot");
    max3_ui8matrix_red(X, 0, h - 1, 0, w0 - 1, Y_red);                          // puts("red");
    max3_ui8matrix_ilu3(X, 0, h - 1, 0, w0 - 1, Y_ilu3);                        // puts("ilu3");
    max3_ui8matrix_ilu3_red(X, 0, h - 1, 0, w0 - 1, Y_ilu3r);                   // puts("ilu3_red");
    max3_ui8matrix_elu2_red(X, 0, h - 1, 0, w0 - 1, Y_elu2r);                   // puts("elu2_red");
    max3_ui8matrix_elu2_red_factor(X, 0, h - 1, 0, w0 - 1, Y_elu2rf);           // puts("elu2_red_factor");
    max3_ui8matrix_ilu3_elu2_red(X, 0, h - 1, 0, w0 - 1, Y_ilu3_elu2r);         // puts("ilu3_elu2_red_factor");
    max3_ui8matrix_ilu3_elu2_red_factor(X, 0, h - 1, 0, w0 - 1, Y_ilu3_elu2rf); // puts("ilu3_elu2_red_factor");
    // puts("done\n");

    // display_ui8matrix(Y_bas,         0, h-1, 0, w0-1, "%5d", "Y bas        ");
    // display_ui8matrix(Y_reg,         0, h-1, 0, w0-1, "%5d", "Y reg        ");
    // display_ui8matrix(Y_reg,         0, h-1, 0, w0-1, "%5d", "Y reg        ");
    // display_ui8matrix(Y_rot,         0, h-1, 0, w0-1, "%5d", "Y rot        ");
    // display_ui8matrix(Y_red,         0, h-1, 0, w0-1, "%5d", "Y red        ");
    // display_ui8matrix(Y_ilu3,        0, h-1, 0, w0-1, "%5d", "Y ilu3       ");
    // display_ui8matrix(Y_ilu3r,       0, h-1, 0, w0-1, "%5d", "Y ilu3r      ");
    // display_ui8matrix(Y_elu2r,       0, h-1, 0, w0-1, "%5d", "Y elu2r      ");
    // display_ui8matrix(Y_elu2rf,      0, h-1, 0, w0-1, "%5d", "Y elu2rf     ");
    // display_ui8matrix(Y_ilu3_elu2r,  0, h-1, 0, w0-1, "%5d", "Y ilu3 elu2r");
    // display_ui8matrix(Y_ilu3_elu2rf, 0, h-1, 0, w0-1, "%5d", "Y ilu3 elu2rf");

    // puts("-- compare max --");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_reg, "Y reg                       ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_rot, "Y rot                       ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_red, "Y red                       ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_ilu3, "Y ilu3                      ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_ilu3r, "Y ilu3 + red                ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_elu2r, "Y elu2 + red                ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_elu2rf, "Y elu2 + red + factor       ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_ilu3_elu2r, "Y ilu3 + elu2 + red         ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_ilu3_elu2rf, "Y ilu3 + elu2 + red + factor");
    putchar('\n'); /**/
    // puts("done\n");

    // free
    free_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);

    free_ui8matrix(Y_bas, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_reg, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_rot, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_red, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_ilu3, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_ilu3r, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_elu2r, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_elu2rf, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_ilu3_elu2r, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_ilu3_elu2rf, 0, h - 1, 0, w1 - 1);
}
// ---------------------------------------
void test_morpho_min_routine(int h, int w0)
// ---------------------------------------
{
    int b = 1; // bord

    int w8 = w0 / 8;
    if (w0 % 8)
        w8 = w8 + 1;
    int w1 = 8 * w8; // w1 >= w

    puts("--------------------------------------------------");
    printf("test_morpho_min_routine h = %d w0 = %d w8 = %d, w1 = %d\n", h, w0, w8, w1);
    if (w1 > w0)
        puts("w1 > w0");
    uint8 **X;
    uint8 **Y_bas, **Y_reg, **Y_rot, **Y_red, **Y_ilu3, **Y_ilu3r, **Y_elu2r, **Y_elu2rf, **Y_ilu3_elu2r, **Y_ilu3_elu2rf;

    int c; // error

    // puts("malloc");
    X = ui8matrix(0 - b, h - 1 + b, 0 - b, w1 - 1 + b);

    Y_bas = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_reg = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_rot = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_ilu3 = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_ilu3r = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_elu2r = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_elu2rf = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_ilu3_elu2r = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_ilu3_elu2rf = ui8matrix(0, h - 1, 0, w1 - 1);

    // puts("zero");
    zero_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);
    zero_ui8matrix(Y_bas, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_reg, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_rot, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_red, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_ilu3, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_ilu3r, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_elu2r, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_elu2rf, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_ilu3_elu2r, 0, h - 1, 0, w1 - 1);
    zero_ui8matrix(Y_ilu3_elu2rf, 0, h - 1, 0, w1 - 1);

    // puts("rand");
    // rand_ui8matrix(X, 0, h-1, 0, w0-1, 1, 255); // niveau de gris [1,255]
    rand1_ui8matrix(X, 0, h - 1, 0, w0 - 1, 80); // binaire [0,1]
    traitement_bordure(X, 0, h - 1,0, w0 - 1, b);

    // display_ui8matrix (X,  0, h-1, 0, w0-1, "%5d", "X0");
    // display_ui8matrix (X,  0 - b, h-1 + b, 0 - b, w0-1 + b, "%5d", "X0");

    // puts("\n-- max3 --");
    min3_ui8matrix_basic(X, 0, h - 1, 0, w0 - 1, Y_bas);                        // puts("bas");
    min3_ui8matrix_reg(X, 0, h - 1, 0, w0 - 1, Y_reg);                          // puts("reg");
    min3_ui8matrix_rot(X, 0, h - 1, 0, w0 - 1, Y_rot);                          // puts("rot");
    min3_ui8matrix_red(X, 0, h - 1, 0, w0 - 1, Y_red);                          // puts("red");
    min3_ui8matrix_ilu3(X, 0, h - 1, 0, w0 - 1, Y_ilu3);                        // puts("ilu3");
    min3_ui8matrix_ilu3_red(X, 0, h - 1, 0, w0 - 1, Y_ilu3r);                   // puts("ilu3_red");
    min3_ui8matrix_elu2_red(X, 0, h - 1, 0, w0 - 1, Y_elu2r);                   // puts("elu2_red");
    min3_ui8matrix_elu2_red_factor(X, 0, h - 1, 0, w0 - 1, Y_elu2rf);           // puts("elu2_red_factor");
    min3_ui8matrix_ilu3_elu2_red(X, 0, h - 1, 0, w0 - 1, Y_ilu3_elu2r);         // puts("ilu3_elu2_red");
    min3_ui8matrix_ilu3_elu2_red_factor(X, 0, h - 1, 0, w0 - 1, Y_ilu3_elu2rf); // puts("ilu3_elu2_red_factor");

    // puts("done\n");
    // display_ui8matrix(Y_bas,         0, h-1, 0, w0-1, "%5d", "Y bas        ");
    // display_ui8matrix(Y_reg,         0, h-1, 0, w0-1, "%5d", "Y reg        ");
    // display_ui8matrix(Y_rot,         0, h-1, 0, w0-1, "%5d", "Y rot        ");
    // display_ui8matrix(Y_red,         0, h-1, 0, w0-1, "%5d", "Y red        ");
    // display_ui8matrix(Y_ilu3,        0, h-1, 0, w0-1, "%5d", "Y ilu3       ");
    // display_ui8matrix(Y_ilu3r,       0, h-1, 0, w0-1, "%5d", "Y ilu3r      ");
    // display_ui8matrix(Y_elu2r,       0, h-1, 0, w0-1, "%5d", "Y elu2r      ");
    // display_ui8matrix(Y_elu2rf,      0, h-1, 0, w0-1, "%5d", "Y elu2rf     ");
    // display_ui8matrix(Y_ilu3_elu2r, 0, h-1, 0, w0-1, "%5d", "Y ilu3 elu2r");
    // display_ui8matrix(Y_ilu3_elu2rf, 0, h-1, 0, w0-1, "%5d", "Y ilu3 elu2rf");

    // puts("-- compare max --");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_reg, "Y reg                       ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_rot, "Y rot                       ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_red, "Y red                       ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_ilu3, "Y ilu3                      ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_ilu3r, "Y ilu3 + red                ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_elu2r, "Y elu2 + red                ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_elu2rf, "Y elu2 + red + factor       ");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_ilu3_elu2r, "Y ilu3 + elu2 + red");
    c = compare_ui8matrix(Y_bas, 0, h - 1, 0, w0 - 1, Y_ilu3_elu2rf, "Y ilu3 + elu2 + red + factor");
    putchar('\n'); /**/
    // puts("done\n");

    free_ui8matrix(X, 0 - b, h - 1 + b, 0 - b, w1 - 1 + b);

    free_ui8matrix(Y_bas, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_reg, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_rot, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_red, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_ilu3, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_ilu3r, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_elu2r, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_elu2rf, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_ilu3_elu2r, 0, h - 1, 0, w1 - 1);
    free_ui8matrix(Y_ilu3_elu2rf, 0, h - 1, 0, w1 - 1);
}
// ----------------------------------------------
void test_morpho_ouverture_routine(int h, int w0)
// ----------------------------------------------
{
    // avec des ADD
    int r = 1; // rayon des operateurs
    int w8 = w0 / 8;
    if (w0 % 8)
        w8 = w8 + 1;
    int w1 = 8 * w8; // w1 >= w

    puts("--------------------------------------------------");
    printf("test_morpho_ouverture_routine h = %d w0 = %d w8 = %d, w1 = %d\n", h, w0, w8, w1);
    if (w1 > w0)
        puts("w1 > w0");

    uint8 **X;
    uint8 **X8, **X1;

    uint8 **T_basic;

    uint8 **Y_basic;
    uint8 **Y_fusion;
    uint8 **Y_fusion_ilu5_red;
    uint8 **Y_fusion_ilu5_elu2_red;
    uint8 **Y_fusion_ilu5_elu2_red_factor;
    uint8 **Y_fusion_ilu15_red;

    uint8 **T_pipeline, **Y_pipeline;
    uint8 **T_pipeline_red, **Y_pipeline_red;
    uint8 **T_pipeline_ilu3_red, **Y_pipeline_ilu3_red;
    uint8 **T_pipeline_elu2_red, **Y_pipeline_elu2_red;
    uint8 **T_pipeline_elu2_red_factor, **Y_pipeline_elu2_red_factor;
    uint8 **T_pipeline_ilu3_elu2_red, **Y_pipeline_ilu3_elu2_red;
    uint8 **T_pipeline_ilu3_elu2_red_factor, **Y_pipeline_ilu3_elu2_red_factor;

    char *format1; // affichage 1 bit/pixel
    char *format8; // affichage 8 bits/pixel

    int c; // error

    format8 = "%3x";
    format1 = "%4d";
    format1 = "%2d";
    // format1 = "%1d";

    puts("malloc");
    // X 2r-border
    X = ui8matrix(0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);
    X8 = ui8matrix(0 - 2 * r, h - 1 + 2 * r, 0 - 1 * r, w8 - 1 + 1 * r);

    // T 1r-border
    T_basic = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_red = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_ilu3_red = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_elu2_red = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_elu2_red_factor = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_ilu3_elu2_red = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_ilu3_elu2_red_factor = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);

    // Y 0r-border
    Y_basic = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion_ilu5_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion_ilu5_elu2_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion_ilu5_elu2_red_factor = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion_ilu15_red = ui8matrix(0, h - 1, 0, w1 - 1);

    Y_pipeline = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_ilu3_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_elu2_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_elu2_red_factor = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_ilu3_elu2_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_ilu3_elu2_red_factor = ui8matrix(0, h - 1, 0, w1 - 1);

    // PUTS("zero X");
    zero_ui8matrix(X, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);

    // PUTS("zero T");
    zero_ui8matrix(T_basic, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_ilu3_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_elu2_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_elu2_red_factor, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_ilu3_elu2_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_ilu3_elu2_red_factor, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);

    // PUTS("zero Y");
    zero_ui8matrix(Y_basic, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion_ilu5_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion_ilu5_elu2_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion_ilu5_elu2_red_factor, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion_ilu15_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);

    zero_ui8matrix(Y_pipeline, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_ilu3_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_elu2_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_elu2_red_factor, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_ilu3_elu2_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_ilu3_elu2_red_factor, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);

    // rand1_ui8matrix(X, 0, h-1, 0, w0-1, 80);
    rand1_ui8matrix(X, 0, h - 1, 0, w0 - 1, 90); // binaire [0,1] (afin d'avoir des '1' apres traitement et pas que des '0');
    traitement_bordure(X, 0, h - 1,0, w0 - 1, r);
    // display_ui8matrix(X,                  0-2*r, h-1+2*r, 0-2*r, w0-1+2*r, format1, "X +2r");

    // --ouverture3_ui8matrix_basic = min puis max
    ouverture3_ui8matrix_basic(X, 0, h - 1, 0, w0 - 1, T_basic, Y_basic);
    ouverture3_ui8matrix_fusion(X, 0, h - 1, 0, w0 - 1, Y_fusion);
    ouverture3_ui8matrix_fusion_ilu5_red(X, 0, h - 1, 0, w0 - 1, Y_fusion_ilu5_red);
    ouverture3_ui8matrix_fusion_ilu5_elu2_red(X, 0, h - 1, 0, w0 - 1, Y_fusion_ilu5_elu2_red);
    ouverture3_ui8matrix_fusion_ilu5_elu2_red_factor(X, 0, h - 1, 0, w0 - 1, Y_fusion_ilu5_elu2_red_factor);
    ouverture3_ui8matrix_fusion_ilu15_red(X, 0, h - 1, 0, w0 - 1, Y_fusion_ilu15_red);

    // display_ui8matrix(Y_basic                      ,  0, h-1, 0, w0-1, format1, "Y_basic                      ");
    // display_ui8matrix(Y_fusion                     ,  0, h-1, 0, w0-1, format1, "Y_fusion                     ");
    // display_ui8matrix(Y_fusion_ilu5_red            ,  0, h-1, 0, w0-1, format1, "Y_fusion_ilu5_red            ");
    // display_ui8matrix(Y_fusion_ilu5_elu2_red       ,  0, h-1, 0, w0-1, format1, "Y_fusion_ilu5_elu2_red       ");
    // display_ui8matrix(Y_fusion_ilu5_elu2_red_factor,  0, h-1, 0, w0-1, format1, "Y_fusion_ilu5_elu2_red_factor");
    // display_ui8matrix(Y_fusion_ilu15_red           ,  0, h-1, 0, w0-1, format1, "Y_fusion_ilu15_red           ");

    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_fusion, "Y_fusion                     ");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_fusion_ilu5_red, "Y_fusion_ilu5_red            ");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_fusion_ilu5_elu2_red, "Y_fusion_ilu5_elu2_red       ");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_fusion_ilu5_elu2_red_factor, "Y_fusion_ilu5_elu2_red_factor");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_fusion_ilu15_red, "Y_fusion_ilu15_red           ");

    puts("--------------");
    puts("-- pipeline --");
    puts("--------------");

    // T: full-matrix
    ouverture3_ui8matrix_pipeline_basic(X, 0, h - 1, 0, w0 - 1, T_pipeline, Y_pipeline);
    ouverture3_ui8matrix_pipeline_red(X, 0, h - 1, 0, w0 - 1, T_pipeline_red, Y_pipeline_red);
    ouverture3_ui8matrix_pipeline_ilu3_red(X, 0, h - 1, 0, w0 - 1, T_pipeline_ilu3_red, Y_pipeline_ilu3_red);
    ouverture3_ui8matrix_pipeline_elu2_red(X, 0, h - 1, 0, w0 - 1, T_pipeline_elu2_red, Y_pipeline_elu2_red);
    ouverture3_ui8matrix_pipeline_elu2_red_factor(X, 0, h - 1, 0, w0 - 1, T_pipeline_elu2_red_factor, Y_pipeline_elu2_red_factor);
    ouverture3_ui8matrix_pipeline_ilu3_elu2_red(X, 0, h - 1, 0, w0 - 1, T_pipeline_ilu3_elu2_red, Y_pipeline_ilu3_elu2_red);
    ouverture3_ui8matrix_pipeline_ilu3_elu2_red_factor(X, 0, h - 1, 0, w0 - 1, T_pipeline_ilu3_elu2_red_factor, Y_pipeline_ilu3_elu2_red_factor);

    // display_ui8matrix(T_basic                        , 0-r, h-1+r, 0-r, w0-1+r, format1, "T_basic                        ");
    // display_ui8matrix(T_pipeline                     , 0-r, h-1+r, 0-r, w0-1+r, format1, "T_pipeline                     ");
    // display_ui8matrix(T_pipeline_red                 , 0-r, h-1+r, 0-r, w0-1+r, format1, "T_pipeline_red                 ");
    // display_ui8matrix(T_pipeline_ilu3_red            , 0-r, h-1+r, 0-r, w0-1+r, format1, "T_pipeline_ilu3_red            ");
    // display_ui8matrix(T_pipeline_elu2_red            , 0-r, h-1+r, 0-r, w0-1+r, format1, "T_pipeline_elu2_red            ");
    // display_ui8matrix(T_pipeline_elu2_red_factor     , 0-r, h-1+r, 0-r, w0-1+r, format1, "T_pipeline_elu2_red_factor     ");
    // display_ui8matrix(T_pipeline_ilu3_elu2_red       , 0-r, h-1+r, 0-r, w0-1+r, format1, "T_pipeline_ilu3_elu2_red       ");
    // display_ui8matrix(T_pipeline_ilu3_elu2_red_factor, 0-r, h-1+r, 0-r, w0-1+r, format1, "T_pipeline_ilu3_elu2_red_factor");

    c = compare_ui8matrix(T_basic, 0 - r, h - 1 + r, 0 - r, w0 - 1 + r, T_pipeline, "T_pipeline                     ");
    c = compare_ui8matrix(T_basic, 0 - r, h - 1 + r, 0 - r, w0 - 1 + r, T_pipeline_red, "T_pipeline_red                 ");
    c = compare_ui8matrix(T_basic, 0 - r, h - 1 + r, 0 - r, w0 - 1 + r, T_pipeline_ilu3_red, "T_pipeline_ilu3_red            ");
    c = compare_ui8matrix(T_basic, 0 - r, h - 1 + r, 0 - r, w0 - 1 + r, T_pipeline_elu2_red, "T_pipeline_elu2_red            ");
    c = compare_ui8matrix(T_basic, 0 - r, h - 1 + r, 0 - r, w0 - 1 + r, T_pipeline_elu2_red_factor, "T_pipeline_elu2_red_factor     ");
    c = compare_ui8matrix(T_basic, 0 - r, h - 1 + r, 0 - r, w0 - 1 + r, T_pipeline_ilu3_elu2_red, "T_pipeline_ilu3_elu2_red       ");
    c = compare_ui8matrix(T_basic, 0 - r, h - 1 + r, 0 - r, w0 - 1 + r, T_pipeline_ilu3_elu2_red_factor, "T_pipeline_ilu3_elu2_red_factor");
    putchar('\n');

    // display_ui8matrix(Y_basic                        , 0, h-1, 0, w0-1, format1, "Y_basic                        ");
    // display_ui8matrix(Y_pipeline                     , 0, h-1, 0, w0-1, format1, "Y_pipeline                     ");
    // display_ui8matrix(Y_pipeline_red                 , 0, h-1, 0, w0-1, format1, "Y_pipeline_red                 ");
    // display_ui8matrix(Y_pipeline_ilu3_red            , 0, h-1, 0, w0-1, format1, "Y_pipeline_ilu3_red            ");
    // display_ui8matrix(Y_pipeline_elu2_red            , 0, h-1, 0, w0-1, format1, "Y_pipeline_elu2_red            ");
    // display_ui8matrix(Y_pipeline_elu2_red_factor     , 0, h-1, 0, w0-1, format1, "Y_pipeline_elu2_red_factor     ");
    // display_ui8matrix(Y_pipeline_ilu3_elu2_red       , 0, h-1, 0, w0-1, format1, "Y_pipeline_ilu3_elu2_red       ");
    // display_ui8matrix(Y_pipeline_ilu3_elu2_red_factor, 0, h-1, 0, w0-1, format1, "Y_pipeline_ilu3_elu2_red_factor");

    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_pipeline, "Y_pipeline                     ");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_pipeline_red, "Y_pipeline_red                 ");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_pipeline_ilu3_red, "Y_pipeline_ilu3_red            ");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_pipeline_elu2_red, "Y_pipeline_elu2_red            ");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_pipeline_elu2_red_factor, "Y_pipeline_elu2_red_factor     ");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_pipeline_ilu3_elu2_red, "Y_pipeline_ilu3_elu2_red       ");
    c = compare_ui8matrix(Y_basic, 0, h - 1, 0, w0 - 1, Y_pipeline_ilu3_elu2_red_factor, "Y_pipeline_ilu3_elu2_red_factor");
    putchar('\n');

    puts("----------");
    puts("-- free --");
    puts("----------");

    free_ui8matrix(X, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);
    free_ui8matrix(X8, 0 - 2 * r, h - 1 + 2 * r, 0 - 1 * r, w8 - 1 + 1 * r);
    free_ui8matrix(X1, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);

    free_ui8matrix(T_basic, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    free_ui8matrix(T_pipeline, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    free_ui8matrix(T_pipeline_ilu3_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    free_ui8matrix(T_pipeline_elu2_red_factor, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);

    free_ui8matrix(Y_basic, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    free_ui8matrix(Y_fusion, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    free_ui8matrix(Y_fusion_ilu5_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    free_ui8matrix(Y_fusion_ilu5_elu2_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    free_ui8matrix(Y_fusion_ilu5_elu2_red_factor, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    free_ui8matrix(Y_fusion_ilu15_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);

    free_ui8matrix(T_pipeline, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    // puts("free done");
}
// -----------------------
void test_morpho_max(void)
// -----------------------
{
    int h0 = 8;
    int w0 = 16;

    int dh = 4;
    int dw = 3;

    for (int h = h0; h <= h0 + dh; h++)
    { // pour tester elu2
        for (int w = w0; w <= w0 + dw; w++)
        { // pour tester ilu3
            test_morpho_max_routine(h, w);
        }
    }
}
// -----------------------
void test_morpho_min(void)
// -----------------------
{
    int h0 = 8;
    int w0 = 16;

    int dh = 4;
    int dw = 3;

    for (int h = h0; h <= h0 + dh; h++)
    { // pour tester elu2
        for (int w = w0; w <= w0 + dw; w++)
        { // pour tester ilu3
            test_morpho_min_routine(h, w);
        }
    }
}
// -----------------------------
void test_morpho_ouverture(void)
// -----------------------------
{
    int h0 = 4;
    int w0 = 16;

    int dh = 4;
    int dw = 3;

    // dh = dw = 0;

    dh = 5;
    dw = 4;
    h0 = 8;
    w0 = 32;

    for (int h = h0; h <= h0 + dh; h++)
    { // pour tester elu2
        for (int w = w0; w <= w0 + dw; w++)
        { // pour tester ilu3
            test_morpho_ouverture_routine(h, w);
        }
    }
}
// ----------------------
void test_wikipedia(void)
// ----------------------
{
    //
    int b = 1; // bord rayon de l'element structurant
    int h = 17;
    int w = 20 + 4;

    // images binaires
    uint8 **X;
    uint8 **E; // erosion
    uint8 **D; // dilatation
    uint8 **O; // ouverture
    uint8 **F; // fermeture

    // images packees
    uint8 **X8;
    uint8 **E8; // erosion
    uint8 **D8; // dilatation
    uint8 **O8; // ouverture
    uint8 **F8; // fermeture

    // images unpackees
    uint8 **E1; // erosion
    uint8 **D1; // dilatation
    uint8 **O1; // ouverture
    uint8 **F1; // fermeture

    X = ui8matrix(0 - 2 * b, h - 1 + 2 * b, 0 - 2 * b, w - 1 + 2 * b);
    E = ui8matrix(0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    D = ui8matrix(0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    O = ui8matrix(0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);
    F = ui8matrix(0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);

    zero_ui8matrix(X, 0 - 2 * b, h - 1 + 2 * b, 0 - 2 * b, w - 1 + 2 * b);
    zero_ui8matrix(E, 0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    zero_ui8matrix(D, 0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    zero_ui8matrix(O, 0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);
    zero_ui8matrix(F, 0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);

    X8 = ui8matrix(0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w / 8 - 1 + 1 * b);
    E8 = ui8matrix(0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w / 8 - 1 + 1 * b);
    D8 = ui8matrix(0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w / 8 - 1 + 1 * b);
    O8 = ui8matrix(0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w / 8 - 1 + 0 * b);
    F8 = ui8matrix(0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w / 8 - 1 + 0 * b);

    zero_ui8matrix(X8, 0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w / 8 - 1 + 1 * b);
    zero_ui8matrix(E8, 0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w / 8 - 1 + 1 * b);
    zero_ui8matrix(D8, 0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w / 8 - 1 + 1 * b);
    zero_ui8matrix(O8, 0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w / 8 - 1 + 0 * b);
    zero_ui8matrix(F8, 0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w / 8 - 1 + 0 * b);

    E1 = ui8matrix(0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    D1 = ui8matrix(0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    O1 = ui8matrix(0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);
    F1 = ui8matrix(0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);

    zero_ui8matrix(E1, 0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    zero_ui8matrix(D1, 0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    zero_ui8matrix(O1, 0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);
    zero_ui8matrix(F1, 0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);

    int i = 0;

    // il faut un bord de 1 dans l'image pour que la modification (dilatation) soit visible *dans* l'image

    set_ui8vector_str(X[i++], 0, w - 1, ".....................");
    set_ui8vector_str(X[i++], 0, w - 1, "....11.........11....");
    set_ui8vector_str(X[i++], 0, w - 1, ".111..1..111..1..111.");
    set_ui8vector_str(X[i++], 0, w - 1, ".111...1111111...111.");
    set_ui8vector_str(X[i++], 0, w - 1, ".111..111111111..111.");
    set_ui8vector_str(X[i++], 0, w - 1, ".....11111111111.....");
    set_ui8vector_str(X[i++], 0, w - 1, ".....11111111111.....");
    set_ui8vector_str(X[i++], 0, w - 1, "....1111.111.1111....");
    set_ui8vector_str(X[i++], 0, w - 1, "....1111111111111....");
    set_ui8vector_str(X[i++], 0, w - 1, "....111111.111111....");
    set_ui8vector_str(X[i++], 0, w - 1, "....1111111111111....");
    set_ui8vector_str(X[i++], 0, w - 1, ".....1111...1111.....");
    set_ui8vector_str(X[i++], 0, w - 1, ".....1111...1111.....");
    set_ui8vector_str(X[i++], 0, w - 1, "......111...111......");
    set_ui8vector_str(X[i++], 0, w - 1, ".......1111111.......");
    set_ui8vector_str(X[i++], 0, w - 1, ".........111.........");
    set_ui8vector_str(X[i++], 0, w - 1, ".....................");

    // display_ui8matrix_positive(X, 0, h - 1, 0, w - 1, 3, "X");

    // ------------------------
    // operateurs 1 bit / pixel
    // ------------------------
    min3_ui8matrix_basic(X, 0, h - 1, 0, w - 1, E);
    max3_ui8matrix_basic(X, 0, h - 1, 0, w - 1, D);

    // max3_ui8matrix_basic(E, 0, h-1, 0, w-1, O);
    // min3_ui8matrix_basic(D, 0, h-1, 0, w-1, F);

    // display_ui8matrix_positive(E, 0, h-1, 0, w-1, 3, "Erosion");
    // display_ui8matrix_positive(D, 0, h-1, 0, w-1, 3, "Dilatation");

    // display_ui8matrix_positive(O, 0, h-1, 0, w-1, 3, "Ouverture");
    // display_ui8matrix_positive(F, 0, h-1, 0, w-1, 3, "Fermeture");

    free_ui8matrix(X, 0 - 2 * b, h - 1 + 2 * b, 0 - 2 * b, w - 1 + 2 * b);
    free_ui8matrix(E, 0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    free_ui8matrix(D, 0 - 1 * b, h - 1 + 1 * b, 0 - 1 * b, w - 1 + 1 * b);
    free_ui8matrix(O, 0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);
    free_ui8matrix(F, 0 - 0 * b, h - 1 + 0 * b, 0 - 0 * b, w - 1 + 0 * b);
}
// ---------------------------------------------------
void bench_morpho_ouverture(int n0, int n1, int nstep)
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

    uint8 **X;
    uint8 **X8, **X1;

    uint8 **T_basic;

    uint8 **Y_basic;
    uint8 **Y_fusion;
    uint8 **Y_fusion_ilu5_red;
    uint8 **Y_fusion_ilu5_elu2_red;
    uint8 **Y_fusion_ilu5_elu2_red_factor;
    uint8 **Y_fusion_ilu15_red;

    uint8 **T_pipeline, **Y_pipeline;
    uint8 **T_pipeline_red, **Y_pipeline_red;
    uint8 **T_pipeline_ilu3_red, **Y_pipeline_ilu3_red;
    uint8 **T_pipeline_elu2_red, **Y_pipeline_elu2_red;
    uint8 **T_pipeline_elu2_red_factor, **Y_pipeline_elu2_red_factor;
    uint8 **T_pipeline_ilu3_elu2_red, **Y_pipeline_ilu3_elu2_red;
    uint8 **T_pipeline_ilu3_elu2_red_factor, **Y_pipeline_ilu3_elu2_red_factor;

    double cpp_basic;

    double cpp_fusion;
    double cpp_fusion_red;
    double cpp_fusion_ilu5_red;
    double cpp_fusion_elu2_red;
    double cpp_fusion_elu2_red_factor;
    double cpp_fusion_ilu5_elu2_red;
    double cpp_fusion_ilu5_elu2_red_factor;
    double cpp_fusion_ilu15_red;

    double cpp_pipeline_basic;
    double cpp_pipeline_red;
    double cpp_pipeline_ilu3_red;
    double cpp_pipeline_elu2_red;
    double cpp_pipeline_elu2_red_factor;
    double cpp_pipeline_ilu3_elu2_red;
    double cpp_pipeline_ilu3_elu2_red_factor;

    char *format = "%8.2f";

    format = "%5.0f";
    format = "%6.1f";

    puts("=============================");
    puts("== bench_morpho_ouverture ===");
    puts("=============================");

    puts("malloc");

    // X 2r-border
    X = ui8matrix(0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);

    // T 1r-border
    T_basic = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_red = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_ilu3_red = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_elu2_red = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_elu2_red_factor = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_ilu3_elu2_red = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    T_pipeline_ilu3_elu2_red_factor = ui8matrix(0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);

    // Y 0r-border
    Y_basic = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion_ilu5_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion_ilu5_elu2_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion_ilu5_elu2_red_factor = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_fusion_ilu15_red = ui8matrix(0, h - 1, 0, w1 - 1);

    Y_pipeline = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_ilu3_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_elu2_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_elu2_red_factor = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_ilu3_elu2_red = ui8matrix(0, h - 1, 0, w1 - 1);
    Y_pipeline_ilu3_elu2_red_factor = ui8matrix(0, h - 1, 0, w1 - 1);

    // PUTS("zero X");
    zero_ui8matrix(X, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);

    // PUTS("zero T");
    zero_ui8matrix(T_basic, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_ilu3_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_elu2_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_elu2_red_factor, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_ilu3_elu2_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
    zero_ui8matrix(T_pipeline_ilu3_elu2_red_factor, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);

    // PUTS("zero Y");
    zero_ui8matrix(Y_basic, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion_ilu5_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion_ilu5_elu2_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion_ilu5_elu2_red_factor, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_fusion_ilu15_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);

    zero_ui8matrix(Y_pipeline, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_ilu3_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_elu2_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_elu2_red_factor, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_ilu3_elu2_red, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);
    zero_ui8matrix(Y_pipeline_ilu3_elu2_red_factor, 0 - 0 * r, h - 1 + 0 * r, 0 - 0 * r, w1 - 1 + 0 * r);

    puts("temps de calcul en ccp (cycle/point)");

    for (int n = n0; n <= n1; n += nstep)
    {

        h = n;
        w8 = n / 8;
        w1 = n / 1;

        // printf("i = %3d\n", n);

        resize_ui8matrix(X, 0 - 2 * r, h - 1 + 2 * r, 0 - 2 * r, w1 - 1 + 2 * r);

        // T 1r-border
        resize_ui8matrix(T_basic, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
        resize_ui8matrix(T_pipeline, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
        resize_ui8matrix(T_pipeline_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
        resize_ui8matrix(T_pipeline_ilu3_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
        resize_ui8matrix(T_pipeline_elu2_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
        resize_ui8matrix(T_pipeline_elu2_red_factor, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
        resize_ui8matrix(T_pipeline_ilu3_elu2_red, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);
        resize_ui8matrix(T_pipeline_ilu3_elu2_red_factor, 0 - 1 * r, h - 1 + 1 * r, 0 - 1 * r, w1 - 1 + 1 * r);

        // Y 0r-border
        resize_ui8matrix(Y_basic, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_fusion, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_fusion_ilu5_red, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_fusion_ilu5_elu2_red, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_fusion_ilu5_elu2_red_factor, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_fusion_ilu15_red, 0, h - 1, 0, w1 - 1);

        resize_ui8matrix(Y_pipeline, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_pipeline_red, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_pipeline_ilu3_red, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_pipeline_elu2_red, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_pipeline_elu2_red_factor, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_pipeline_ilu3_elu2_red, 0, h - 1, 0, w1 - 1);
        resize_ui8matrix(Y_pipeline_ilu3_elu2_red_factor, 0, h - 1, 0, w1 - 1);

        BENCH(ouverture3_ui8matrix_basic(X, 0, h - 1, 0, w1 - 1, T_basic, Y_basic), n, cpp_basic);
        BENCH(ouverture3_ui8matrix_fusion(X, 0, h - 1, 0, w1 - 1, Y_fusion), n, cpp_fusion);
        BENCH(ouverture3_ui8matrix_fusion_ilu5_red(X, 0, h - 1, 0, w1 - 1, Y_fusion_ilu5_red), n, cpp_fusion_ilu5_red);
        BENCH(ouverture3_ui8matrix_fusion_ilu5_elu2_red(X, 0, h - 1, 0, w1 - 1, Y_fusion_ilu5_elu2_red), n, cpp_fusion_ilu5_elu2_red);
        BENCH(ouverture3_ui8matrix_fusion_ilu5_elu2_red_factor(X, 0, h - 1, 0, w1 - 1, Y_fusion_ilu5_elu2_red_factor), n, cpp_fusion_ilu5_elu2_red_factor);
        BENCH(ouverture3_ui8matrix_fusion_ilu15_red(X, 0, h - 1, 0, w1 - 1, Y_fusion_ilu15_red), n, cpp_fusion_ilu15_red);

        BENCH(ouverture3_ui8matrix_pipeline_basic(X, 0, h - 1, 0, w1 - 1, T_pipeline, Y_pipeline), n, cpp_pipeline_basic);
        BENCH(ouverture3_ui8matrix_pipeline_red(X, 0, h - 1, 0, w1 - 1, T_pipeline_red, Y_pipeline_red), n, cpp_pipeline_red);
        BENCH(ouverture3_ui8matrix_pipeline_ilu3_red(X, 0, h - 1, 0, w1 - 1, T_pipeline_ilu3_red, Y_pipeline_ilu3_red), n, cpp_pipeline_ilu3_red);
        BENCH(ouverture3_ui8matrix_pipeline_elu2_red(X, 0, h - 1, 0, w1 - 1, T_pipeline_elu2_red, Y_pipeline_elu2_red), n, cpp_pipeline_elu2_red);
        BENCH(ouverture3_ui8matrix_pipeline_elu2_red_factor(X, 0, h - 1, 0, w1 - 1, T_pipeline_elu2_red_factor, Y_pipeline_elu2_red_factor), n, cpp_pipeline_elu2_red_factor);
        BENCH(ouverture3_ui8matrix_pipeline_ilu3_elu2_red(X, 0, h - 1, 0, w1 - 1, T_pipeline_ilu3_elu2_red, Y_pipeline_ilu3_elu2_red), n, cpp_pipeline_ilu3_elu2_red);
        BENCH(ouverture3_ui8matrix_pipeline_ilu3_elu2_red_factor(X, 0, h - 1, 0, w1 - 1, T_pipeline_ilu3_elu2_red_factor, Y_pipeline_ilu3_elu2_red_factor), n, cpp_pipeline_ilu3_elu2_red_factor);

        /**/
        printf("i = %4d", n);
        printf("   ");
        printf(format, cpp_basic);
        printf(format, cpp_fusion);
        printf(format, cpp_fusion_ilu5_red);
        printf(format, cpp_fusion_ilu5_elu2_red);
        printf(format, cpp_fusion_ilu5_elu2_red_factor);
        printf(format, cpp_fusion_ilu15_red);

        printf("   ");

        printf(format, cpp_pipeline_basic);
        printf(format, cpp_pipeline_red);
        printf(format, cpp_pipeline_ilu3_red);
        printf(format, cpp_pipeline_elu2_red);
        printf(format, cpp_pipeline_elu2_red_factor);
        printf(format, cpp_pipeline_ilu3_elu2_red);
        printf(format, cpp_pipeline_ilu3_elu2_red_factor);

        putchar('\n');
    }
}
// ====================================
int test_morpho(int argc, char *argv[])
// ====================================
{
    puts("=== test_morpho ===");

    test_morpho_min();
    test_morpho_max();

    test_morpho_ouverture();

    // bench_morpho_ouverture(128, 512, 8);
    // bench_morpho_ouverture(128, 1024, 8);

    // test_set_str();
    // test_wikipedia();

    return 0;
}
