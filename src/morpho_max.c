/* -------------------- */
/* --- morpho_max.c --- */
/* -------------------- */

/*
 * Copyright (c) 2004 - 2013, Lionel Lacassagne, All rights reserved
 * University of Paris Sud, Laboratoire de Recherche en Informatique 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"


// ------------------------------------------------------------------------
void line_max3_ui8matrix_basic(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ------------------------------------------------------------------------
{
    uint8 max;

    for (int j = j0 ; j <= j1 ; j += 1) {
        
        max = max3( 
                    max3(load2(X,i-1,j-1) , load2(X,i-1,j) , load2(X,i-1,j+1)),
                    max3(load2(X,i,j-1) , load2(X, i, j), load2(X,i,j+1)), 
                    max3(load2(X,i+1, j-1) , load2(X,i+1,j) , load2(X,i+1,j+1))
                  );
        store2(Y,i,j,max);

    }
}
// ----------------------------------------------------------------------
void line_max3_ui8matrix_reg(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------
{
    uint8 max;
    uint8 b0,b1,b2,b3,b4,b5,b6,b7,b8;

    for (int j = j0 ; j <= j1    ; j += 1) {
        b0 = load2(X,i,j);
        b1 = load2(X,i-1,j-1)   ; b2 = load2(X,i-1,j)   ; b3 = load2(X,i-1,j+1);
        b4 = load2(X,i,j-1)     ; b5 = load2(X,i,j+1)   ;
        b6 = load2(X,i+1, j-1)  ; b7 = load2(X,i+1,j)   ; b8 = load2(X,i+1,j+1);
        
        max = max3( 
                    max3(b1,b2,b3),
                    max3(b4,b0,b5), 
                    max3(b6,b7,b8)
                  );

        store2(Y,i,j, max);

    }
}
// ----------------------------------------------------------------------
void line_max3_ui8matrix_rot(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------
{

    uint8 max;
    uint8 b0,b1,b2,b3,b4,b5,b6,b7,b8;

    b0 = load2(X,i,j0)       ; b1 = load2(X,i-1,j0-1)   ; 
    b2 = load2(X,i-1,j0)     ; b4 = load2(X,i,j0-1)     ; 
    b6 = load2(X,i+1, j0-1)  ; b7 = load2(X,i+1,j0)   ; 

    for (int j = j0 ; j <= j1 ; j += 1) {
        b3 = load2(X,i-1,j+1)  ;
        b5 = load2(X,i,j+1)    ;
        b8 = load2(X,i+1,j+1)  ;

        max = max3( 
                    max3(b1,b2,b3),
                    max3(b4, b0, b5), 
                    max3(b6,b7,b8)
                  );

        store2(Y,i,j,max);

        //ROTATION
        b1 = b2; b4 = b0; b6 = b7;
        b2 = b3; b0 = b5; b7 = b8;

    }
}
// ----------------------------------------------------------------------
void line_max3_ui8matrix_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------
{
    uint8 max_1, max_2, max_3;

    max_1 = max3(load2(X, i-1, j0-1), load2(X, i, j0-1), load2(X, i+1, j0-1));
    max_2 = max3(load2(X, i-1, j0), load2(X, i, j0), load2(X, i+1, j0));

    for (int j = j0 ; j<=j1 ; j++) {
        max_3 = max3(load2(X, i-1, j+1), load2(X, i, j+1), load2(X, i+1, j+1));
        store2(Y, i, j, max3(max_1, max_2, max_3));
        max_1 = max_2; max_2 = max_3;
    }
}
// -----------------------------------------------------------------------
void line_max3_ui8matrix_ilu3(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------
// b1 b2 b3  b4   b5
// b6 X1 X2  X3   b7
// b8 b9 b10 b11  b12
{
    uint8 b1,b6,b8,b2,x1,b9,b3,x2,b10,b4,x3,b11;
    uint8 max_1, max_2, max_3;
    int r = (j1 - j0 + 1) % 3;


    b1 = load2(X, i-1, j0-1); b6 = load2(X, i, j0-1); b8 = load2(X, i+1, j0-1); 
    b2 = load2(X, i-1, j0); x1 = load2(X, i, j0); b9 = load2(X, i+1, j0);
    max_1 = max3(b1,b6,b8); max_2 = max3(b2,x1,b9);

     for (int j = j0 ; j <= j1-r ; j+=3) {

        b3 = load2(X, i-1, j+1); x2 = load2(X, i, j+1); b10 = load2(X, i+1, j+1);
        max_3 = max3(b3, x2, b10);
        store2(Y, i, j, max3(max_1, max_2, max_3));
        
        b1 = load2(X, i-1, j+2); b6 = load2(X, i, j+2); b8 = load2(X, i+1, j+2);
        max_1 = max3(b1, b6, b8);
        store2(Y, i, j+1, max3(max_2, max_3, max_1));

        b2 = load2(X, i-1, j+3);  x1 = load2(X, i, j+3);  b9 = load2(X, i+1, j+3);
        max_2 = max3(b2,x1,b9);
        store2(Y, i, j+2, max3(max_3, max_1, max_2));

    }

     if (r) {
        line_max3_ui8matrix_red(X,i,j1-r+1, j1, Y);
     }
}
// ---------------------------------------------------------------------------
void line_max3_ui8matrix_ilu3_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------
// b1 b2 b3  b4   b5
// b6 X1 X2  X3   b7
// b8 b9 b10 b11  b12
// ---------------------------------------------------------------------------
{

    uint8 max_1,max_2,max_3;

    max_1 = max3(load2(X, i-1, j0-1), load2(X, i, j0-1), load2(X,i+1, j0-1));
    max_2 = max3(load2(X,i-1,j0), load2(X,i,j0), load2(X,i+1,j0));
   
    //Prologue
    int r = (j1 - j0 + 1) % 3;

    for (int j = j0 ; j <= j1-r ; j+=3) {
        
        max_3 = max3(load2(X,i-1, j+1), load2(X,i,j+1), load2(X,i+1,j+1));
        store2(Y,i,j, max3(max_1, max_2, max_3));

        max_1 = max3(load2(X,i-1,j+2), load2(X,i,j+2), load2(X,i+1, j+2));
        store2(Y,i,j+1, max3(max_2, max_3, max_1));

        max_2 = max3(load2(X,i-1, j+3), load2(X, i, j+3), load2(X, i+1, j+3));
        store2(Y,i,j+2,max3(max_3, max_1, max_2));
    }

    // Epilogue
    if (r) { 
            line_max3_ui8matrix_red(X,i,j1-r+1, j1, Y);
    }   

}
// ---------------------------------------------------------------------------
void line_max3_ui8matrix_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------
// b1 b2 b3
// b4 X1 b5
// b6 X2 b7
// b8 b9 b10
{
uint8 b4,b5,b6,b7,x1,x2;
uint8 max_11, max_12, max_13, max_21, max_22, max_23;

    b4 = load2(X, i, j0-1); b6 = load2(X, i+1, j0-1); 
    x1 = load2(X, i, j0); x2 = load2(X, i+1, j0);

    max_11 = max3(load2(X, i-1, j0-1), b4, b6); 
    max_12 = max3(load2(X, i-1, j0), x1, x2); 
    max_21 = max3(b4, b6, load2(X, i+2, j0-1)); 
    max_22 = max3(x1, x2, load2(X, i+2, j0)); 

    for (int j = j0 ; j <= j1 ; j++) {

        b5 = load2(X, i, j+1); b7 = load2(X, i+1, j+1);
        max_13 = max3(load2(X, i-1, j+1), b5, b7);
        store2(Y, i, j, max3(max_11, max_12, max_13));
        max_11 = max_12; max_12 = max_13;

        max_23 = max3(b5, b7, load2(X, i+2, j+1));
        store2(Y, i+1, j, max3(max_21, max_22, max_23));
        max_21 = max_22; max_22 = max_23;

    }
}


// ----------------------------------------------------------------------------------
void line_max3_ui8matrix_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------
{

// b1 b2 b3
// b4 X1 b5
// b6 X2 b7
// b8 b9 b10

uint8 max_11, max_12, max_13, max_21, max_22, max_23;
uint8 fact;

    fact = max2(load2(X, i, j0-1),load2(X, i+1, j0-1));
    max_11 = max2(load2(X, i-1, j0-1), fact); 
    max_21 = max2(fact, load2(X, i+2, j0-1)); 

    fact = max2(load2(X, i, j0),load2(X, i+1, j0));
    max_12 = max2(load2(X, i-1, j0), fact); 
    max_22 = max2(fact , load2(X, i+2, j0)); 

    for (int j = j0 ; j <= j1 ; j++) {

        fact = max2(load2(X, i, j+1),load2(X, i+1, j+1));
        max_13 = max2(load2(X, i-1, j+1), fact);
        max_23 = max2(fact, load2(X, i+2, j+1));
        fact = max2(max_11, max_12);
        store2(Y, i, j, max2(fact, max_13));
        max_11 = max_12; max_12 = max_13;
        fact = max2(max_21, max_22);
        store2(Y, i+1, j, max2(fact, max_23));
        max_21 = max_22; max_22 = max_23;    
    }
}
// --------------------------------------------------------------------------------
void line_max3_ui8matrix_ilu3_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------
//b1  b2  b3  b4  b5
//b6  x1  x2  x3  b7
//b8  x4  x5  x6  b9
//b10 b11 b12 b13 b14
{
    uint8 b6, x1, x2, b8, x4, x5;
    uint8 c1, c2, c3, c4, c5, c6;

    b6 = load2(X, i, j0 - 1);
    b8 = load2(X, i + 1, j0 - 1);

    x1 = load2(X, i, j0);
    x4 = load2(X, i + 1, j0);


    c1 = max3(load2(X, i - 1, j0 - 1), b6, b8);
    c2 = max3(load2(X, i - 1, j0), x1, x4);

    c4 = max3(load2(X, i + 2, j0 - 1), b6, b8);
    c5 = max3(load2(X, i + 2, j0), x1, x4);


    int r = (j1 - j0 + 1) % 3;

    for (int j = j0; j <= j1 - r; j += 3)
    {
        x2 = load2(X, i, j + 1); 
        x5 = load2(X, i + 1, j + 1);

        c3 = max3(load2(X, i - 1, j + 1), x2, x5);
        c6 = max3(load2(X, i + 2, j + 1), x2, x5);

        store2(Y, i, j, max3(c1, c2, c3)); //write y1
        store2(Y, i + 1, j, max3(c4, c5, c6)); //write y4

        b6 = load2(X, i, j + 2); //x3
        b8 = load2(X, i + 1, j + 2); //x6
        c1 = max3(load2(X, i - 1, j + 2), b6, b8);
        c4 = max3(load2(X, i + 2, j + 2), b6, b8); 

        store2(Y, i, j + 1, max3(c2, c3, c1)); //write y2
        store2(Y, i + 1, j + 1, max3(c5, c6, c4)); //write y5

        x1 = load2(X, i, j + 3); //b7
        x4 = load2(X, i + 1, j + 3); //b9
        c2 = max3(load2(X ,i - 1, j + 3), x1, x4);
        c5 = max3(load2(X, i + 2, j + 3), x1, x4);

        store2(Y, i, j + 2, max3(c3, c1, c2)); //write y3
        store2(Y, i + 1, j + 2, max3(c6, c4, c5)); //write y6
    }

    if (r)
    {
        line_max3_ui8matrix_elu2_red_factor(X, i, j1 - r + 1, j1, Y);
    }
}

// ---------------------------------------------------------------------------------------
void line_max3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------
//b1  b2  b3  b4  b5
//b6  x1  x2  x3  b7
//b8  x4  x5  x6  b9
//b10 b11 b12 b13 b14
{
    uint8 b6, x1, x2, b8, x4, x5;
    uint8 c1, c2, c3, c4, c5, c6;
    uint8 tmp;

    b6 = load2(X, i, j0 - 1);
    b8 = load2(X, i + 1, j0 - 1);

    x1 = load2(X, i, j0);
    x4 = load2(X, i + 1, j0);

    tmp = max2(b6, b8);
    c1 = max2(load2(X, i - 1, j0 - 1), tmp);
    c4 = max2(load2(X, i + 2, j0 - 1), tmp);

    tmp = max2(x1, x4); 
    c2 = max2(load2(X, i - 1, j0    ), tmp);
    c5 = max2(load2(X, i + 2, j0    ), tmp);

    int r = (j1 - j0 + 1) % 3;

    for (int j = j0; j <= j1 - r; j += 3)
    {
        x2 = load2(X, i, j + 1); 
        x5 = load2(X, i + 1, j + 1);

        tmp = max2(x2, x5);

        c3 = max2(load2(X, i - 1, j + 1), tmp);
        c6 = max2(load2(X, i + 2, j + 1), tmp);

        store2(Y, i, j, max3(c1, c2, c3)); //write y1
        store2(Y, i + 1, j, max3(c4, c5, c6)); //write y4

        b6 = load2(X, i, j + 2); //x3
        b8 = load2(X, i + 1, j + 2); //x6

        tmp = max2(b6, b8);

        c1 = max2(load2(X, i - 1, j + 2), tmp);
        c4 = max2(load2(X, i + 2, j + 2), tmp); 

        store2(Y, i, j + 1, max3(c2, c3, c1)); //write y2
        store2(Y, i + 1, j + 1, max3(c5, c6, c4)); //write y5

        x1 = load2(X, i, j + 3); //b7
        x4 = load2(X, i + 1, j + 3); //b9

        tmp = max2(x1, x4);

        c2 = max2(load2(X ,i - 1, j + 3), tmp);
        c5 = max2(load2(X, i + 2, j + 3), tmp);

        store2(Y, i, j + 2, max3(c3, c1, c2)); //write y3
        store2(Y, i + 1, j + 2, max3(c6, c4, c5)); //write y6
    }

    if (r)
    {
        line_max3_ui8matrix_elu2_red_factor(X, i, j1 - r + 1, j1, Y);
    }
}
// ----------------------------------------------------------------------------
void max3_ui8matrix_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------
{
    for (int i = i0 ; i <= i1 ; i++) {
        line_max3_ui8matrix_basic(X,i,j0,j1,Y);
    }
}
// --------------------------------------------------------------------------
void max3_ui8matrix_reg(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------
{
    for (int i = i0 ; i <= i1 ; i++) {
        line_max3_ui8matrix_reg(X,i,j0,j1,Y);
    }
}
// --------------------------------------------------------------------------
void max3_ui8matrix_rot(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------
{
    
    for (int i = i0 ; i <= i1 ; i++) {
        line_max3_ui8matrix_rot(X,i,j0,j1,Y);
    }
}
// --------------------------------------------------------------------------
void max3_ui8matrix_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------
{
    
    for (int i = i0 ; i <= i1 ; i++) {
        line_max3_ui8matrix_red(X,i,j0,j1,Y);
    }
}
// ---------------------------------------------------------------------------
void max3_ui8matrix_ilu3(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------
{
    
    for (int i = i0 ; i <= i1 ; i++) {
        line_max3_ui8matrix_ilu3(X,i,j0,j1,Y);
    }
}
// -------------------------------------------------------------------------------
void max3_ui8matrix_ilu3_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------
{
    for (int i = i0 ; i <= i1 ; i++) {
        line_max3_ui8matrix_ilu3_red(X,i,j0,j1,Y);
    }
}
// -------------------------------------------------------------------------------
void max3_ui8matrix_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------
{
    int r = (i1 - i0 + 1) % 2;

    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_max3_ui8matrix_elu2_red(X,i,j0,j1,Y);
    }

    if (r) line_max3_ui8matrix_red(X,i1-r+1,j0,j1,Y);
    
}

// --------------------------------------------------------------------------------------
void max3_ui8matrix_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------
{

    int r = (i1 - i0 + 1) % 2;
    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_max3_ui8matrix_elu2_red_factor(X,i,j0,j1,Y);
    }
    if (r) line_max3_ui8matrix_red(X,i1-r+1,j0,j1,Y);
}
// ------------------------------------------------------------------------------------
void max3_ui8matrix_ilu3_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ------------------------------------------------------------------------------------
{
    int r = (i1 - i0 + 1) % 2;
    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_max3_ui8matrix_ilu3_elu2_red(X,i,j0,j1,Y);
    }
    if (r) line_max3_ui8matrix_red(X,i1-r+1,j0,j1,Y); 
}
// -------------------------------------------------------------------------------------------
void max3_ui8matrix_ilu3_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------------------
{
    int r = (i1 - i0 + 1) % 2;
    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_max3_ui8matrix_ilu3_elu2_red_factor(X,i,j0,j1,Y);
    }
    if (r) line_max3_ui8matrix_red(X,i1-r+1,j0,j1,Y); 
}
