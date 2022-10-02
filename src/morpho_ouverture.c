/* -------------------------- */
/* --- morpho_ouverture.c --- */
/* -------------------------- */

/*
 * Copyright (c) 2020 - 2021, Lionel Lacassagne, All rights reserved
  
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"
//#include "sequence.h"

#include "swp.h"
#include "morpho_min.h"
#include "morpho_max.h"


// -------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------
{  

    //  c1  c2  c3  c4   c5
    //  c6  b1  b2  b3   c7
    //  c8  b4  X   b5   c9
    //  c10 b6  b7  b8   c11
    //  c12 c13 c14 c15  c16

    uint8 c1, c6, c8, c10, c12, c2, b1, b4, b6, c13, c3, b2, x, b7, c14;
    uint8 c4,b3,b5,b8,c15;
    uint8 c5,c7,c9,c11,c16;
    uint8 min_1, min_2, min_3, min_4, min_5, min_6, min_7, min_8, min_9, max_1, max_2, max_3;
    uint8 tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7,tmp8, tmp9;
    uint8 looptmp1, looptmp2, looptmp3;
    uint8 tmp10, tmp11, tmp12;

    c1 = load2(X,i-2, j0-2);    c6 = load2(X,i-1,j0-2);     c8 = load2(X, i, j0-2); c10 = load2(X, i+1, j0-2);  c12 = load2(X, i+2, j0-2);
    c2 = load2(X, i-2, j0-1);   b1 = load2(X, i-1, j0-1);   b4 = load2(X, i, j0-1); b6 = load2(X, i+1, j0-1);   c13 = load2(X, i+2, j0-1);
    c3 = load2(X, i-2, j0);     b2 = load2(X, i-1, j0);     x = load2(X, i, j0);    b7 = load2(X, i+1, j0);     c14 = load2(X, i+2, j0);
    c4 = load2(X,i-2, j0+1);    b3 = load2(X,i-1,j0+1);     b5 = load2(X, i, j0+1); b8 = load2(X, i+1, j0+1);   c15 = load2(X, i+2, j0+1);

    tmp1 = min3(c1,c6,c8);      tmp2 = min3(c2,b1,b4);   tmp3 = min3(c3,b2,x);      tmp4 = min3(c4,b3,b5);
    tmp5 = min3(c6,c8,c10);     tmp6 = min3(b1, b4, b6); tmp7 = min3(b2,x,b7);      tmp8 = min3(b3,b5,b8);
    tmp9 = min3(c8,c10,c12);    tmp10 = min3(b4,b6,c13); tmp11 = min3(x,b7,c14);    tmp12 = min3(b5,b8,c15);

    min_1 = min3(tmp1,tmp2,tmp3);     min_4 = min3(tmp2,tmp3,tmp4);
    min_2 = min3(tmp5,tmp6,tmp7);     min_5 = min3(tmp6,tmp7,tmp8);
    min_3 = min3(tmp9,tmp10,tmp11);   min_6 = min3(tmp10,tmp11,tmp12);

    max_1 = max3(min_1,min_2,min_3);
    max_2 = max3(min_4,min_5,min_6);

    for (int j = j0 ; j <= j1 ; j += 1) {
        c5 = load2(X, i-2, j+2); c7 = load2(X, i-1, j+2); c9 = load2(X, i, j+2); c11 = load2(X, i+1, j+2); c16 = load2(X, i+2, j+2);

        looptmp1 = min3(c5,c7,c9);  looptmp2 = min3(c7,c9,c11); looptmp3 = min3(c9,c11,c16);
        min_7 = min3(tmp3, tmp4, looptmp1);   min_8 = min3(tmp7,tmp8,looptmp2); min_9 = min3(tmp11,tmp12,looptmp3);
        

        max_3 = max3(min_7,min_8,min_9);

        store2(Y,i,j,max3(max_1, max_2, max_3));   

        tmp3 = tmp4; tmp4 = looptmp1; tmp7 = tmp8; tmp8 = looptmp2; tmp11 = tmp12; tmp12 = looptmp3;
        max_1 = max_2; max_2 = max_3;
    }
}
// ----------------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion_ilu5_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------------
{   //                                  
    //  c01  c02  c03  c04   c05   c06  c08   c09   c10
    //  c11  c12  c13  c14   c15   c16  c17   c18   c19
    //  c20  c21  x01  x02   x03   x04  x05   c22   c23
    //  c24  c25  c26  c27   c28   c29  c30   c31   c32
    //  c33  c34  c35  c36   c37   c38  c39   c40   c41

    int r = (j1 - j0 + 1) % 5;

    uint8 c01,c11,c20,c24,c33;
    uint8 c02,c12,c21,c25,c34;
    uint8 c03,c13,x01,c26,c35;
    uint8 c04,c14,x02,c27,c36;
    uint8 c05,c15,x03,c28,c37;

    uint8 max_1, max_2, max_3;

    uint8 tmp1, tmp2, tmp3, tmp4, tmp5;
    uint8 tmp6, tmp7, tmp8, tmp9, tmp10;
    uint8 tmp11, tmp12, tmp13, tmp14, tmp15;

    uint8 min_1, min_2, min_3, min_4, min_5, min_6, min_7, min_8, min_9;

    // Prologue
    c01 = load2(X, i-2, j0-2); c11 = load2(X, i-1, j0-2); c20 = load2(X, i, j0-2); c24 = load2(X, i+1, j0-2); c33 = load2(X, i+2, j0-2);
    c02 = load2(X, i-2, j0-1); c12 = load2(X, i-1, j0-1); c21 = load2(X, i, j0-1); c25 = load2(X, i+1, j0-1); c34 = load2(X, i+2, j0-1);
    c03 = load2(X, i-2, j0);   c13 = load2(X, i-1, j0);   x01 = load2(X, i, j0);   c26 = load2(X, i+1, j0);   c35 = load2(X, i+2, j0);
    c04 = load2(X, i-2, j0+1); c14 = load2(X, i-1, j0+1); x02 = load2(X, i, j0+1); c27 = load2(X, i+1, j0+1); c36 = load2(X, i+2, j0+1);

    tmp1 = min3(c01,c11,c20)    ; tmp2 = min3(c02,c12,c21)  ;   tmp3 = min3(c03,c13,x01)    ;   tmp4 = min3(c04, c14, x02);
    tmp6 = min3(c11, c20, c24)  ; tmp7 = min3(c12,c21,c25)  ;   tmp8 = min3(c13,x01,c26)    ;   tmp9 = min3(c14,x02,c27);
    tmp11 = min3(c20,c24,c33)   ; tmp12 = min3(c21,c25,c34) ;   tmp13 = min3(x01,c26,c35)   ;   tmp14 = min3(x02,c27,c36);

    min_1 = min3(tmp1,tmp2,tmp3)    ;   min_4 = min3(tmp2,tmp3,tmp4)    ;
    min_2 = min3(tmp6,tmp7,tmp8)    ;   min_5 = min3(tmp7,tmp8,tmp9)    ;
    min_3 = min3(tmp11,tmp12,tmp13) ;   min_6 = min3(tmp12,tmp13,tmp14) ;
    
    max_1 = max3(min_1,min_2,min_3);
    max_2 = max3(min_4,min_5,min_6);

    for (int j = j0 ; j<=j1-r; j+=5) {
        
        // -------------------- 1st element --------------------
        c05 = load2(X, i-2, j+2); c15 = load2(X, i-1, j+2); x03 = load2(X, i, j+2); c28 = load2(X, i+1, j+2); c37 = load2(X, i+2, j+2);
        tmp5 = min3(c05,c15,x03); tmp10 = min3(c15,x03,c28); tmp15 = min3(x03,c28,c37);

        min_7 = min3(tmp3,tmp4,tmp5);   min_8 = min3(tmp8,tmp9,tmp10);  min_9 = min3(tmp13,tmp14,tmp15);
        max_3 = max3(min_7,min_8,min_9);

        store2(Y,i,j,max3(max_1,max_2,max_3));

        // -------------------- 2nd element --------------------
        c05 = load2(X, i-2, j+3); c15 = load2(X, i-1, j+3); x03 = load2(X, i, j+3); c28 = load2(X, i+1, j+3); c37 = load2(X, i+2, j+3);
        tmp1 = min3(c05,c15,x03); tmp6 = min3(c15,x03,c28); tmp11 = min3(x03,c28,c37);

        min_7 = min3(tmp4,tmp5,tmp1);   min_8 = min3(tmp9,tmp10,tmp6);  min_9 = min3(tmp14,tmp15,tmp11);
        max_1 = max3(min_7,min_8,min_9);

        store2(Y,i,j+1,max3(max_2,max_3,max_1));

        // -------------------- 3rd element --------------------
        c05 = load2(X, i-2, j+4); c15 = load2(X, i-1, j+4); x03 = load2(X, i, j+4); c28 = load2(X, i+1, j+4); c37 = load2(X, i+2, j+4);
        tmp2 = min3(c05,c15,x03); tmp7 = min3(c15,x03,c28); tmp12 = min3(x03,c28,c37);
        
        min_7 = min3(tmp5,tmp1,tmp2);   min_8 = min3(tmp10,tmp6,tmp7);  min_9 = min3(tmp15,tmp11,tmp12);
        max_2 = max3(min_7,min_8,min_9);

        store2(Y,i,j+2,max3(max_3,max_1,max_2));

        // -------------------- 4th element --------------------
        c05 = load2(X, i-2, j+5); c15 = load2(X, i-1, j+5); x03 = load2(X, i, j+5); c28 = load2(X, i+1, j+5); c37 = load2(X, i+2, j+5);
        tmp3 = min3(c05,c15,x03); tmp8 = min3(c15,x03,c28); tmp13 = min3(x03,c28,c37);

        min_7 = min3(tmp1,tmp2,tmp3);   min_8 = min3(tmp6,tmp7,tmp8);  min_9 = min3(tmp11,tmp12,tmp13);
        max_3 = max3(min_7,min_8,min_9);

        store2(Y,i,j+3,max3(max_1,max_2,max_3));

        // -------------------- 5th element --------------------
        c05 = load2(X, i-2, j+6); c15 = load2(X, i-1, j+6); x03 = load2(X, i, j+6); c28 = load2(X, i+1, j+6); c37 = load2(X, i+2, j+6);
        tmp4 = min3(c05,c15,x03); tmp9 = min3(c15,x03,c28); tmp14 = min3(x03,c28,c37);

        min_7 = min3(tmp2,tmp3,tmp4);   min_8 = min3(tmp7,tmp8,tmp9);  min_9 = min3(tmp12,tmp13,tmp14);
        max_1 = max3(min_7,min_8,min_9);


        store2(Y,i,j+4,max3(max_2,max_3, max_1));

        // ---------- Rotation pour l'iteration suivante -------------------
        min_1 = min3(tmp1,tmp2,tmp3)    ;   min_4 = min3(tmp2,tmp3,tmp4)     ;   
        min_2 = min3(tmp6,tmp7,tmp8)    ;   min_5 = min3(tmp7,tmp8,tmp9)     ;
        min_3 = min3(tmp11,tmp12,tmp13) ;   min_6 = min3(tmp12,tmp13,tmp14)  ;
        
        max_1 = max3(min_1, min_2, min_3);
        max_2 = max3(min_4,min_5,min_6);
    }

    if (r) line_ouverture3_ui8matrix_fusion(X,i,j1-r+1,j1,Y);

}
// ---------------------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion_ilu5_elu2_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------------
{
//  c01  c02  c03  c04   c05   c06  c08   c09   c10
//  c11  c12  c13  c14   c15   c16  c17   c18   c19
//  c20  c21  x01  x02   x03   x04  x05   c22   c23
//  c24  c25  x11  x12   x13   x14  x15   c26   c27
//  c28  c29  c30  c31   c32   c33  c34   c35   c36
//  c37  c38  c39  c40   c41   c42  c43   c44   c45

    uint8 c01,c11,c20,c24,c28,c37;
    uint8 c02,c12,c21,c25,c29,c38;
    uint8 c03,c13,x01,x11,c30,c39;
    uint8 c04,c14,x02,x12,c31,c40;

    uint8 c05,c15,x03,x13,c32,c41;

    uint8 tmp1,tmp2,tmp3,tmp4,tmp5;
    uint8 tmp6,tmp7,tmp8,tmp9,tmp10;
    uint8 tmp11,tmp12,tmp13,tmp14,tmp15;
    uint8 tmp16,tmp17,tmp18,tmp19,tmp20;

    uint8 min_1,min_2,min_3,min_4,min_5,min_6,min_7,min_8,min_9,min_10,min_11,min_12;
    uint8 max_1,max_2,max_3,max_4,max_5,max_6;

    // loading of pixels
    c01 = load2(X, i-2, j0-2);  c11 = load2(X,i-1,j0-2);    c20 = load2(X,i,j0-2);  c24 = load2(X,i+1,j0-2);    c28 = load2(X,i+2,j0-2);    c37 = load2(X,i+3,j0-2);
    c02 = load2(X, i-2, j0-1);  c12 = load2(X,i-1,j0-1);    c21 = load2(X,i,j0-1);  c25 = load2(X,i+1,j0-1);    c29 = load2(X,i+2,j0-1);    c38 = load2(X,i+3,j0-1);
    c03 = load2(X, i-2, j0  );  c13 = load2(X,i-1,j0  );    x01 = load2(X,i,j0  );  x11 = load2(X,i+1,j0  );    c30 = load2(X,i+2,j0  );    c39 = load2(X,i+3,j0  );
    c04 = load2(X, i-2, j0+1);  c14 = load2(X,i-1,j0+1);    x02 = load2(X,i,j0+1);  x12 = load2(X,i+1,j0+1);    c31 = load2(X,i+2,j0+1);    c40 = load2(X,i+3,j0+1);

    // calculation of min values
    tmp1 = min3(c01,c11,c20);   tmp2 = min3(c02,c12,c21);   tmp3 = min3(c03,c13,x01);   tmp4 = min3(c04,c14,x02);
    tmp6 = min3(c11,c20,c24);   tmp7 = min3(c12,c21,c25);   tmp8 = min3(c13,x01,x11);   tmp9 = min3(c14,x02,x12);
    tmp11 = min3(c20,c24,c28);  tmp12 = min3(c21,c25,c29);  tmp13 = min3(x01,x11,c30);  tmp14 = min3(x02,x12,c31);
    tmp16 = min3(c24,c28,c37);  tmp17 = min3(c25,c29,c38);  tmp18 = min3(x11,c30,c39);  tmp19 = min3(x12,c31,c40);
    
    min_1 = min3(tmp1,tmp2,tmp3);       min_5 = min3(tmp2,tmp3,tmp4);   
    min_2 = min3(tmp6,tmp7,tmp8);       min_6 = min3(tmp7,tmp8,tmp9);
    min_3 = min3(tmp11,tmp12,tmp13);    min_7 = min3(tmp12,tmp13,tmp14);
    min_4 = min3(tmp16,tmp17,tmp18);    min_8 = min3(tmp17,tmp18,tmp19);

    // calculation of max values
    max_1 = max3(min_1,min_2,min_3);    max_2 = max3(min_5,min_6,min_7);
    max_4 = max3(min_2,min_3,min_4);    max_5 = max3(min_6,min_7,min_8);

    int r = (j1 - j0 + 1) % 5;

    for (int j = j0 ; j <= j1-r ; j+=5) {
        //---------------------- 1st two elements ------------------------
        c05 = load2(X, i-2, j+2);  c15 = load2(X,i-1,j+2);    x03 = load2(X,i,j+2);  x13 = load2(X,i+1,j+2);    c32 = load2(X,i+2,j+2);    c41 = load2(X,i+3,j+2);
        tmp5 = min3(c05,c15,x03);   tmp10 = min3(c15,x03,x13);  tmp15 = min3(x03,x13,c32);  tmp20 = min3(x13,c32,c41);

        min_9 = min3(tmp3,tmp4,tmp5);   
        min_10 = min3(tmp8,tmp9,tmp10); 
        min_11 = min3(tmp13,tmp14,tmp15);
        min_12 = min3(tmp18,tmp19,tmp20);

        max_3 = max3(min_9,min_10,min_11);
        max_6 = max3(min_10,min_11,min_12);

        store2(Y,i,j,max3(max_1,max_2,max_3));
        store2(Y,i+1,j,max3(max_4,max_5,max_6));

        //---------------------- 2nd two elements ------------------------
        c05 = load2(X, i-2, j+3);  c15 = load2(X,i-1,j+3);    x03 = load2(X,i,j+3);  x13 = load2(X,i+1,j+3);    c32 = load2(X,i+2,j+3);    c41 = load2(X,i+3,j+3);
        tmp1 = min3(c05,c15,x03);   tmp6 = min3(c15,x03,x13);  tmp11 = min3(x03,x13,c32);  tmp16 = min3(x13,c32,c41);

        min_9 = min3(tmp4,tmp5,tmp1);   
        min_10 = min3(tmp9,tmp10,tmp6); 
        min_11 = min3(tmp14,tmp15,tmp11);
        min_12 = min3(tmp19,tmp20,tmp16);

        max_1 = max3(min_9,min_10,min_11);
        max_4 = max3(min_10,min_11,min_12);

        store2(Y,i,j+1,max3(max_2,max_3,max_1));
        store2(Y,i+1,j+1,max3(max_5,max_6,max_4));

        //---------------------- 3rd two elements ------------------------
        c05 = load2(X, i-2, j+4);  c15 = load2(X,i-1,j+4);    x03 = load2(X,i,j+4);  x13 = load2(X,i+1,j+4);    c32 = load2(X,i+2,j+4);    c41 = load2(X,i+3,j+4);
        tmp2 = min3(c05,c15,x03);   tmp7 = min3(c15,x03,x13);  tmp12 = min3(x03,x13,c32);  tmp17 = min3(x13,c32,c41);

        min_9 = min3(tmp5,tmp1,tmp2);   
        min_10 = min3(tmp10,tmp6,tmp7); 
        min_11 = min3(tmp15,tmp11,tmp12);
        min_12 = min3(tmp20,tmp16,tmp17);

        max_2 = max3(min_9,min_10,min_11);
        max_5 = max3(min_10,min_11,min_12);

        store2(Y,i,j+2,max3(max_3,max_1, max_2));
        store2(Y,i+1,j+2,max3(max_6,max_4,max_5));

        //---------------------- 4th two elements ------------------------
        c05 = load2(X, i-2, j+5);  c15 = load2(X,i-1,j+5);    x03 = load2(X,i,j+5);  x13 = load2(X,i+1,j+5);    c32 = load2(X,i+2,j+5);    c41 = load2(X,i+3,j+5);
        tmp3 = min3(c05,c15,x03);   tmp8 = min3(c15,x03,x13);  tmp13 = min3(x03,x13,c32);  tmp18 = min3(x13,c32,c41);

        min_9 = min3(tmp1,tmp2,tmp3);   
        min_10 = min3(tmp6,tmp7,tmp8); 
        min_11 = min3(tmp11,tmp12,tmp13);
        min_12 = min3(tmp16,tmp17,tmp18);


        max_3 = max3(min_9,min_10,min_11);
        max_6 = max3(min_10,min_11,min_12);

        store2(Y,i,j+3,max3(max_1,max_2,max_3));
        store2(Y,i+1,j+3,max3(max_4,max_5,max_6));

        //---------------------- 5th two elements ------------------------
        c05 = load2(X, i-2, j+6);  c15 = load2(X,i-1,j+6);    x03 = load2(X,i,j+6);  x13 = load2(X,i+1,j+6);    c32 = load2(X,i+2,j+6);    c41 = load2(X,i+3,j+6);
        tmp4 = min3(c05,c15,x03);   tmp9 = min3(c15,x03,x13);  tmp14 = min3(x03,x13,c32);  tmp19 = min3(x13,c32,c41);

        min_9 = min3(tmp2,tmp3,tmp4);   
        min_10 = min3(tmp7,tmp8,tmp9); 
        min_11 = min3(tmp12,tmp13,tmp14);
        min_12 = min3(tmp17,tmp18,tmp19);

        max_1 = max3(min_9,min_10,min_11);
        max_4 = max3(min_10,min_11,min_12);

        store2(Y,i,j+4,max3(max_2,max_3,max_1));
        store2(Y,i+1,j+4,max3(max_5,max_6,max_4));


        // ---------- Rotation pour l'iteration suivante -------------------

        min_1 = min3(tmp1,tmp2,tmp3);       min_5 = min3(tmp2,tmp3,tmp4);   
        min_2 = min3(tmp6,tmp7,tmp8);       min_6 = min3(tmp7,tmp8,tmp9);
        min_3 = min3(tmp11,tmp12,tmp13);    min_7 = min3(tmp12,tmp13,tmp14);
        min_4 = min3(tmp16,tmp17,tmp18);    min_8 = min3(tmp17,tmp18,tmp19);

        max_1 = max3(min_1,min_2,min_3);    max_2 = max3(min_5,min_6,min_7);
        max_4 = max3(min_2,min_3,min_4);    max_5 = max3(min_6,min_7,min_8);
    }

    if (r) {
            line_ouverture3_ui8matrix_fusion(X,i,j1-r+1,j1,Y);
            line_ouverture3_ui8matrix_fusion(X,i+1,j1-r+1,j1,Y);

           }

}

// ----------------------------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion_ilu5_elu2_red_factor(uint8 **X, int i, int j0, int j1, uint8 **Y)
// ----------------------------------------------------------------------------------------------------
{
//  c01  c02  c03  c04   c05   c06  c08   c09   c10
//  c11  c12  c13  c14   c15   c16  c17   c18   c19
//  c20  c21  x01  x02   x03   x04  x05   c22   c23
//  c24  c25  x11  x12   x13   x14  x15   c26   c27
//  c28  c29  c30  c31   c32   c33  c34   c35   c36
//  c37  c38  c39  c40   c41   c42  c43   c44   c45

    uint8 c01,c11,c20,c24,c28,c37;
    uint8 c02,c12,c21,c25,c29,c38;
    uint8 c03,c13,x01,x11,c30,c39;
    uint8 c04,c14,x02,x12,c31,c40;

    uint8 c05,c15,x03,x13,c32,c41;

    uint8 tmp1,tmp2,tmp3,tmp4,tmp5;
    uint8 tmp6,tmp7,tmp8,tmp9,tmp10;
    uint8 tmp11,tmp12,tmp13,tmp14,tmp15;
    uint8 tmp16,tmp17,tmp18,tmp19,tmp20;

    uint8 min_1,min_2,min_3,min_4,min_5,min_6,min_7,min_8,min_9,min_10,min_11,min_12;
    uint8 fact1,fact2,fact3,fact4,fact5,fact6,fact7,fact8,fact9,fact10,fact11,fact12;
    uint8 max_1,max_2,max_3,max_4,max_5,max_6;
    uint8 facttmp1,facttmp2,facttmp3,facttmp4;
    uint8 factmin1,factmin2;
    uint8 factmax1,factmax2;

    // loading of pixels
    c01 = load2(X, i-2, j0-2);  c11 = load2(X,i-1,j0-2);    c20 = load2(X,i,j0-2);  c24 = load2(X,i+1,j0-2);    c28 = load2(X,i+2,j0-2);    c37 = load2(X,i+3,j0-2);
    c02 = load2(X, i-2, j0-1);  c12 = load2(X,i-1,j0-1);    c21 = load2(X,i,j0-1);  c25 = load2(X,i+1,j0-1);    c29 = load2(X,i+2,j0-1);    c38 = load2(X,i+3,j0-1);
    c03 = load2(X, i-2, j0  );  c13 = load2(X,i-1,j0  );    x01 = load2(X,i,j0  );  x11 = load2(X,i+1,j0  );    c30 = load2(X,i+2,j0  );    c39 = load2(X,i+3,j0  );
    c04 = load2(X, i-2, j0+1);  c14 = load2(X,i-1,j0+1);    x02 = load2(X,i,j0+1);  x12 = load2(X,i+1,j0+1);    c31 = load2(X,i+2,j0+1);    c40 = load2(X,i+3,j0+1);


    // calculation of factors
    fact1 = min2(c11,c20);  fact2 = min2(c12,c21);  fact3 = min2(c13,x01)   ;  fact4 = min2(c14,x02);
    fact5 = min2(c20,c24);  fact6 = min2(c21,c25);  fact7 = min2(x01,x11)   ;  fact8 = min2(x02,x12);
    fact9 = min2(c24,c28);  fact10 = min2(c25,c29); fact11 = min2(x11,c30)  ; fact12 = min2(x12,c31);

    // calculation of min values
    tmp1 = min2(c01,fact1);   tmp2 = min2(c02,fact2);   tmp3 = min2(c03,fact3);   tmp4 = min2(c04,fact4);
    tmp6 = min2(fact1,c24);   tmp7 = min2(fact2,c25);   tmp8 = min2(fact3,x11);   tmp9 = min2(fact4,x12);
    tmp11 = min2(fact5,c28);  tmp12 = min2(fact6,c29);  tmp13 = min2(fact7,c30);  tmp14 = min2(fact8,c31);
    tmp16 = min2(fact9,c37);  tmp17 = min2(fact10,c38);  tmp18 = min2(fact11,c39);  tmp19 = min2(fact12,c40);
    

    facttmp1 = min2(tmp2,tmp3);    facttmp2 = min2(tmp7,tmp8);    facttmp3 = min2(tmp12,tmp13);  facttmp4 = min2(tmp17,tmp18);

    min_1 = min2(tmp1,facttmp1);       min_5 = min2(facttmp1,tmp4);   
    min_2 = min2(tmp6,facttmp2);       min_6 = min2(facttmp2,tmp9);
    min_3 = min2(tmp11,facttmp3);      min_7 = min2(facttmp3,tmp14);
    min_4 = min2(tmp16,facttmp4);      min_8 = min2(facttmp4,tmp19);

    // calculation of max values
    factmin1 = max2(min_2,min_3);  factmin2 = max2(min_6,min_7);

    max_1 = max2(min_1,factmin1);    max_2 = max2(min_5,factmin2);
    max_4 = max2(factmin1,min_4);    max_5 = max2(factmin2,min_8);

    int r = (j1 - j0 + 1) % 5;

    for (int j = j0 ; j <= j1-r ; j+=5) {
        //---------------------- 1st two elements ------------------------
        c05 = load2(X, i-2, j+2);  c15 = load2(X,i-1,j+2);    x03 = load2(X,i,j+2);  x13 = load2(X,i+1,j+2);    c32 = load2(X,i+2,j+2);    c41 = load2(X,i+3,j+2);
        fact7 = min2(c15,x03);  fact8 = min2(x13,c32);
        tmp5 = min2(c05,fact7);   tmp10 = min2(fact7,x13);  tmp15 = min2(x03,fact8);  tmp20 = min2(fact8,c41);

        facttmp1 = min2(tmp4,tmp5);    facttmp2 = min2(tmp9,tmp10);   facttmp3 = min2(tmp14,tmp15);  facttmp4 = min2(tmp19,tmp20);

        min_9 = min2(tmp3,facttmp1);   
        min_10 = min2(tmp8,facttmp2); 
        min_11 = min2(tmp13,facttmp3);
        min_12 = min2(tmp18,facttmp4);

        factmin1 = max2(min_10,min_11);

        max_3 = max2(min_9,factmin1);
        max_6 = max2(factmin1,min_12);

        factmax1 = max2(max_2,max_3);
        factmax2 = max2(max_5,max_6);

        store2(Y,i,j,max2(max_1,factmax1));
        store2(Y,i+1,j,max2(max_4,factmax2));

        //---------------------- 2nd two elements ------------------------
        c05 = load2(X, i-2, j+3);  c15 = load2(X,i-1,j+3);    x03 = load2(X,i,j+3);  x13 = load2(X,i+1,j+3);    c32 = load2(X,i+2,j+3);    c41 = load2(X,i+3,j+3);
        fact7 = min2(c15,x03);  fact8 = min2(x13,c32);
        tmp1 = min2(c05,fact7);   tmp6 = min2(fact7,x13);  tmp11 = min2(x03,fact8);  tmp16 = min2(fact8,c41);

        min_9 = min2(facttmp1,tmp1);   
        min_10 = min2(facttmp2,tmp6); 
        min_11 = min2(facttmp3,tmp11);
        min_12 = min2(facttmp4,tmp16);

        factmin1 = max2(min_10,min_11);

        max_1 = max2(min_9,factmin1);
        max_4 = max2(factmin1,min_12);

        store2(Y,i,j+1,max2(factmax1,max_1));
        store2(Y,i+1,j+1,max2(factmax2,max_4));

        //---------------------- 3rd two elements ------------------------
        c05 = load2(X, i-2, j+4);  c15 = load2(X,i-1,j+4);    x03 = load2(X,i,j+4);  x13 = load2(X,i+1,j+4);    c32 = load2(X,i+2,j+4);    c41 = load2(X,i+3,j+4);
        fact7 = min2(c15,x03);  fact8 = min2(x13,c32);
        tmp2 = min2(c05,fact7);   tmp7 = min2(fact7,x13);  tmp12 = min2(x03,fact8);  tmp17 = min2(fact8,c41);


        facttmp1 = min2(tmp1,tmp2);    facttmp2 = min2(tmp6,tmp7);   facttmp3 = min2(tmp11,tmp12);  facttmp4 = min2(tmp16,tmp17);


        min_9 = min2(tmp5,facttmp1);   
        min_10 = min2(tmp10,facttmp2); 
        min_11 = min2(tmp15,facttmp3);
        min_12 = min2(tmp20,facttmp4);

        factmin1 = max2(min_10,min_11);

        max_2 = max2(min_9,factmin1);
        max_5 = max2(factmin1,min_12);

        factmax1 = max2(max_1,max_2);
        factmax2 = max2(max_4,max_5);

        store2(Y,i,j+2,max2(max_3,factmax1));
        store2(Y,i+1,j+2,max2(max_6,factmax2));

        //---------------------- 4th two elements ------------------------
        c05 = load2(X, i-2, j+5);  c15 = load2(X,i-1,j+5);    x03 = load2(X,i,j+5);  x13 = load2(X,i+1,j+5);    c32 = load2(X,i+2,j+5);    c41 = load2(X,i+3,j+5);
        fact7 = min2(c15,x03);  fact8 = min2(x13,c32);
        tmp3 = min2(c05,fact7);   tmp8 = min2(fact7,x13);  tmp13 = min2(x03,fact8);  tmp18 = min2(fact8,c41);

        min_9 = min2(facttmp1,tmp3);   
        min_10 = min2(facttmp2,tmp8); 
        min_11 = min2(facttmp3,tmp13);
        min_12 = min2(facttmp4,tmp18);

        factmin1 = max2(min_10,min_11);

        max_3 = max2(min_9,factmin1);
        max_6 = max2(factmin1,min_12);

        store2(Y,i,j+3,max2(factmax1,max_3));
        store2(Y,i+1,j+3,max2(factmax2,max_6));

        //---------------------- 5th two elements ------------------------
        c05 = load2(X, i-2, j+6);  c15 = load2(X,i-1,j+6);    x03 = load2(X,i,j+6);  x13 = load2(X,i+1,j+6);    c32 = load2(X,i+2,j+6);    c41 = load2(X,i+3,j+6);
        fact7 = min2(c15,x03);  fact8 = min2(x13,c32);
        tmp4 = min2(c05,fact7);   tmp9 = min2(fact7,x13);  tmp14 = min2(x03,fact8);  tmp19 = min2(fact8,c41);

        facttmp1 = min2(tmp3,tmp4);    facttmp2 = min2(tmp8,tmp9);   facttmp3 = min2(tmp13,tmp14);  facttmp4 = min2(tmp18,tmp19);

        min_9 = min2(tmp2,facttmp1);   
        min_10 = min2(tmp7,facttmp2); 
        min_11 = min2(tmp12,facttmp3);
        min_12 = min2(tmp17,facttmp4);

        factmin1 = max2(min_10,min_11);

        max_1 = max2(min_9,factmin1);
        max_4 = max2(factmin1,min_12);

        factmax1 = max2(max_3,max_1);
        factmax2 = max2(max_6,max_4);

        store2(Y,i,j+4,max2(max_2,factmax1));
        store2(Y,i+1,j+4,max2(max_5,factmax2));


        // ---------- Rotation pour l'iteration suivante -------------------

        facttmp1 = min2(tmp2,tmp3); facttmp2 = min2(tmp7,tmp8); facttmp3 = min2(tmp12,tmp13);   facttmp4 = min2(tmp17,tmp18);
        min_1 = min2(tmp1,facttmp1);       min_5 = min2(facttmp1,tmp4);   
        min_2 = min2(tmp6,facttmp2);       min_6 = min2(facttmp2,tmp9);
        min_3 = min2(tmp11,facttmp3);    min_7 = min2(facttmp3,tmp14);
        min_4 = min2(tmp16,facttmp4);    min_8 = min2(facttmp4,tmp19);

        factmax1 = max2(min_2,min_3);   factmax2 = max2(min_6,min_7);
        max_1 = max2(min_1,factmax1);    max_2 = max2(min_5,factmax2);
        max_4 = max2(factmax1,min_4);    max_5 = max2(factmax2,min_8);
    }

    if (r) {
            line_ouverture3_ui8matrix_fusion_ilu5_elu2_red(X,i,j1-r+1,j1,Y);
           }
}
// -----------------------------------------------------------------------------------------
void line_ouverture3_ui8matrix_fusion_ilu15_red(uint8 **X, int i, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------------
{

//                                  
    //  c01  c02  c03  c04   c05   c06  c08   c09   c10
    //  c11  c12  c13  c14   c15   c16  c17   c18   c19
    //  c20  c21  x01  x02   x03   x04  x05   c22   c23
    //  c24  c25  c26  c27   c28   c29  c30   c31   c32
    //  c33  c34  c35  c36   c37   c38  c39   c40   c41

    int r = (j1 - j0 + 1) % 15;

    uint8 c01,c11,c20,c24,c33;
    uint8 c02,c12,c21,c25,c34;
    uint8 c03,c13,x01,c26,c35;
    uint8 c04,c14,x02,c27,c36;
    uint8 c05,c15,x03,c28,c37;

    uint8 max_1, max_2, max_3;

    uint8 tmp1, tmp2, tmp3, tmp4, tmp5;
    uint8 tmp6, tmp7, tmp8, tmp9, tmp10;
    uint8 tmp11, tmp12, tmp13, tmp14, tmp15;

    uint8 min_1, min_2, min_3, min_4, min_5, min_6, min_7, min_8, min_9;

    // Prologue
    c01 = load2(X, i-2, j0-2); c11 = load2(X, i-1, j0-2); c20 = load2(X, i, j0-2); c24 = load2(X, i+1, j0-2); c33 = load2(X, i+2, j0-2);
    c02 = load2(X, i-2, j0-1); c12 = load2(X, i-1, j0-1); c21 = load2(X, i, j0-1); c25 = load2(X, i+1, j0-1); c34 = load2(X, i+2, j0-1);
    c03 = load2(X, i-2, j0);   c13 = load2(X, i-1, j0);   x01 = load2(X, i, j0);   c26 = load2(X, i+1, j0);   c35 = load2(X, i+2, j0);
    c04 = load2(X, i-2, j0+1); c14 = load2(X, i-1, j0+1); x02 = load2(X, i, j0+1); c27 = load2(X, i+1, j0+1); c36 = load2(X, i+2, j0+1);

    tmp1 = min3(c01,c11,c20)    ; tmp2 = min3(c02,c12,c21)  ;   tmp3 = min3(c03,c13,x01)    ;   tmp4 = min3(c04, c14, x02);
    tmp6 = min3(c11, c20, c24)  ; tmp7 = min3(c12,c21,c25)  ;   tmp8 = min3(c13,x01,c26)    ;   tmp9 = min3(c14,x02,c27);
    tmp11 = min3(c20,c24,c33)   ; tmp12 = min3(c21,c25,c34) ;   tmp13 = min3(x01,c26,c35)   ;   tmp14 = min3(x02,c27,c36);

    min_1 = min3(tmp1,tmp2,tmp3)    ;   min_4 = min3(tmp2,tmp3,tmp4)    ;
    min_2 = min3(tmp6,tmp7,tmp8)    ;   min_5 = min3(tmp7,tmp8,tmp9)    ;
    min_3 = min3(tmp11,tmp12,tmp13) ;   min_6 = min3(tmp12,tmp13,tmp14) ;
    
    max_1 = max3(min_1,min_2,min_3);
    max_2 = max3(min_4,min_5,min_6);

    for (int j = j0 ; j<=j1-r; j+=15) {
        
        // -------------------- 1st element --------------------
        c05 = load2(X, i-2, j+2); c15 = load2(X, i-1, j+2); x03 = load2(X, i, j+2); c28 = load2(X, i+1, j+2); c37 = load2(X, i+2, j+2);
        tmp5 = min3(c05,c15,x03); tmp10 = min3(c15,x03,c28); tmp15 = min3(x03,c28,c37);

        min_7 = min3(tmp3,tmp4,tmp5);   min_8 = min3(tmp8,tmp9,tmp10);  min_9 = min3(tmp13,tmp14,tmp15);
        max_3 = max3(min_7,min_8,min_9);

        store2(Y,i,j,max3(max_1,max_2,max_3));

        // -------------------- 2nd element --------------------
        c05 = load2(X, i-2, j+3); c15 = load2(X, i-1, j+3); x03 = load2(X, i, j+3); c28 = load2(X, i+1, j+3); c37 = load2(X, i+2, j+3);
        tmp1 = min3(c05,c15,x03); tmp6 = min3(c15,x03,c28); tmp11 = min3(x03,c28,c37);

        min_7 = min3(tmp4,tmp5,tmp1);   min_8 = min3(tmp9,tmp10,tmp6);  min_9 = min3(tmp14,tmp15,tmp11);
        max_1 = max3(min_7,min_8,min_9);

        store2(Y,i,j+1,max3(max_2,max_3,max_1));

        // -------------------- 3rd element --------------------
        c05 = load2(X, i-2, j+4); c15 = load2(X, i-1, j+4); x03 = load2(X, i, j+4); c28 = load2(X, i+1, j+4); c37 = load2(X, i+2, j+4);
        tmp2 = min3(c05,c15,x03); tmp7 = min3(c15,x03,c28); tmp12 = min3(x03,c28,c37);
        
        min_7 = min3(tmp5,tmp1,tmp2);   min_8 = min3(tmp10,tmp6,tmp7);  min_9 = min3(tmp15,tmp11,tmp12);
        max_2 = max3(min_7,min_8,min_9);

        store2(Y,i,j+2,max3(max_3,max_1,max_2));

        // -------------------- 4th element --------------------
        c05 = load2(X, i-2, j+5); c15 = load2(X, i-1, j+5); x03 = load2(X, i, j+5); c28 = load2(X, i+1, j+5); c37 = load2(X, i+2, j+5);
        tmp3 = min3(c05,c15,x03); tmp8 = min3(c15,x03,c28); tmp13 = min3(x03,c28,c37);

        min_7 = min3(tmp1,tmp2,tmp3);   min_8 = min3(tmp6,tmp7,tmp8);  min_9 = min3(tmp11,tmp12,tmp13);
        max_3 = max3(min_7,min_8,min_9);

        store2(Y,i,j+3,max3(max_1,max_2,max_3));

        // -------------------- 5th element --------------------
        c05 = load2(X, i-2, j+6); c15 = load2(X, i-1, j+6); x03 = load2(X, i, j+6); c28 = load2(X, i+1, j+6); c37 = load2(X, i+2, j+6);
        tmp4 = min3(c05,c15,x03); tmp9 = min3(c15,x03,c28); tmp14 = min3(x03,c28,c37);

        min_7 = min3(tmp2,tmp3,tmp4);   min_8 = min3(tmp7,tmp8,tmp9);  min_9 = min3(tmp12,tmp13,tmp14);
        max_1 = max3(min_7,min_8,min_9);


        store2(Y,i,j+4,max3(max_2,max_3, max_1));

        // -------------------- 6th element --------------------
        c05 = load2(X, i-2, j+7); c15 = load2(X, i-1, j+7); x03 = load2(X, i, j+7); c28 = load2(X, i+1, j+7); c37 = load2(X, i+2, j+7);
        tmp5 = min3(c05,c15,x03); tmp10 = min3(c15,x03,c28); tmp15 = min3(x03,c28,c37);

        min_7 = min3(tmp3,tmp4,tmp5);   min_8 = min3(tmp8,tmp9,tmp10);  min_9 = min3(tmp13,tmp14,tmp15);
        max_2 = max3(min_7,min_8,min_9);

        store2(Y,i,j+5,max3(max_3, max_1,max_2));

        // -------------------- 7th element --------------------
        c05 = load2(X, i-2, j+8); c15 = load2(X, i-1, j+8); x03 = load2(X, i, j+8); c28 = load2(X, i+1, j+8); c37 = load2(X, i+2, j+8);
        tmp1 = min3(c05,c15,x03); tmp6 = min3(c15,x03,c28); tmp11 = min3(x03,c28,c37);

        min_7 = min3(tmp4,tmp5,tmp1);   min_8 = min3(tmp9,tmp10,tmp6);  min_9 = min3(tmp14,tmp15,tmp11);
        max_3 = max3(min_7,min_8,min_9);

        store2(Y,i,j+6,max3(max_2,max_3,max_1));

        // -------------------- 8th element --------------------
        c05 = load2(X, i-2, j+9); c15 = load2(X, i-1, j+9); x03 = load2(X, i, j+9); c28 = load2(X, i+1, j+9); c37 = load2(X, i+2, j+9);
        tmp2 = min3(c05,c15,x03); tmp7 = min3(c15,x03,c28); tmp12 = min3(x03,c28,c37);
        
        min_7 = min3(tmp5,tmp1,tmp2);   min_8 = min3(tmp10,tmp6,tmp7);  min_9 = min3(tmp15,tmp11,tmp12);
        max_1 = max3(min_7,min_8,min_9);

        store2(Y,i,j+7,max3(max_3,max_1,max_2));

        // -------------------- 9th element --------------------
        c05 = load2(X, i-2, j+10); c15 = load2(X, i-1, j+10); x03 = load2(X, i, j+10); c28 = load2(X, i+1, j+10); c37 = load2(X, i+2, j+10);
        tmp3 = min3(c05,c15,x03); tmp8 = min3(c15,x03,c28); tmp13 = min3(x03,c28,c37);

        min_7 = min3(tmp1,tmp2,tmp3);   min_8 = min3(tmp6,tmp7,tmp8);  min_9 = min3(tmp11,tmp12,tmp13);
        max_2 = max3(min_7,min_8,min_9);

        store2(Y,i,j+8,max3(max_1,max_2,max_3));

        // -------------------- 10th element --------------------

        c05 = load2(X, i-2, j+11); c15 = load2(X, i-1, j+11); x03 = load2(X, i, j+11); c28 = load2(X, i+1, j+11); c37 = load2(X, i+2, j+11);
        tmp4 = min3(c05,c15,x03); tmp9 = min3(c15,x03,c28); tmp14 = min3(x03,c28,c37);

        min_7 = min3(tmp2,tmp3,tmp4);   min_8 = min3(tmp7,tmp8,tmp9);  min_9 = min3(tmp12,tmp13,tmp14);
        max_3 = max3(min_7,min_8,min_9);

        store2(Y,i,j+9,max3(max_2,max_3,max_1));

        // -------------------- 11th element --------------------
        c05 = load2(X, i-2, j+12); c15 = load2(X, i-1, j+12); x03 = load2(X, i, j+12); c28 = load2(X, i+1, j+12); c37 = load2(X, i+2, j+12);
        tmp5 = min3(c05,c15,x03); tmp10 = min3(c15,x03,c28); tmp15 = min3(x03,c28,c37);

        min_7 = min3(tmp3,tmp4,tmp5);   min_8 = min3(tmp8,tmp9,tmp10);  min_9 = min3(tmp13,tmp14,tmp15);
        max_1 = max3(min_7,min_8,min_9);

        store2(Y,i,j+10,max3(max_3, max_1,max_2));

        // -------------------- 12th element --------------------

        c05 = load2(X, i-2, j+13); c15 = load2(X, i-1, j+13); x03 = load2(X, i, j+13); c28 = load2(X, i+1, j+13); c37 = load2(X, i+2, j+13);
        tmp1 = min3(c05,c15,x03); tmp6 = min3(c15,x03,c28); tmp11 = min3(x03,c28,c37);

        min_7 = min3(tmp4,tmp5,tmp1);   min_8 = min3(tmp9,tmp10,tmp6);  min_9 = min3(tmp14,tmp15,tmp11);
        max_2 = max3(min_7,min_8,min_9);

        store2(Y,i,j+11,max3(max_2,max_3,max_1));

        // -------------------- 13th element --------------------

        c05 = load2(X, i-2, j+14); c15 = load2(X, i-1, j+14); x03 = load2(X, i, j+14); c28 = load2(X, i+1, j+14); c37 = load2(X, i+2, j+14);
        tmp2 = min3(c05,c15,x03); tmp7 = min3(c15,x03,c28); tmp12 = min3(x03,c28,c37);
        
        min_7 = min3(tmp5,tmp1,tmp2);   min_8 = min3(tmp10,tmp6,tmp7);  min_9 = min3(tmp15,tmp11,tmp12);
        max_3 = max3(min_7,min_8,min_9);

        store2(Y,i,j+12,max3(max_3,max_1,max_2));

        // -------------------- 14th element --------------------

        c05 = load2(X, i-2, j+15); c15 = load2(X, i-1, j+15); x03 = load2(X, i,j+15); c28 = load2(X, i+1, j+15); c37 = load2(X, i+2, j+15);
        tmp3 = min3(c05,c15,x03); tmp8 = min3(c15,x03,c28); tmp13 = min3(x03,c28,c37);

        min_7 = min3(tmp1,tmp2,tmp3);   min_8 = min3(tmp6,tmp7,tmp8);  min_9 = min3(tmp11,tmp12,tmp13);
        max_1 = max3(min_7,min_8,min_9);

        store2(Y,i,j+13,max3(max_1,max_2,max_3));

        // -------------------- 15th element --------------------

        c05 = load2(X, i-2, j+16); c15 = load2(X, i-1, j+16); x03 = load2(X, i, j+16); c28 = load2(X, i+1, j+16); c37 = load2(X, i+2, j+16);
        tmp4 = min3(c05,c15,x03); tmp9 = min3(c15,x03,c28); tmp14 = min3(x03,c28,c37);

        min_7 = min3(tmp2,tmp3,tmp4);   min_8 = min3(tmp7,tmp8,tmp9);  min_9 = min3(tmp12,tmp13,tmp14);
        max_2 = max3(min_7,min_8,min_9);

        store2(Y,i,j+14,max3(max_2,max_3,max_1));
        

        // ---------- Rotation pour l'iteration suivante -------------------
        min_1 = min3(tmp1,tmp2,tmp3)    ;   min_4 = min3(tmp2,tmp3,tmp4)     ;   
        min_2 = min3(tmp6,tmp7,tmp8)    ;   min_5 = min3(tmp7,tmp8,tmp9)     ;
        min_3 = min3(tmp11,tmp12,tmp13) ;   min_6 = min3(tmp12,tmp13,tmp14)  ;
        
        max_1 = max3(min_1, min_2, min_3);
        max_2 = max3(min_4,min_5,min_6);
    }

    if (r) line_ouverture3_ui8matrix_fusion(X,i,j1-r+1,j1,Y);

}
// ---------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y, uint8 **Z)
// ---------------------------------------------------------------------------------------------
{
    min3_ui8matrix_basic(X, i0-1, i1+1, j0-1, j1+1, Y);
    max3_ui8matrix_basic(Y, i0,   i1,   j0,   j1,   Z);
}
// -----------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------
{
    for (int i = i0 ; i <= i1 ; i++) {
        line_ouverture3_ui8matrix_fusion(X,i,j0,j1,Y);
    }
}
// --------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion_ilu5_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------------
{
    for (int i = i0 ; i <= i1 ; i++) {
        line_ouverture3_ui8matrix_fusion_ilu5_red(X,i,j0,j1,Y);
    }

}
// -------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion_ilu5_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -------------------------------------------------------------------------------------------------
{
    int r = (i1 - i0 + 1) % 2;
    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_ouverture3_ui8matrix_fusion_ilu5_elu2_red(X,i,j0,j1,Y);
    }
    if (r) ouverture3_ui8matrix_fusion_ilu5_red(X,i1-r+1,i1,j0,j1,Y);
}
// --------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion_ilu5_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// --------------------------------------------------------------------------------------------------------
{
    int r = (i1 - i0 + 1) % 2;
    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_ouverture3_ui8matrix_fusion_ilu5_elu2_red_factor(X,i,j0,j1,Y);
    }
    if (r) ouverture3_ui8matrix_fusion_ilu5_red(X,i1-r+1,i1,j0,j1,Y);
}
// ---------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_fusion_ilu15_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// ---------------------------------------------------------------------------------------------
{
        for (int i = i0 ; i <= i1 ; i++) {
        line_ouverture3_ui8matrix_fusion_ilu15_red(X,i,j0,j1,Y);
    }
}
// ------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_basic(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ------------------------------------------------------------------------------------------------------
{
    //Prologue
    line_min3_ui8matrix_basic(X,i0 - 1,j0,j1,T);
    line_min3_ui8matrix_basic(X,i0,j0,j1,T);

    for (int i = i0 ; i <= i1 ; i++) {
        line_min3_ui8matrix_basic(X,i+1,j0,j1,T);  
        line_max3_ui8matrix_basic(T,i,j0,j1,Y);
    }
}
// ----------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ----------------------------------------------------------------------------------------------------
{
    //Prologue
    line_min3_ui8matrix_red(X,i0 - 1,j0,j1,T);
    line_min3_ui8matrix_red(X,i0,j0,j1,T);

    for (int i = i0 ; i <= i1 ; i++) {
        line_min3_ui8matrix_red(X,i+1,j0,j1,T);  
        line_max3_ui8matrix_red(T,i,j0,j1,Y);
    }
}
// ---------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_ilu3_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ---------------------------------------------------------------------------------------------------------
{
    //Prologue
    line_min3_ui8matrix_ilu3_red(X,i0 - 1,j0,j1,T);
    line_min3_ui8matrix_ilu3_red(X,i0,j0,j1,T);

    for (int i = i0 ; i <= i1 ; i++) {
        line_min3_ui8matrix_ilu3_red(X,i+1,j0,j1,T);  
        line_max3_ui8matrix_ilu3_red(T,i,j0,j1,Y);
    }
}
// ---------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ---------------------------------------------------------------------------------------------------------
{
    //Prologue
    int r = (i1 - i0 + 1) % 2;
    line_min3_ui8matrix_elu2_red(X,i0 - 1,j0,j1,T);


    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_min3_ui8matrix_elu2_red(X,i+1,j0,j1,T);  
        line_max3_ui8matrix_elu2_red(T,i,j0,j1,Y);
    }

    if (r) {
        line_max3_ui8matrix_red(T, i1-r+1, j0, j1, Y);
    }
}
// ----------------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ----------------------------------------------------------------------------------------------------------------
{
    //Prologue
    int r = (i1 - i0 + 1) % 2;
    line_min3_ui8matrix_elu2_red_factor(X,i0 - 1,j0,j1,T);


    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_min3_ui8matrix_elu2_red_factor(X,i+1,j0,j1,T);  
        line_max3_ui8matrix_elu2_red_factor(T,i,j0,j1,Y);
    }

    if (r) {
        line_max3_ui8matrix_red(T, i1-r+1, j0, j1, Y);
    }
}
// --------------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_ilu3_elu2_red(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// --------------------------------------------------------------------------------------------------------------
{
    //Prologue
    int r = (i1 - i0 + 1) % 2;
    line_min3_ui8matrix_ilu3_elu2_red(X,i0 - 1,j0,j1,T);


    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_min3_ui8matrix_ilu3_elu2_red(X,i+1,j0,j1,T);  
        line_max3_ui8matrix_ilu3_elu2_red(T,i,j0,j1,Y);
    }

    if (r) {
        line_max3_ui8matrix_ilu3_red(T, i1-r+1, j0, j1, Y);
    }
}
// ---------------------------------------------------------------------------------------------------------------------
void ouverture3_ui8matrix_pipeline_ilu3_elu2_red_factor(uint8 **X, int i0, int i1, int j0, int j1, uint8 **T, uint8 **Y)
// ---------------------------------------------------------------------------------------------------------------------
{
    //Prologue
    int r = (i1 - i0 + 1) % 2;
    line_min3_ui8matrix_ilu3_elu2_red_factor(X,i0 - 1,j0,j1,T);


    for (int i = i0 ; i <= i1-r ; i+=2) {
        line_min3_ui8matrix_ilu3_elu2_red_factor(X,i+1,j0,j1,T);  
        line_max3_ui8matrix_ilu3_elu2_red_factor(T,i,j0,j1,Y);
    }

    if (r) {
        line_max3_ui8matrix_ilu3_red(T, i1-r+1, j0, j1, Y);
    }
}
