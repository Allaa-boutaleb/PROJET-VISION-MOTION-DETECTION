/* ---------------- */
/* --- motion.c --- */
/* ---------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "nrtype.h"
#include "nrdef.h"
#include "nrutil.h"

#include "sigmadelta.h"
#include "morpho_min.h"
#include "morpho_max.h"
#include "morpho_ouverture.h"
#include "motion.h"

// -------------------------------------------------------
void init_image(uint8 **X, int i0, int i1, int j0, int j1)
// -------------------------------------------------------
{
    // for debugging
    for(int i = i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            X[i][j] = i+j;
        }
    }
}
// ----------------------------------------------------------
void reverse_video(uint8 **X, int i0, int i1, int j0, int j1)
// ----------------------------------------------------------
{
    // for debugging
    for(int i = i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            X[i][j] = 255 - X[i][j];
        }
    }
}
// --------------------------------------------------------------------------------------------------
void threshold_ui8matrix(uint8 **X, uint8 t, uint8 valeur, uint8 **Y, int i0, int i1, int j0, int j1)
// --------------------------------------------------------------------------------------------------
{
    for(int i=i0; i<=i1; i++) {
        for(int j=j0; j<=j1; j++) {
            
            if(X[i][j] < t)
                Y[i][j] = 0;
            else
                Y[i][j] = valeur;
        }
    }
}
// ----------------
void test_PGM(void)
// ----------------
{
    char *src_path = SEQUENCE_SRC_PATH;
    char *dst_path = SEQUENCE_DST_PATH;
    char *filename = SEQUENCE_FILENAME;
    int ndigit     = SEQUENCE_NDIGIT;
    
    char complete_filenameI[1024];
    char complete_filenameM[1024];
    
    int tstart = SEQUENCE_TSTART;
    int tstop  = SEQUENCE_TSTOP;
    int tstep  = SEQUENCE_TSTEP;
      
    int h = SEQUENCE_HEIGHT;
    int w = SEQUENCE_WIDTH;
    
    int i0, i1, j0, j1;
    long li0, li1, lj0, lj1;
    
    uint8 **I, **M;
    
    puts("--------------");
    puts("-- test_PGM --");
    puts("--------------");
    
    i0 = 0; i1 = h-1; j0 = 0; j1 = w-1;
    
    src_path = "C:/Users/Studies/Desktop/Vision/LO3IN401_projet_etu/projet_etu/results/";
    
    generate_path_filename_k_ndigit_extension(dst_path, "I_", tstart, ndigit, "pgm", complete_filenameI);
    generate_path_filename_k_ndigit_extension(dst_path, "M_", tstart, ndigit, "pgm", complete_filenameM);
    
    puts(complete_filenameI);
    puts(complete_filenameM);
    
    I  = ui8matrix(0, h-1, 0, w-1);
    M  = ui8matrix(0, h-1, 0, w-1);
    
    init_image(I, 0, h-1, 0, w-1);
    
    SavePGM_ui8matrix(I, i0, i1, j0, j1, complete_filenameI);
    MLoadPGM_ui8matrix(complete_filenameI, i0, i1, j0, j1, M);
    reverse_video(M, i0+h/4, i1-h/4, j0+w/8, j1-w/8);
    SavePGM_ui8matrix(M, i0, i1, j0, j1, complete_filenameM);
    
    free_ui8matrix(I, 0, h-1, 0, w-1);
    free_ui8matrix(M, 0, h-1, 0, w-1);
}
// ==================================
void motion_detection_morpho_v1(void)
// ==================================
{
    // version sans optimisation servant de base apres validation des operateurs morphologiques sans optimisation
    
    // sequence
    char *src_path;
    char *dst_path;
    char *filename;
    int ndigit;
    int tstart, tstop, tstep;
    
    // image
    int h, w, b; // height, width, border
    int i0, i1, j0, j1;
    long li0, li1, lj0, lj1;
    
    // sigma-delta
    uint8 **I, **M, **O, **V, **E, **E_8;
    
    // morpho
    uint8 **Erosion1, **Dilatation1, **Dilatation2, **Erosion2;
    uint8 **Erosion1_8, **Dilatation1_8, **Dilatation2_8, **Erosion2_8;
    // filename for saving data
    
    char complete_filename_I[1024];
    char complete_filename_M[1024];
    char complete_filename_O[1024];
    char complete_filename_V[1024];
    char complete_filename_E[1024];
    
    char complete_filename_erosion1[1024];
    char complete_filename_erosion2[1024];
    char complete_filename_dilatation1[1024];
    char complete_filename_dilatation2[1024];
    
    // -------
    // -- init
    // -------
    //puts("[motion_detection_morpho_v1]: les parametres de la sequence sont dans motion.h");
    //puts("[motion_detection_morpho_v1]: supprimer ce message une fois les parametres configures");
    //puts("[motion_detection_morpho_v1]: bye...");
    //return;
    
    src_path = SEQUENCE_SRC_PATH;
    filename = SEQUENCE_FILENAME;
    ndigit   = SEQUENCE_NDIGIT;
    
    tstart = SEQUENCE_TSTART;
    tstop  = SEQUENCE_TSTOP;
    tstep  = SEQUENCE_TSTEP;
    
    h = SEQUENCE_HEIGHT;
    w = SEQUENCE_WIDTH;
    b = 2;
    dst_path = SEQUENCE_DST_PATH;
    
    i0 = 0; i1 = h-1; j0 = 0; j1 = w-1;
    
    // ----------------
    // -- allocation --
    // ----------------
    
    I   = ui8matrix(i0,   i1,   j0,   j1);
    M   = ui8matrix(i0,   i1,   j0,   j1);
    O   = ui8matrix(i0,   i1,   j0,   j1);
    V   = ui8matrix(i0,   i1,   j0,   j1);
    E   = ui8matrix(i0-b, i1+b, j0-b, j1+b); // image 1 bit en binaire {0,1}
    E_8 = ui8matrix(i0,   i1,   j0,   j1); // image 8 bits en, niveau de gris {0,255}
    
    // image 1bit/pixel
    Erosion1    = ui8matrix(i0-b, i1+b, j0-b, j1+b);
    Erosion2    = ui8matrix(i0-b, i1+b, j0-b, j1+b);
    Dilatation1 = ui8matrix(i0-b, i1+b, j0-b, j1+b);
    Dilatation2 = ui8matrix(i0-b, i1+b, j0-b, j1+b);
    
    // image 8bits pour la visualisation
    Erosion1_8    = ui8matrix(i0, i1, j0, j1);
    Erosion2_8    = ui8matrix(i0, i1, j0, j1);
    Dilatation1_8 = ui8matrix(i0, i1, j0, j1);
    Dilatation2_8 = ui8matrix(i0, i1, j0, j1);
    // --------------------
    // -- initialisation --
    // --------------------
    
    // bords exterieurs a zero
    
    zero_ui8matrix(I  , i0,   i1,   j0,   j1);
    zero_ui8matrix(M  , i0,   i1,   j0,   j1);
    zero_ui8matrix(O  , i0,   i1,   j0,   j1);
    zero_ui8matrix(V  , i0,   i1,   j0,   j1);
    zero_ui8matrix(E  , i0-b, i1+b, j0-b, j1+b); // image 1 bit en binaire {0,1}
    zero_ui8matrix(E_8, i0,   i1,   j0,   j1); // image 8 bits en, niveau de gris {0,255}
    
    zero_ui8matrix(Erosion1   , i0-b, i1+b, j0-b, j1+b);
    zero_ui8matrix(Erosion2   , i0-b, i1+b, j0-b, j1+b);
    zero_ui8matrix(Dilatation1, i0-b, i1+b, j0-b, j1+b);
    zero_ui8matrix(Dilatation2, i0-b, i1+b, j0-b, j1+b);
    
    zero_ui8matrix(Erosion1_8   , i0, i1, j0, j1);
    zero_ui8matrix(Erosion2_8   , i0, i1, j0, j1);
    zero_ui8matrix(Dilatation1_8, i0, i1, j0, j1);
    zero_ui8matrix(Dilatation2_8, i0, i1, j0, j1);
    
    // ----------------
    // -- traitement --
    // ----------------
    
    // -- prologue --
    
    generate_path_filename_k_ndigit_extension(src_path, filename, tstart, ndigit, "pgm", complete_filename_I);
    //puts(complete_filenameI);
        
    MLoadPGM_ui8matrix(complete_filename_I, i0, i1, j0, j1, I);
    
    // fragment de code pour verifier qu'on lit bien les bonnes images
    //SavePGM_ui8matrix(I, i0, i1, j0, j1, "verif.pgm");
    //reverse_video(I, i0+h/4, i1-h/4, j0+w/8, j1-w/8);
    //SavePGM_ui8matrix(I, i0, i1, j0, j1, complete_filename_M);
    
    SigmaDelta_Step0(I, M, O, V, E, i0, i1, j0, j1);
    
    // -- boucle --
    for(int t=tstart; t<=tstop; t+=tstep) {
        
        printf("-- i = %3d ----------\n", t);
        
        generate_path_filename_k_ndigit_extension(src_path, filename, t, ndigit, "pgm", complete_filename_I);
        MLoadPGM_ui8matrix(complete_filename_I, i0, i1, j0, j1, I);

        // N = 3 ecart type autour de la moyenne
        SigmaDelta_1Step(I, M, O, V, E, 3, i0, i1, j0, j1);
    
        // morpho en niveau de gris fonctionnant aussi sur des images 1 bit / pixel
        min3_ui8matrix_basic(E,           i0, i1, j0, j1, Erosion1);


        // HELLOOOOOOOOOOOOOOO
        min3_ui8matrix_ilu3(E,           i0, i1, j0, j1, Erosion1);
        max3_ui8matrix_ilu3(Erosion1,    i0, i1, j0, j1, Dilatation1);

        // ouverture3_ui8matrix_pipeline_ilu3_elu2_red(E, i0, i1, j0, j1,Erosion1, Dilatation1);
        max3_ui8matrix_ilu3(Dilatation1, i0, i1, j0, j1, Dilatation2);
        min3_ui8matrix_ilu3(Dilatation2, i0, i1, j0, j1, Erosion2);
    
        // traitement pour visualisation
        threshold_ui8matrix(E,           1, 255, E_8,           i0, i1, j0, j1);
        threshold_ui8matrix(Erosion1,    1, 255, Erosion1_8,    i0, i1, j0, j1);
        threshold_ui8matrix(Erosion2,    1, 255, Erosion2_8,    i0, i1, j0, j1);
        threshold_ui8matrix(Dilatation1, 1, 255, Dilatation1_8, i0, i1, j0, j1);
        threshold_ui8matrix(Dilatation2, 1, 255, Dilatation2_8, i0, i1, j0, j1);
        
        generate_path_filename_k_ndigit_extension(dst_path, "I_",        t, ndigit, "pgm", complete_filename_I);
        generate_path_filename_k_ndigit_extension(dst_path, "M_",        t, ndigit, "pgm", complete_filename_M);
        generate_path_filename_k_ndigit_extension(dst_path, "O_",        t, ndigit, "pgm", complete_filename_O);
        generate_path_filename_k_ndigit_extension(dst_path, "V_",        t, ndigit, "pgm", complete_filename_V);
        generate_path_filename_k_ndigit_extension(dst_path, "E_",        t, ndigit, "pgm", complete_filename_E);
        generate_path_filename_k_ndigit_extension(dst_path, "Ero1_",     t, ndigit, "pgm", complete_filename_erosion1);
        generate_path_filename_k_ndigit_extension(dst_path, "Ero2_",     t, ndigit, "pgm", complete_filename_erosion2);
        generate_path_filename_k_ndigit_extension(dst_path, "Dil1_",     t, ndigit, "pgm", complete_filename_dilatation1);
        generate_path_filename_k_ndigit_extension(dst_path, "Dil2_",     t, ndigit, "pgm", complete_filename_dilatation2);
        
        //puts(complete_filename_I);
        //puts(complete_filename_M);
        //puts(complete_filename_O);
        //puts(complete_filename_V);
        //puts(complete_filename_E);
        
        puts(complete_filename_erosion1);
        puts(complete_filename_dilatation1);
        puts(complete_filename_dilatation2);
        puts(complete_filename_erosion2);
        
        SavePGM_ui8matrix(I, i0, i1, j0, j1, complete_filename_I);
        SavePGM_ui8matrix(M, i0, i1, j0, j1, complete_filename_M);
        SavePGM_ui8matrix(O, i0, i1, j0, j1, complete_filename_O);
        SavePGM_ui8matrix(V, i0, i1, j0, j1, complete_filename_V);
        SavePGM_ui8matrix(E, i0, i1, j0, j1, complete_filename_E);
        
        SavePGM_ui8matrix(Erosion1_8   , i0, i1, j0, j1, complete_filename_erosion1   );
        SavePGM_ui8matrix(Dilatation1_8, i0, i1, j0, j1, complete_filename_dilatation1);
        SavePGM_ui8matrix(Dilatation2_8, i0, i1, j0, j1, complete_filename_dilatation2);
        SavePGM_ui8matrix(Erosion2_8   , i0, i1, j0, j1, complete_filename_erosion2   );
    } // t
    
    // ----------
    // -- free --
    // ----------
    
    free_ui8matrix(I  , i0,   i1,   j0,   j1);
    free_ui8matrix(M  , i0,   i1,   j0,   j1);
    free_ui8matrix(O  , i0,   i1,   j0,   j1);
    free_ui8matrix(V  , i0,   i1,   j0,   j1);
    free_ui8matrix(E  , i0-1, i1+1, j0-1, j1+1);
    free_ui8matrix(E_8, i0,   i1,   j0,   j1);
    
    free_ui8matrix(Erosion1   , i0-b, i1+b, j0-b, j1+b);
    free_ui8matrix(Erosion2   , i0-b, i1+b, j0-b, j1+b);
    free_ui8matrix(Dilatation1, i0-b, i1+b, j0-b, j1+b);
    free_ui8matrix(Dilatation2, i0-b, i1+b, j0-b, j1+b);
    
    free_ui8matrix(Erosion1_8   , i0, i1, j0, j1);
    free_ui8matrix(Erosion2_8   , i0, i1, j0, j1);
    free_ui8matrix(Dilatation1_8, i0, i1, j0, j1);
    free_ui8matrix(Dilatation2_8, i0, i1, j0, j1);
}
// ===============================
void motion_detection_morpho(void)
// ===============================
{
    motion_detection_morpho_v1(); // version basique sans optimisation
    test_PGM();
}
