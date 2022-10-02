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

//#include "parser.h"
//#include "sequence.h"
//#include "sigmadelta.h"

//#include "features.h"
//#include "ecc_filter.h"
//#include "morpho.h"

//#include "nrutil_ext_test.h"
#include "swp_test.h"

#include "morpho_test.h"
#include "motion_test.h"

/* --------------- */
int getIter(int size)
/* --------------- */
{
    if(size <   32) return 10240;
    if(size <   64) return 2560 ;
    if(size <  128) return 640;
    if(size <  256) return 160;
    if(size <  512) return 40;
    if(size < 1024) return 20;
    if(size < 2048) return 10;
    return 2;
}
/* --------------- */
int getIterAV(int size)
/* --------------- */
{
    return 2*getIter(size);
}
/*
#include <stdint.h>
#include <string.h> // memset
#include <mach/mach_time.h>
double chrono(void)
{
static int first_time = 1;
static double conversion = 0.0;

static mach_timebase_info_data_t timebase;
static kern_return_t err;

uint64_t t = mach_absolute_time();

if(first_time) {
memset( &timebase, 0, sizeof( timebase ));
err = mach_timebase_info( &timebase );
//if( 0 != err ) return err;
if(err != 0) return 0.0;
conversion = 1e-9 * (double) timebase.numer / (double) timebase.denom;
//printf("numer = %d denom = %d\n", timebase.numer, timebase.denom);
first_time = 0;
}
return conversion * (double) t;
}
*/
// =============================
int main(int argc, char *argv[])
// =============================
{
    // test_swp(argc, argv);
    
    test_morpho(argc, argv);
    test_motion(argc, argv);
    return 0;
}
