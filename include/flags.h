/*
preprocessor flag list for expansions

MSP 24 Dec 2021 original commit
MSP 11 Jan 2022 set eigen flag
MSP 14 Jan 2022 set yaml flag

*/
#ifndef FLAGS_H
#define FLAGS_H

// flags for external libraries

// eigen compilation
#ifndef HAVEEIGEN
#define HAVEEIGEN 0
#endif

// yaml compilation
#ifndef HAVEYAML
#define HAVEYAML 0
#endif

// flags for options in different headers

// use splines in orient files (orient.h)
#ifndef SPLINEORIENT
#define SPLINEORIENT 0
#endif

// use splines in coefficient files (sphcoefs.h)
#ifndef SPLINECOEFS
#define SPLINECOEFS 0
#endif

// use splines in model files (sphmodel)
#ifndef SPLINEMODEL
#define SPLINEMODEL 1
#endif

// additional diagnostic output in coefficient files (sphcoefs, cylcoefs)
#ifndef DEBUGCOEFS
#define DEBUGCOEFS 0
#endif

// flags for deep debugging (expansion.h)
#ifndef DEEPDEBUGCOEFS
#define DEEPDEBUGCOEFS 0
#endif

#ifndef DEEPDEBUGTIME
#define DEEPDEBUGTIME 0
#endif

// diagnostics for debugging cylexpansion.h
#ifndef FORMAT
#define FORMAT 0
#endif

#ifndef ORDERS
#define ORDERS 0
#endif


#endif
