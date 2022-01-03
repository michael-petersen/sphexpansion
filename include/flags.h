/*
preprocessor flag list for expansions

MSP 24 Dec 2021 original commit

*/
#ifndef FLAGS_H
#define FLAGS_H


// use splines in orient files (orient.h)
#define SPLINEORIENT 0

// use splines in coefficient files (sphcoefs)
#define SPLINECOEFS 0

// use splines in model files (sphmodel)
#define SPLINEMODEL 1

// additional diagnostic output in coefficient files (sphcoefs, cylcoefs)
#define DEBUGCOEFS 0

// flags for deep debugging (expansion)
#define DEEPDEBUGCOEFS 0
#define DEEPDEBUGTIME 0

#endif

