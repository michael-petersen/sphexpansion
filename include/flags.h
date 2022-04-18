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
#define HAVEEIGEN 1
#endif

// yaml compilation
#ifndef HAVEYAML
#define HAVEYAML 0
#endif

// flags for options in different headers

// additional diagnostic output in coefficient files (sphcoefs.h, cylcoefs.h)
#ifndef DEBUGCOEFS
#define DEBUGCOEFS 0
#endif

// flags for deep debugging (sphexpansion.h, cylexpansion.h)
#ifndef DEEPDEBUGCOEFS
#define DEEPDEBUGCOEFS 0
#endif

#ifndef DEEPDEBUGTIME
#define DEEPDEBUGTIME 0
#endif

#ifndef DEBUGCACHE
#define DEBUGCACHE 0
#endif

// diagnostics for debugging cylexpansion.h
#ifndef FORMAT
#define FORMAT 0
#endif

#ifndef ORDERS
#define ORDERS 0
#endif


#endif
