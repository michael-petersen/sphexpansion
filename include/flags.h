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
#define HAVEYAML 1
#endif

// flags for options in different headers

// additional diagnostic output in coefficient files (sphcoefs, cylcoefs)
#ifndef DEBUGCOEFS
#define DEBUGCOEFS 1
#endif

// flags for deep debugging (expansion.h)
#ifndef DEEPDEBUGCOEFS
#define DEEPDEBUGCOEFS 1
#endif

#ifndef DEEPDEBUGTIME
#define DEEPDEBUGTIME 1
#endif

#ifndef DEBUGCACHE
#define DEBUGCACHE 1
#endif

// diagnostics for debugging cylexpansion.h
#ifndef FORMAT
#define FORMAT 1
#endif

#ifndef ORDERS
#define ORDERS 1
#endif


#endif
