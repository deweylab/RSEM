/*************************** RANDOMC.H ***************** 2007-09-22 Agner Fog *
*
* This file contains class declarations and other definitions for the C++ 
* library of uniform random number generators.
*
* Overview of classes:
* ====================
*
* class CRandomMersenne:
* Random number generator of type Mersenne twister.
* Source file mersenne.cpp
*
* class CRandomMother:
* Random number generator of type Mother-of-All (Multiply with carry).
* Source file mother.cpp
*
*
* Member functions (methods):
* ===========================
*
* All these classes have identical member functions:
*
* Constructor(uint32 seed):
* The seed can be any integer. Usually the time is used as seed.
* Executing a program twice with the same seed will give the same sequence of
* random numbers. A different seed will give a different sequence.
*
* void RandomInit(uint32 seed);
* Re-initializes the random number generator with a new seed.
*
* void RandomInitByArray(uint32 seeds[], int length);
* In CRandomMersenne only: Use this function if you want to initialize with
* a seed with more than 32 bits. All bits in the seeds[] array will influence
* the sequence of random numbers generated. length is the number of entries
* in the seeds[] array.
*
* double Random();
* Gives a floating point random number in the interval 0 <= x < 1.
* The resolution is 32 bits in CRandomMother and CRandomMersenne.
*
* int IRandom(int min, int max);
* Gives an integer random number in the interval min <= x <= max.
* (max-min < MAXINT).
* The precision is 2^-32 (defined as the difference in frequency between 
* possible output values). The frequencies are exact if max-min+1 is a
* power of 2.
*
* int IRandomX(int min, int max);
* Same as IRandom, but exact. In CRandomMersenne only.
* The frequencies of all output values are exactly the same for an 
* infinitely long sequence. (Only relevant for extremely long sequences).
*
* uint32 BRandom();
* Gives 32 random bits. 
*
*
* Example:
* ========
* The file EX-RAN.CPP contains an example of how to generate random numbers.
*
*
* Optimized version:
* ==================
* Faster versions of these random number generators are provided as function
* libraries in asmlib.zip. These function libraries are coded in assembly
* language and support only x86 platforms, including 32-bit and 64-bit
* Windows, Linux, BSD, Mac OS-X (Intel based). Use asmlibran.h from asmlib.zip
*
*
* Non-uniform random number generators:
* =====================================
* Random number generators with various non-uniform distributions are available
* in stocc.zip (www.agner.org/random).
*
*
* Further documentation:
* ======================
* The file randomc.htm contains further documentation on these random number
* generators.
*
*
* Copyright:
============
* © 1997 - 2007 Agner Fog. All software in this library is published under the
* GNU General Public License with the further restriction that it cannot be 
* used for gambling applications. See licence.htm
*******************************************************************************/

#ifndef RANDOMC_H
#define RANDOMC_H


// Define 32 bit signed and unsigned integers.
// Change these definitions, if necessary, to match a particular platform
#if defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
   // 16 bit systems use long int for 32 bit integer
   typedef long int           int32;   // 32 bit signed integer
   typedef unsigned long int  uint32;  // 32 bit unsigned integer
#else
   // Most other systems use int for 32 bit integer
   typedef int                int32;   // 32 bit signed integer
   typedef unsigned int       uint32;  // 32 bit unsigned integer
#endif

// Define 64 bit signed and unsigned integers, if possible
#if (defined(__WINDOWS__) || defined(_WIN32)) && (defined(_MSC_VER) || defined(__INTEL_COMPILER))
   // Microsoft and other compilers under Windows use __int64
   typedef __int64            int64;   // 64 bit signed integer
   typedef unsigned __int64   uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#elif defined(__unix__) && (defined(_M_IX86) || defined(_M_X64))
   // Gnu and other compilers under Linux etc. use long long
   typedef long long          int64;   // 64 bit signed integer
   typedef unsigned long long uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#else
   // 64 bit integers not defined
   // You may include definitions for other platforms here
#endif


/***********************************************************************
System-specific user interface functions
***********************************************************************/

void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

void FatalError(char * ErrorText);     // System-specific error reporting (userintf.cpp)


/***********************************************************************
Define random number generator classes
***********************************************************************/

class CRandomMersenne {                // Encapsulate random number generator
#if 0
   // Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else    
   // or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
public:
   CRandomMersenne(uint32 seed) {      // Constructor
      RandomInit(seed); LastInterval = 0;}
   void RandomInit(uint32 seed);       // Re-seed
   void RandomInitByArray(uint32 seeds[], int length); // Seed by more than 32 bits
   int IRandom (int min, int max);     // Output random integer
   int IRandomX(int min, int max);     // Output random integer, exact
   double Random();                    // Output random float
   uint32 BRandom();                   // Output random bits
private:
   void Init0(uint32 seed);            // Basic initialization procedure
   uint32 mt[MERS_N];                  // State vector
   int mti;                            // Index into mt
   uint32 LastInterval;                // Last interval length for IRandomX
   uint32 RLimit;                      // Rejection limit used by IRandomX
   enum TArch {LITTLE_ENDIAN1, BIG_ENDIAN1, NONIEEE}; // Definition of architecture
   TArch Architecture;                 // Conversion to float depends on architecture
};    


class CRandomMother {             // Encapsulate random number generator
public:
   void RandomInit(uint32 seed);       // Initialization
   int IRandom(int min, int max);      // Get integer random number in desired interval
   double Random();                    // Get floating point random number
   uint32 BRandom();                   // Output random bits
   CRandomMother(uint32 seed) {   // Constructor
      RandomInit(seed);}
protected:
   uint32 x[5];                        // History buffer
};

#endif
