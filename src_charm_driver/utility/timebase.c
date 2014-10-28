
#define SPRN_TBRL       0x10C       // Time Base Read Lower Register (user & sup R/O)
#define SPRN_TBRU       0x10D       // Time Base Read Upper Register (user & sup R/O)

#define _bgp_mfspr( SPRN )\
  ({\
   unsigned int tmp;\
   do {\
   asm volatile ("mfspr %0,%1" : "=&r" (tmp) : "i" (SPRN) : "memory" );\
   }\
   while(0);\
   tmp;\
   })\

unsigned long long timebase( void )
{
  union {
    unsigned int ul[2];
    unsigned long long  ull;
  }
  hack;
  unsigned int utmp;

  do {
    utmp       = _bgp_mfspr( SPRN_TBRU );
    hack.ul[1] = _bgp_mfspr( SPRN_TBRL );
    hack.ul[0] = _bgp_mfspr( SPRN_TBRU );
  }
  while( utmp != hack.ul[0] );

  return( hack.ull );
}

double getTime( void )
{
  union {
    unsigned int ul[2];
    unsigned long long  ull;
  }
  hack;
  unsigned int utmp;

  do {
    utmp       = _bgp_mfspr( SPRN_TBRU );
    hack.ul[1] = _bgp_mfspr( SPRN_TBRL );
    hack.ul[0] = _bgp_mfspr( SPRN_TBRU );
  }
  while( utmp != hack.ul[0] );

  return( (1.0/850000000.0) * ((double) hack.ull)  );
}

