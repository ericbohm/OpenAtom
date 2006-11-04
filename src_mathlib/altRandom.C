#define MODULUS_R    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER_R 48271      /* DON'T CHANGE THIS VALUE                  */

//=================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================
// Random returns a pseudo-random real number uniformly distributed 
// between 0.0 and 1.0. 
//=================================================================
   double altRandom(long *seed){
    long t;
    const long Q = MODULUS_R / MULTIPLIER_R;
    const long R = MODULUS_R % MULTIPLIER_R;

    t = MULTIPLIER_R * (seed[0] % Q) - R * (seed[0] / Q);
    if(t > 0){
      seed[0] = t;
    }else {
      seed[0] = t + MODULUS_R;
    }//endif
    return ((double) seed[0] / MODULUS_R);
}
//=================================================================
