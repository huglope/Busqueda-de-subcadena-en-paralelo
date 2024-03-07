/*
 * Simple random generator
 * LCG (Linear Congruential Generator)
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2023/2024
 *
 * v1.0
 *
 * (c) 2024, Arturo Gonzalez-Escribano
 */
#include<stdint.h>
#include<math.h>

/*
 * Constants
 */
#define RNG_MULTIPLIER 6364136223846793005ULL
#define RNG_INCREMENT  1442695040888963407ULL
#define RNG_MODULUS    18446744073709551615ULL // 2^64 - 1

/*
 * Type for random sequences state
 */
typedef uint64_t	rng_t;

/*
 * Constructor: Create a new state from a seed
 */
rng_t rng_new(uint64_t seed) {
    uint64_t hash = seed;
    hash = (hash ^ (hash >> 30)) * 0xbf58476d1ce4e5b9ULL;
    hash = (hash ^ (hash >> 27)) * 0x94d049bb133111ebULL;
    hash = hash ^ (hash >> 31);
    return hash; // initial state
}

/*
 * Next: Advance state and return a double number uniformely distributed
 */
double rng_next(rng_t *seq) {
    *seq = ( *seq * RNG_MULTIPLIER + RNG_INCREMENT) % RNG_MODULUS;
    return (double) *seq / RNG_MODULUS;
}

/*
 * Next Normal: Advance state and return a double number distributed with a normal(mu,sigma)
 */
double rng_next_normal( rng_t *seq, double mu, double sigma) {
    double u1 = rng_next(seq);
    double u2 = rng_next(seq);

    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    // double z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
    
    return mu + sigma * z0;
}
