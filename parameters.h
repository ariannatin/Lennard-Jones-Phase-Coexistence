// the simulation is runned at fixed temperature, chemical potential and volume

// parameters of the simulation
#define T 1.1                                             // temperature
#define zeta 0.0473                                       // fugacity
#define N_SAMPLE_MAX 500                                  // maximum value of N to sample with SUS
#define nSteps 200                                       // number of Monte Carlo steps
#define W 1                                               // dimension of the window in Successive Umbrella Sampling


#define L 9.0755                                          // size of the simulation box
#define V (L*L*L)                                         // volume
#define k_B 1.0                                           // Boltzmann constant
#define beta (1/T)                                        // beta
#define r_cut 2.5                                         // cut-off of Lennard-Jones potential
#define r2_cut (r_cut*r_cut)                              // squared cut-off
#define r6_cut (r_cut*r_cut*r_cut*r_cut*r_cut*r_cut)      // sixth power of cut-off
#define Li (1./L)                                         // inverse of the side size
#define L_2 (L*L)                                         // squared side size
#define doubleL (2*L)                                     // twice the side size
#define NMAX 1000                                         // maximum number of particles
#define delta 0.35                                        // magnitude of displacement
#define lambda 1                                          // thermal De Broglie wavelength
#define zz (zeta/(lambda*lambda*lambda))                  // activity
#define nDispl 1000                                       // number of trial displacement for every step
#define nEx 10                                            // number of trial exchange of particles
#define N_WINDOWS (N_SAMPLE_MAX/W)                        // number of windows for SUS
#define INIT_FROM_FILE 0                                  // decide to begin from the last simulation

// values for argon
#define MASS 6.6335209e-26
#define EPSILON 1.65e-21
#define SIGMA 3.405e-10
#define TT (T*EPSILON/k_b)                                // absolute temperature
#define TAU (sqrt(MASS/EPSILON)*SIGMA)                    // time unity

#define h_bar (1.054571817e-34/TAU/EPSILON)               // reduced Planck constant in dimensionless units
