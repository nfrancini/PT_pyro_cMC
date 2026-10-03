// Librerie
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>
#include <errno.h>
#include <fenv.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <omp.h>
#include "../dSFMT/dSFMT.h"
#include "../dSFMT/dSFMT-jump.h"

// Macro
#define STD_STRING_LENGTH 50    // Lunghezza standard per le stringhe
#define PI 3.14159265358979     // Pi greco
#define NUM_COMPS 3             // Numero di componenti dello spin
#define NUM_NEAREST_NEIGHBOR 6  // Numero di primi vicini per ogni spin
#define NUM_SUBLATTICES 4       // Number of sublattices
#define INORM 10000             // Numero di misure tra una rinormalizzazione e l'altra
#define IOVERR 3
// #define DEBUG

// Variabili globali
// jump per il generatore random
static const char * jump = "f4dfa6c62049d0776e0bf6f1e953f3aa38abb113df86be024eab3773ad5f2b82ead936022e656dff7e562691c59dd5f7d2566b78d9669002503c4ddb1888a49f32333f515e6c60c4ecd221078ec6f26f0a90f4875067ca1f399a99775037adf905566e2c7e6b42131420f8f04f112c92621c9b1502f2a8aefad6c"
    "667904af62f0d55e02d396902d3b89450103c5ce5fe0408d97cbb864861b49e4e42048ff3310b48faac55095a7f422eea4aade752f947f947c6be0a0c665bdea099246ab9eff658ea8ca468bf49d0227748367878de06d7bd86ea6708fcac6e252f5f00f04309b2aac3036b64afb39d990427c6c9f03477cc7e935c43c0e61bc161db8eb15516"
    "eee8cb377ecbc1849207990fb6778721b29bfe0d89bfda1b3772fa5b0b1f7ec3daf36052032285898c6f6396f55010c31f8201b7e2e51d94f920bfe57684c5415cc342cb39a0045d9793d13cf8646096daeb8bb9bfc20a90de8f2426da8733267a9b9674f32154e8f84a9932223a2ca3c787d0b66df6675febbdfcba2f9cef09c621c57e11098b"
    "3289c77397aaae8b104642ffe0c4b75598efbc53745984d68b4d6656cae299ae2be55217a9a02b009ca7be32f47fbe434bce4914a34d0c9b0085bede9b8a99319c34660d66f0124b5a7714c4bf3cbfec3ee43ed817087168bad80133bebaeeb68cf7929a24d1bb3de831a8340d220906ab04159cf94b21d5ee813bd7c80f10f01b43052af5309"
    "17513b169254c25d6fcfe6cb420d6ce92f54886ef6eaf9a5ba35e893ff593834d05ddf28899e42d729c7df3d21ef036020789739366f0c11ec52ff92a0bfd8ba69508e27b20fabb8217bd36b90e5aa918159ac87913bc7b46c04e366c23c92807fbe9c6a407e6a4db0b4fc23c3b6c706b5ca058fe8c190f849f18d16d6b48b5ed760eb202fd566"
    "291a799420b9654e08b8118bcbfead8e9dd2fdb9b053e9bdfb665285c78718f726d0b3d6c37e116428ec9ac9db2637259e4e8d6402bbada46c6bdb03985e19a82e9b4e57de1b025a3cb1f850beae7e8da9941655825bce0e89d536b6ee9064865b1a85c185e9fc9cb7f435de13d44773c00eed442a286e4ab807e3cab4dc3441d1b7d2af693812a"
    "e8b39652bb8c835fc895d13d6da93541afeadeee450475c29f3b2dfa8ef1c1e2547463b2cc2f0ff7a42ac4dd35e25c4fa030d2d2766fbe9f2d04c1304671747bace2f7dd55142bfa60f8cbc968bfc3d7a342152dc684a0fb5a32c0962a62b5220ac0f72add9d8b84d6cc76b97d03245e01fc8da3414a49bb4075d3488f29b56dc42ba69e3b58529"
    "448c943ecfd98b3784a39d0b8609a8fb945e757f4569f53bd2cf80f7f638acf5b67fe9c560a3b7b0cf7e0398f31aa8b03cf9c62b24296b6d8596b694469a02686c38daa16a1ef86e012d61a2f7de1693a5c00b3685175caec3c67146477eba54830f1d546cb18a553779aa46adb4f2010e33f3def847c7d89b51a8462b227605f6c920fd558a6da"
    "f64bc98682e508ae960c0c571870e603ba1fce0c13d53176f353fd319959e13db93eae1359f06e3dd4767c04f824cf34ec7bf8f60161ba1a615db82852eca9e3869afa711ab9a090660b0dc6cfbea310dda77e02310fbaeacd2636f975838c2dbcdbe9ac2cd85cee28f5e3f0c73abf62f9fa02cd79a7606b7ba855db68a07848b057c3aaf38f1a"
    "70086e14616f6f88305a1f9ce6b41378a620d4db3e0e7e1d421590dccaeff86212e232eeb5eb8a8d33a8c9b25ae88f3a7bd5032b4efa68f8af3186a02ffcbf5456f12beccace94c81c360cc4a0dcc642b59f991eec68c59af78139ca60b96d6a18e9535f8995e89bd2cf6a0aef3acffd33d1c0c1b79b66414a91d9f65b2b4ec65844b96f725d2"
    "b4b0c309f3eb9d714e9dd939bbdfd85ce8fb43679aeab13f6c29549949503c9466dbd337c4cdde46d6eacd15f21f4d8fdeaa627a47884c88a9c85f0b731d271a8ea7cb9e04a4a149c23c10f56b3a0476dc77a999d6e4f813e4b0f805e2a693e2ae4ae0ecc423c9ba5d17b42e691abf83784a582f2b1fd85d1e0a27ba38a500963568b2450363d"
    "2c5e3f7b8ba3e5b56e4e9f745a3a710bf2ae233c303068c532ce78ff031e6ab28b705dd94d7db4500909edb5626b8c9bd5ff4f0b4741388f0b91563ee516934c013e901572cba005ac5c535f4f107903be9af7b2793dfb61b5070facbe71eefe1b5600f975c8c38c3a2350d78beadfecb78e981164ae8bc866e732972d3ceef4aac68e15861"
    "f9b881d9b51b4edece150bc124b07645defb4202ef5d0e0962db98cae6ed459561c93c74c20bd64362e4f4fffc389a6cd80514604ff22eecc10c9cbc7981d19a8102b24146354c463107c9dc070e29e70df3578022acf72289ef071ab9f9402a544d0399f1b1e5f206b6d46d445f6d612a490e72918e00c853eda8493bef511149e80c9ab56"
    "e8b4b8cba3987249f77d060e61760e5792ac321c987c03c2606e9393a7970212992cdbd16448078d5039d4c2c3199714f53278f4f7b1d2e514cf95bdfc078b8bb0db659cb2c3f5cc02890ea84f05d414c88d2db9e9f8455659b9fa6254405317245fa070d6970cafb4dadb2522b490a5c8e02fe973a8cdbfbfbdbfb01535099ffba3d3896bc4"
    "d1189fc570c3e6fdc6469265b8da912772e75dd62ab71be507f700d56cac5e68fd6b57ec166168ab5258a69625c142a5b1b3519f94be1bde5e51d3bd8ea0c12d5af2fe4615b1b7bd4a96628a4fabc65925ff09718f63bbebaad98f89bd9543a27b3ff3b5d8bfa89f941a5eb8cc005ccd4a705190e1c9dc6a9f4264e5ee658520a4438e92de854"
    "bffc39f8dc7dfbb5de4f14ba63ea16a37d14a7b4610f95b6cffd55e4679b29cedbdf20e7bd16da822fad910c359ee3a68e48aae6e769b0e291d5d3aa3e2ca9d8d23abe8a1d5349f4991e9300852cc0befb20c2fc0d169306b260763344024f8092cbcc24c6807363e9fc548a30d5faab3a94b2af0782a2942be80c45d8b0587efd587394ef33c"
    "33022436e285806ddffdd32fe36345c3c38ed8d680abeb7a028b44ee6f94d060a14c7019bb6af1f1b5f0a562957d19826d8cc216f9b908c989ccd5415e3525dfe9422ffb5b50b7cc3083dc325544751e5683535d7439d3da2b0bb73bea551dd99e04e0e793804f4774eb6b1daf781d9caa5128274e599e847862fe309027813d3e4eda0bbeb72"
    "01856a5c5d8370e44dabff0bb229c723ba0a6bcf29c44536147de11b7835991018100105bd4329217f7386903fe8e7363cd7b3e893244e245e0a187467664c05b0be1fd429722b9b9a5e3198147fad72776e8a63aab9054fa9d259af0198d088d71d132e6068676a8e9ebb0f616b51ee34aac39c2c2221c71124017270d75ff4a048363c389e04"
    "e9b440ad2032a381ac2cfc54f409caa791e65ee4f5d6cd035008f219b88a803a7382ae447bf65a3df2176b25b3b7b67dabe34decd9a1384dc7a003916ca8fbcb29b3ad6fd8eac5bbbaa3bdfa6c6a3ad9427c4f3ed79fea26e14c8ce5fa3b4f82c5f7b6d2125916753a7b92ce9b46d45";

extern int (*nn)[NUM_NEAREST_NEIGHBOR];            // Variabile globale per i primi vicini
extern double (*rr)[NUM_NEAREST_NEIGHBOR][NUM_COMPS];

extern double acc_min;          // Accettanza minima per il metropolis
extern double acc_max;          // Accettanza massima per il metropolis

extern clock_t start, end;      // Variabili per tempo di esecuzione

extern double local_x[NUM_SUBLATTICES][NUM_COMPS];
extern double local_y[NUM_SUBLATTICES][NUM_COMPS];
extern double local_z[NUM_SUBLATTICES][NUM_COMPS];    // Variabile per salvare local_z axis

extern double DM_vectors[NUM_SUBLATTICES][NUM_SUBLATTICES][NUM_COMPS];

extern double vv[NUM_SUBLATTICES][NUM_COMPS][NUM_COMPS][NUM_COMPS];

// Structures
typedef struct {                // Structure for system parameter
  int size;                     // System size
  int N;                        // Number of sites

  int iStart;                   // Starting index: random (0), ordered (1), da file precedente (2)

  int n_temp;                   // Total temperatures
  double t_min;                 // Minimum temperature
  double t_max;                 // Maximum temperature

  double heis;                  // Heisenberg coupling
  double ising;                 // Kitaev coupling
  double PD;                    // Pseudodipolar coupling
  double DM;                    // DM coupling
  double single_ion;            // Single ion anisotropy

  int iSet;                     // Number of iterations to set the temperature
  int iPT;                      // Swaps number in the temperature setting phase for each iteration

  int iTerm;                    // thermalization loops, burn-in sweeps
                                // Used also to fix the Metropolis acceptance
                                // after the temperature setting

  int iMis;                     // Number of measures
  int iDec;                     // Decorrelation number of sweeps between two measures

  int iMinimize;                // Minimization index, if 1 iterative minimization is executed

  int iSave;                    // Every iSave append spin configs on file for each T
  int iBackup;                  // Every iBackup configs backup

  unsigned int dSFMT_seed;      // Seed for random generator

  char filenames_file[STD_STRING_LENGTH];   // Nome del file su cui salvo i nomi dei file
  char log_file[STD_STRING_LENGTH];         // Nome del file di log
} param_t;

typedef struct {                // Struttura della singola replica
  double (*lattice)[NUM_COMPS]; // Reticolo lattice[N][K], il primo indice è un indice di sito, il secondo di componente di spin
  double T;                     // Temperatura della replica
  unsigned long long acc;       // Numero di step Metropolis accettati
  unsigned long long count;     // Numero di step Metropolis provati
  double acc_Parameter;         // Parametro regolatore dell'accettanza Metropolis
} chain_t;

typedef struct {                // Struttura delle osservabili
  double ene_density;      // Densità di energia ene_density = ene_tot/N_siti
  double mA2;                   // All-in-all-out parameter
  double mT1_par[NUM_COMPS];    // Vettore di magnetizzazione  magn^alpha = sum_i S_i^alpha /N_siti
  double mT2[NUM_COMPS];
  double mE[NUM_COMPS-1];
  double mT1_perp[NUM_COMPS];
  double mSub;
  double nematic_op[NUM_COMPS-1];
  double bond_nematic[NUM_COMPS-1];
} obs_t;

typedef struct {                            // Struttura che contiene i nomi dei file per singola temperatura
  char data_file[STD_STRING_LENGTH];        // Nome del file su cui scrivo le misure
  char backup_file[STD_STRING_LENGTH];      // Nome del file su cui scrivo le configurazioni delle repliche
  char aux_file[STD_STRING_LENGTH];         // Nome dei file su cui scrivo tempToreplica
  char conf_file[STD_STRING_LENGTH];        // Nome dei file su cui appendo le config di spin per il fattore di struttura
} filenames_t;

// Funzioni e procedure, divise per utilizzo
// Inizializzazione (initialize.c)
void geometry(param_t *myParams);
void initialize_auxiliary_array(param_t *myParams, int *tempToreplica, int *exchAcc);
void initialize_local_axes();
void initialize_DM_vectors();
void initialize_v_vectors();
dsfmt_t * initialize_system(chain_t **myChains, param_t *myParams, obs_t **myObs, filenames_t **myFilenames, int **tempToreplica, int **exchAcc, dsfmt_t *dsfmt, char const *finput, chain_t *zero_T_replica);
void set_temperatures(chain_t *myChains, param_t *myParams, dsfmt_t *dsfmt, int *tempToreplica, int *exchAcc);
void check_metro_acc(chain_t *single_chain, param_t *myParams);
void set_acceptances(chain_t *myChains, param_t *myParams, dsfmt_t *dsfmt, int *tempToreplica);
void set_filenames(chain_t *myChains, param_t *myParams, filenames_t *myFilenames, int *tempToreplica);
void burn_in(chain_t *myChains, param_t *myParams, dsfmt_t *dsfmt, int *tempToreplica, int *exchAcc);
void deallocation(chain_t **myChains, param_t *myParams, obs_t **myObs, filenames_t **myFilenames, int **tempToreplica, int **exchAcc, dsfmt_t *dsfmt, chain_t *zero_T_replica);

// Gestione input e output (in_out.c)
void remove_white_line_and_comments(FILE *input);
void read_from_input_Param(param_t *Par, char const *finput);
void read_from_filenames(filenames_t *myFilenames, param_t *myParams);
void write_Status(param_t *myParams, chain_t *myChains, filenames_t *myFilenames, int *tempToreplica);
void read_Status(param_t *myParams, chain_t *myChains, filenames_t *myFilenames, int *tempToreplica);
void write_Log(param_t *myParams, chain_t *myChains, int *tempToreplica, int *exchAcc);
void save_Configs(chain_t *myChains, param_t *myParams, filenames_t *myFilenames, int *tempToreplica);
void save_zeroT_config(param_t *myParams, chain_t *single_chain);

// Gestione degli spin (utility.c)
void renormalize_spin(chain_t *single_chain, param_t *myParams);
void spin_difference(double *spin_1, double *spin_2, double *difference);
double spin_product(double *spin_1, double *spin_2);
void two_spin_sum(double *spin_1, double *spin_2, double *sum);
void three_spin_sum(double *spin_1, double *spin_2, double *spin_3, double *sum);
void scalar_X_vector(double scalar, double *spin_in, double *spin_out);
void spin_copy(double *spin_in, double *spin_out);
void resetErr();
void replica_copy(param_t *myParams, chain_t *single_chain, chain_t *copy);
void get_coordinates(param_t *myParams, int iSite, int *m1, int *m2, int *m3, int *mu);
int get_index(param_t *myParams, int m1, int m2, int m3, int mu);
void cross_product(double *spin_1, double *spin_2, double *result);

// Update (update.c)
void spin_trial(chain_t *single_chain, dsfmt_t *dsfmt, int rank, int iSite, double *trial);
void mean_field(chain_t *single_chain, param_t *myParams, int iSite, double *mean_field);
void update_metropolis(param_t *myParams, chain_t *single_chain, dsfmt_t *dsfmt, int rank);
void swap_temperatures(param_t *myParams, chain_t *myChains, dsfmt_t *dsfmt, int rank, int *tempToreplica, int *exchAcc);
void update_overrelax(param_t *myParams, chain_t *single_chain, dsfmt_t *dsfmt, int rank);
void iterative_minimization(param_t *myParams, chain_t *single_chain);
void push_to_zero(param_t *myParams, chain_t *low_T_chain, chain_t *zero_T_replica);

// Misure (measure.c)
double get_energy(chain_t *single_chain, param_t *myParams);
double get_energy_2(chain_t *single_chain, param_t *myParams);
void get_measures(chain_t *myChains, param_t *myParams, filenames_t *myFilenames, obs_t *myObs, int *tempToreplica, int *exchAcc, dsfmt_t *dsfmt);
void get_magns(chain_t *single_chain, param_t *myParams, obs_t *single_obs);
