// Includo l'header
#include "../head/header.h"

// Procedura per creare lo spin di prova
// La procedura prende in ingresso una singola catena
// il sito in cui viene creato lo spin di prova
// ed ovviamente i parametri di sistema ed un vettore trial
// dove salvo lo spin di prova.
// L'idea è di usare gli angoli di Eulero:
// Lo spin scelto può essere visto come combinazione di due
// rotazioni del vettore (001). Quindi prima genero un
// cono di punti uniformi attorno ad (001), poi lo ruoto sfruttando
// le rotazioni del vettore scelto. Così ottengo un cono
// di punti casuali attorno al punto scelto.
void spin_trial(chain_t *single_chain, dsfmt_t *dsfmt, int rank, int iSite, double *trial){
  // Metodo Angoli di Eulero
  double z_delta, phi_delta;
  double x_delta, y_delta;
  double c1, s1;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  // Mi servono due angoli per generare un punto random
  // nel cono sferico attorno ad (001). Qui il parametro acc_Parameter
  // regola l'ampiezza del cono da un punto (acc_Parameter=0) a tutta la sfera (acc_Parameter=1.0)
  // x_delta, y_delta, z_delta sono le coordinate del punto random nel cono attorno a (001)
  // 1-2*acc_Parameter < z_delta < 1, phi_delta \in [0,2PI)
  z_delta = 1.0 - 2*single_chain->acc_Parameter*dsfmt_genrand_close_open(&dsfmt[rank]);
  phi_delta = 2*PI*dsfmt_genrand_close_open(&dsfmt[rank]);
  x_delta=sin(acos(z_delta))*cos(phi_delta);
  y_delta=sin(acos(z_delta))*sin(phi_delta);
  #ifdef DEBUG
  if(errno == EDOM){
    printf("Funzione renormalize, sin o cos\n");
    perror("    errno == EDOM");
  }
  if(errno == ERANGE){
    printf("Funzione renormalize, sin o cos\n");
    perror("    errno == ERANGE");
  }
  if(fetestexcept(FE_INVALID)){
    puts("    FE_INVALID was raised");
    exit(EXIT_FAILURE);
  }
  #endif
  // Coordinate "sferiche" del punto iniziale
  c1=cos(atan2(single_chain->lattice[iSite][1], single_chain->lattice[iSite][0]));
  s1=sin(atan2(single_chain->lattice[iSite][1], single_chain->lattice[iSite][0]));
  #ifdef DEBUG
  if(errno == EDOM){
    printf("Funzione renormalize, sin o cos\n");
    perror("    errno == EDOM");
  }
  if(errno == ERANGE){
    printf("Funzione renormalize, sin o cos\n");
    perror("    errno == ERANGE");
  }
  if(fetestexcept(FE_INVALID)){
    puts("    FE_INVALID was raised");
    exit(EXIT_FAILURE);
  }
  #endif
  // Creo l'array di prova applicando le rotazioni con gli angoli di Eulero ottenuti del vettore
  // in iSite
  trial[0]=single_chain->lattice[iSite][2]*c1*x_delta - s1*y_delta + single_chain->lattice[iSite][0]*z_delta;
  trial[1]=single_chain->lattice[iSite][2]*s1*x_delta + c1*y_delta + single_chain->lattice[iSite][1]*z_delta;
  trial[2]=-sqrt(1-pow(single_chain->lattice[iSite][2], 2))*x_delta + single_chain->lattice[iSite][2]*z_delta;
  #ifdef DEBUG
  if(errno == EDOM){
    printf("Funzione renormalize, sqrt\n");
    perror("    errno == EDOM");
  }
  if(errno == ERANGE){
    printf("Funzione renormalize, sqrt\n");
    perror("    errno == ERANGE");
  }
  if(fetestexcept(FE_INVALID)){
    puts("    FE_INVALID was raised");
    exit(EXIT_FAILURE);
  }
  #endif
}

// Procedura per creare il campo medio
// Prima calcolo la parte relativa al modello HKGG'
void mean_field(chain_t *single_chain, param_t *myParams, int iSite, double *mean_field){
  int m1,m2,m3,mu;
  int n1,n2,n3,nu;
  double aux_spin[NUM_COMPS];
  // Inizializzo a zero
  for(int i=0;i<NUM_COMPS;i++){
    mean_field[i]=0.0;
  }
  get_coordinates(myParams, iSite, &m1, &m2, &m3, &mu);     // Coordinate iSite
  for(int j=0;j<NUM_NEAREST_NEIGHBOR;j++){
    get_coordinates(myParams, nn[iSite][j], &n1,&n2,&n3,&nu);
    cross_product(single_chain->lattice[nn[iSite][j]], DM_vectors[mu][nu], aux_spin);
    for(int i=0;i<NUM_COMPS;i++){
      mean_field[i]+=(myParams->heis)*(single_chain->lattice[nn[iSite][j]][i]) + (myParams->ising)*(local_z[mu][i])*spin_product(single_chain->lattice[nn[iSite][j]],local_z[nu]) + (myParams->PD)*(rr[iSite][j][i])*(spin_product(single_chain->lattice[nn[iSite][j]],rr[iSite][j])) + (myParams->DM)*aux_spin[i];
    }
  }
}

// Update Metropolis
void update_metropolis(param_t *myParams, chain_t *single_chain, dsfmt_t *dsfmt, int rank){
  double trial[NUM_COMPS], diff[NUM_COMPS], mean_f[NUM_COMPS];
  double c=0;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  // Ciclo sui tutti i siti
  for(int iSite=0;iSite<myParams->N;iSite++){
    // Aumento il contatore
    single_chain->count+=1;
    // Creo lo spin di prova
    spin_trial(single_chain, dsfmt, rank, iSite, trial);
    // Differenza tra spin di prova e spin iniziale nel sito iSite
    spin_difference(trial, single_chain->lattice[iSite], diff);
    // Calcolo il campo medio nel sito iSite
    mean_field(single_chain, myParams, iSite, mean_f);
    // Differenza di energia
    c=-spin_product(diff, mean_f)/single_chain->T;
    #ifdef DEBUG
    double ene1;
    ene1=get_energy(single_chain, myParams);
    double dE=spin_product(diff, mean_f);
    #endif
    // Test dell'accettanza
    if((c>=0) || (exp(c)>dsfmt_genrand_close_open(&dsfmt[rank]))){
      #ifdef DEBUG
      if(errno == ERANGE){
        printf("Funzione update_metro\n");
        perror("    errno == ERANGE");
      }
      if(fetestexcept(FE_OVERFLOW)){
        printf("    FE_OVERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_UNDERFLOW)){
        printf("    FE_UNDERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_DIVBYZERO)){
        printf("    FE_DIVBYZERO was raised\n");
        exit(EXIT_FAILURE);
      }
      #endif
      // Accetto il passo e copio lo spin di prova
      spin_copy(trial, single_chain->lattice[iSite]);
      single_chain->acc+=1;
      #ifdef DEBUG
      double ene2;
      ene2=get_energy(single_chain, myParams);
      double delta_glob=ene2-ene1;
      if(fabs(delta_glob-dE)>1.0e-12){
        printf("Errore, ene_locale %lf\tene_globale %lf\n", dE, delta_glob);
        exit(EXIT_FAILURE);
      }
      #endif
    }
  }
}

// Swap delle temperature
void swap_temperatures(param_t *myParams, chain_t *myChains, dsfmt_t *dsfmt, int rank, int *tempToreplica, int *exchAcc){
  double delta;
  double tempT;
  double tempAccparam;
  int tempindex;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  // Ciclo sulle temperature, m è un indice di temperatura
  for(int m=0; m<myParams->n_temp-1;m++){
    // Condizione di swapping tra T_m e T_{m+1}
    delta = (1.0/myChains[tempToreplica[m+1]].T - 1.0/myChains[tempToreplica[m]].T) * (get_energy((myChains + tempToreplica[m]), myParams) - get_energy((myChains + tempToreplica[m+1]), myParams));
    if((delta<0) || (exp(-delta)>dsfmt_genrand_close_open(&dsfmt[rank]))){
      #ifdef DEBUG
      if(errno == ERANGE){
        printf("Funzione swap_temperatures\n");
        perror("    errno == ERANGE");
      }
      if(fetestexcept(FE_OVERFLOW)){
        printf("    FE_OVERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_UNDERFLOW)){
        printf("    FE_UNDERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_DIVBYZERO)){
        printf("    FE_DIVBYZERO was raised\n");
        exit(EXIT_FAILURE);
      }
      #endif
      // Swap avvenuto
      exchAcc[m]+=1;
      // Swap delle temperature
      tempT=myChains[tempToreplica[m]].T;
      myChains[tempToreplica[m]].T=myChains[tempToreplica[m+1]].T;
      myChains[tempToreplica[m+1]].T=tempT;
      // Swap del parametro di accettanza
      tempAccparam=myChains[tempToreplica[m]].acc_Parameter;
      myChains[tempToreplica[m]].acc_Parameter=myChains[tempToreplica[m+1]].acc_Parameter;
      myChains[tempToreplica[m+1]].acc_Parameter = tempAccparam;
      // Tracking dello swap per tempToreplica
      tempindex=tempToreplica[m];
      tempToreplica[m]=tempToreplica[m+1];
      tempToreplica[m+1]=tempindex;
    }
    // Swap rifiutato
  }
}

// Overrelaxation
void update_overrelax(param_t *myParams, chain_t *single_chain, dsfmt_t *dsfmt, int rank){
  int iSite, j;
  double mean_f[NUM_COMPS], trial[NUM_COMPS];
  double norm;

  for(iSite=0; iSite<myParams->N; iSite++){
    #ifdef DEBUG
    double ene1;
    ene1=get_energy(single_chain, myParams);
    #endif
    // Calculate mean field
    mean_field(single_chain, myParams, iSite, mean_f);
    norm=spin_product(mean_f, mean_f);
    // If the norm is too small
    if(norm<1.0e-12){
      for(j=0;j<NUM_COMPS;j++){
        trial[j] = 2*dsfmt_genrand_close_open(&dsfmt[rank])-1.0;
      }
      for(j=0;j<NUM_COMPS;j++){
        trial[j] = trial[j]/sqrt(spin_product(trial, trial));
      }
      spin_copy(trial, single_chain->lattice[iSite]);
      return;
    }
    else{
      for(j=0;j<NUM_COMPS;j++){
        trial[j]=2.0*spin_product(single_chain->lattice[iSite], mean_f)*mean_f[j]/norm - single_chain->lattice[iSite][j];
      }
      spin_copy(trial, single_chain->lattice[iSite]);
      //
      #ifdef DEBUG
      double ene2;
      ene2=get_energy(single_chain, myParams);
      if(fabs(ene1-ene2)>1e-12){
        if((fabs(ene1-ene2)>1.0e-12)&(fabs(ene1-ene2)<1.0e-11)){
          err1 = err1+1;
        }
        else if((fabs(ene1-ene2)>1.0e-11)){
          err2=err2+1;
        }
      }
      #endif
    }
  }
}

// Minimizzazione iterativa
void iterative_minimization(param_t *myParams, chain_t *single_chain){
  double mean_f[NUM_COMPS];
  double trial[NUM_COMPS];
  double diff[NUM_COMPS];
  double norm=0;
  double temp;
  double max_diff=0;
  int exit=0;
  unsigned long iter=0;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  while(exit==0){
    max_diff=0;
    iter+=1;
    // Spazzata sul reticolo
    for(int iSite=0;iSite<myParams->N;iSite++){
      norm=0;
      // Campo medio in iSite
      mean_field(single_chain, myParams, iSite, mean_f);
      // Calcolo la norma del campo medio e lo rinormalizzo
      for(int j=0;j<NUM_COMPS;j++){
        norm+=pow(mean_f[j],2);
      }
      norm=sqrt(norm);
      // Creo lo spin di prova
      for(int j=0;j<NUM_COMPS;j++){
        trial[j]=-mean_f[j]/norm;
      }
      // Calcolo la differenza tra spin di prova e spin originale
      spin_difference(single_chain->lattice[iSite], trial, diff);
      temp=sqrt(spin_product(diff,diff));
      if(temp>max_diff){
        max_diff=temp;
      }
      // Copio lo spin di prova nel sito iSite
      spin_copy(trial, single_chain->lattice[iSite]);
    }
    if(fabs(max_diff)<=1.0e-12){
      exit=1;
    }
  }
}

// Procedura di minimizzazione
void push_to_zero(param_t *myParams, chain_t *low_T_chain, chain_t *zero_T_replica){
  double min_energy=0;
  double temp;
  FILE *fptr;
  char *file_name[STD_STRING_LENGTH];
  // Creo il nome del file per le misurazioni
  for(int i=0;i<myParams->iTerm;i++){
    // Prima copio la replica a temperatura zero
    replica_copy(myParams, low_T_chain, zero_T_replica);
    // Minimizzazione sulla copia
    iterative_minimization(myParams, zero_T_replica);
    // Misuro l'energia
    temp=get_energy(zero_T_replica, myParams);
    if(temp<min_energy){
      min_energy=temp;
      // Scrivo la configurazione di bassa temperatura
      // Comprese le misure
      save_zeroT_config(myParams, zero_T_replica);
    }
  }
}
