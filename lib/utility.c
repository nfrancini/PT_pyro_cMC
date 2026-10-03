// Includo l'header
#include "../head/header.h"

// Procedura che rinormalizza gli spin di una singola replica ad 1
// Conviene una singola replica per parallelizzare meglio
void renormalize_spin(chain_t *single_chain, param_t *myParams){
  double modulus=0;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  for(int i=0;i<myParams->N;i++){
    modulus=0;
    // Calcolo la norma
    // dello spin nel sito i
    for(int j=0;j<NUM_COMPS;j++){
      modulus+=pow(single_chain->lattice[i][j],2);
      #ifdef DEBUG
      if(errno == EDOM){
        printf("Funzione renormalize, pow\n");
        perror("    errno == EDOM");
      }
      if(errno == ERANGE){
        printf("Funzione renormalize, pow\n");
        perror("    errno == ERANGE");
        if(fetestexcept(FE_OVERFLOW)){
          puts("    FE_OVERFLOW was raised");
          exit(EXIT_FAILURE);
        }
        else if(fetestexcept(FE_UNDERFLOW)){
          puts("    FE_UNDERFLOW was raised");
          exit(EXIT_FAILURE);
        }
      }
      if(fetestexcept(FE_INVALID)){
        puts("    FE_INVALID was raised");
        exit(EXIT_FAILURE);
      }
      if(fetestexcept(FE_DIVBYZERO)){
        puts("    FE_DIVBYZERO was raised");
        exit(EXIT_FAILURE);
      }
      #endif
    }
    modulus=sqrt(modulus);
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
    // Divido lo spin nel sito i
    // per la sua norma
    for(int j=0;j<NUM_COMPS;j++){
      single_chain->lattice[i][j]=single_chain->lattice[i][j]/modulus;
      #ifdef DEBUG
      if(errno == EDOM){
        printf("Funzione renormalize, divisione\n");
        perror("    errno == EDOM");
      }
      if(errno == ERANGE){
        printf("Funzione renormalize, divisione\n");
        perror("    errno == ERANGE");
        if(fetestexcept(FE_OVERFLOW)){
          puts("    FE_OVERFLOW was raised");
          exit(EXIT_FAILURE);
        }
        else if(fetestexcept(FE_UNDERFLOW)){
          puts("    FE_UNDERFLOW was raised");
          exit(EXIT_FAILURE);
        }
      }
      if(fetestexcept(FE_INVALID)){
        puts("    FE_INVALID was raised");
        exit(EXIT_FAILURE);
      }
      if(fetestexcept(FE_DIVBYZERO)){
        puts("    FE_DIVBYZERO was raised");
        exit(EXIT_FAILURE);
      }
      #endif
    }
  }
}

// Procedura che fa la differenza tra due vettori di spin \vec{d}=\vec{s1} - \vec{s2}
void spin_difference(double *spin_1, double *spin_2, double *difference){
  for(int i=0;i<NUM_COMPS;i++){
    difference[i] = spin_1[i]-spin_2[i];
  }
}

// Procedura che fa il prodotto scalare di due vettori di spin p= \vec{s1}_i * \vec{s2}_i
double spin_product(double *spin_1, double *spin_2){
  double sum=0;
  for(int i=0;i<NUM_COMPS;i++){
    sum+=spin_1[i]*spin_2[i];
  }
  return sum;
}

// Procedure for spin cross product
void cross_product(double *spin_1, double *spin_2, double *result){
  result[0]=spin_1[1]*spin_2[2] - spin_1[2]*spin_2[1];
  result[1]=spin_1[2]*spin_2[0] - spin_1[0]*spin_2[2];
  result[2]=spin_1[0]*spin_2[1] - spin_1[1]*spin_2[0];
}

// Procedura che somma due spin
void two_spin_sum(double *spin_1, double *spin_2, double *sum){
  for(int i=0;i<NUM_COMPS;i++){
    sum[i] = spin_1[i] + spin_2[i];
  }
}

// Procedura che somma tre spin
void three_spin_sum(double *spin_1, double *spin_2, double *spin_3, double *sum){
  for(int i=0;i<NUM_COMPS;i++){
    sum[i] = spin_1[i] + spin_2[i] + spin_3[i];
  }
}

// Procedura che moltiplica un vettore per uno scalare
void scalar_X_vector(double scalar, double *spin_in, double *spin_out){
  for(int i=0;i<NUM_COMPS;i++){
    spin_out[i]=scalar*spin_in[i];
  }
}

// Funzione che copia uno spin in un altro
void spin_copy(double *spin_in, double *spin_out){
  for(int i=0;i<NUM_COMPS;i++){
    spin_out[i]=spin_in[i];
  }
}

// Procedure that gets the coordinates (m1,m2,m3;mu) from
// the super index iSite.
// The position of a site is then given by the expression
// R=m1*t1 + m2*t2 + m3*t3 + d_mu
// where the first 3 coefficients give the FCC structure (aka the sublattice 0)
// and mu specifies the sublattice
// t1, t2, t3 and d_mu are known
// t1=(110)/2, t2=(101)/2, t3=(011)/2 and
// d_0=0, d_1=t1/2, d_2=t2/2. d3=t3/2
void get_coordinates(param_t *myParams, int iSite, int *m1, int *m2, int *m3, int *mu){
  *m1=iSite%myParams->size;
  *m2=(iSite/myParams->size)%myParams->size;
  *m3=(iSite/((int)pow(myParams->size,2)))%myParams->size;
  *mu=(iSite/((int)pow(myParams->size,3)));
}

// Procedure to get the super index iSite starting
// from the tuple (m1,m2,m3,mu)
int get_index(param_t *myParams, int m1, int m2, int m3, int mu){
  int iSite;
  iSite=m1 + m2*myParams->size + m3*((int)pow(myParams->size, 2)) + mu*((int)pow(myParams->size,3));
  return iSite;
}

// Procedura che resetta gli errori
void resetErr(){
  errno = 0;
  feclearexcept(FE_ALL_EXCEPT);
}

// Procedura che copia una replica in un'altra
void replica_copy(param_t *myParams, chain_t *single_chain, chain_t *copy){
  copy->acc=0;
  copy->acc_Parameter=0;
  copy->count=0;
  copy->T=0;
  for(int i=0;i<myParams->N;i++){
    for(int j=0; j<NUM_COMPS;j++){
      copy->lattice[i][j]=single_chain->lattice[i][j];
    }
  }
}
