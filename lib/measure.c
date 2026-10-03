// Includo l'header
#include "../head/header.h"

// Funzione che calcola l'energia di una singola replica
double get_energy_2(chain_t *single_chain, param_t *myParams){
  double aux_spin[NUM_COMPS];
  double heis_term = 0.0;
  double ising_term = 0.0;
  double PD_term = 0.0;
  double DM_term = 0.0;
  double ene=0.0;
  // sum on every site
  for(int i=0;i<myParams->N;i++){
    int m1,m2,m3,mu;
    get_coordinates(myParams, i, &m1, &m2, &m3, &mu);
    for(int j=0;j<NUM_NEAREST_NEIGHBOR;j++){
      int n1,n2,n3,nu;
      get_coordinates(myParams, nn[i][j], &n1, &n2, &n3, &nu);
      cross_product(single_chain->lattice[i], single_chain->lattice[nn[i][j]], aux_spin);
      heis_term+=spin_product(single_chain->lattice[i], single_chain->lattice[nn[i][j]]);
      ising_term +=spin_product(single_chain->lattice[i], local_z[mu])*spin_product(single_chain->lattice[nn[i][j]], local_z[nu]);
      PD_term +=spin_product(single_chain->lattice[i], rr[i][j])*(spin_product(single_chain->lattice[nn[i][j]], rr[i][j]));
      DM_term +=spin_product(DM_vectors[mu][nu], aux_spin);
    }
  }
  ene = (myParams->heis)*heis_term/2.0 + (myParams->ising)*ising_term/2.0 + (myParams->PD)*PD_term/2.0 + (myParams->DM)*DM_term/2.0;
  return ene;
}

double get_energy(chain_t *single_chain, param_t *myParams){
  double ene=0.0;
  double mean_f[NUM_COMPS];
  // Energy with mean_field
  for(int iSite=0;iSite<myParams->N;iSite++){
    mean_field(single_chain, myParams, iSite, mean_f);
    ene+=spin_product(single_chain->lattice[iSite], mean_f);
  }
  ene=ene/2.0;
  return ene;
}

// Procedura che misura le varie magnetizzazioni del reticolo
// Utilizzo una supercella da 8 siti, quindi devo tenerne conto nella
// misura delle varie magnetizzazioni
void get_magns(chain_t *single_chain, param_t *myParams, obs_t *single_obs){
  double aux_spin[NUM_COMPS];
  double ms[NUM_SUBLATTICES][NUM_COMPS];
  double sum_Ex=0;
  double sum_Ey=0;
  double sum_nemx=0;
  double sum_nemy=0;
  double bond_nemx=0;
  double bond_nemy=0;
  //Azzero gli errori se DEBUG è definito
  #ifdef DEBUG
  resetErr();
  #endif
  // Azzero tutte le misure
  single_obs->mA2=0;
  single_obs->mSub=0;
  for(int i=0;i<NUM_COMPS;i++){
    single_obs->mT1_par[i]=0;
    single_obs->mT2[i]=0;
    single_obs->mT1_perp[i]=0;
    for(int k=0; k<NUM_SUBLATTICES; k++){
      ms[k][i]=0.0;
    }
    if(i<NUM_COMPS-1){
      single_obs->mE[i]=0;
      single_obs->nematic_op[i]=0.0;
      single_obs->bond_nematic[i]=0.0;
    }
  }

  for(int i=0;i<myParams->N;i++){
    int m1,m2,m3,mu;
    get_coordinates(myParams,i,&m1,&m2,&m3,&mu);
    cross_product(local_z[mu], single_chain->lattice[i], aux_spin);
    single_obs->mA2 += spin_product(single_chain->lattice[i], local_z[mu])/myParams->N;
    sum_Ex += spin_product(local_x[mu], single_chain->lattice[i])/myParams->N;
    sum_Ey += spin_product(local_y[mu], single_chain->lattice[i])/myParams->N;
    sum_nemx += (pow(spin_product(local_x[mu], single_chain->lattice[i]),2) - pow(spin_product(local_y[mu], single_chain->lattice[i]),2))/myParams->N;
    sum_nemy += 2.0*(spin_product(local_x[mu], single_chain->lattice[i]))*(spin_product(local_y[mu], single_chain->lattice[i]))/myParams->N;
    for(int k=0;k<NUM_NEAREST_NEIGHBOR;k++){
      bond_nemx += ((single_chain->lattice[i][0]*single_chain->lattice[nn[i][k]][0]) - (single_chain->lattice[i][1]*single_chain->lattice[nn[i][k]][1]))/(6*myParams->N);
      bond_nemy += ((single_chain->lattice[i][0]*single_chain->lattice[nn[i][k]][1]) + (single_chain->lattice[i][1]*single_chain->lattice[nn[i][k]][0]))/(6*myParams->N);
    }
    for(int j=0;j<NUM_COMPS;j++){
      single_obs->mT1_par[j] += single_chain->lattice[i][j]/myParams->N;
      single_obs->mT2[j] += aux_spin[j]*sqrt(3.0/2.0)/myParams->N;
      ms[mu][j] += single_chain->lattice[i][j]*4.0/myParams->N;
      for(int k=0;k<NUM_COMPS;k++){
        for(int l=0;l<NUM_COMPS;l++){
          if((k<l)&&(j!=k)&&(j!=l)){
            single_obs->mT1_perp[j] += local_z[mu][j]*spin_product(vv[mu][k][l], single_chain->lattice[i])*sqrt(3.0)/myParams->N;
          }
        }
      }
    }
  }
  single_obs->mE[0]=sum_Ex;
  single_obs->mE[1]=sum_Ey;
  single_obs->nematic_op[0]=sum_nemx;
  single_obs->nematic_op[1]=sum_nemy;
  single_obs->bond_nematic[0]=bond_nemx;
  single_obs->bond_nematic[1]=bond_nemy;
  for(int mu=0;mu<NUM_SUBLATTICES;mu++){
    single_obs->mSub += spin_product(ms[mu], ms[mu])/4.0;
  }
}

// Procedura per le misure
void get_measures(chain_t *myChains, param_t *myParams, filenames_t *myFilenames, obs_t *myObs, int *tempToreplica, int *exchAcc, dsfmt_t *dsfmt){
  FILE *fOut[myParams->n_temp];
  int rank;
  int flag;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  // Reset delle accettanze degli swap
  for(int i=0;i<myParams->n_temp-1;i++){
    exchAcc[i]=0;
  }
  // Apro i file delle misure
  for(int i=0;i<myParams->n_temp;i++){
    fOut[i]=fopen((myFilenames + i)->data_file, "a");
    if(fOut[i]==NULL){
      perror("Errore in apertura in get_measures");
      exit(EXIT_FAILURE);
    }
  }
  // Ciclo di misurazioni
  for(int j=0;j<myParams->iMis;j++){
    for(int k=0;k<myParams->iDec;k++){
      // Parallelizzo per il Metropolis
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          rank=omp_get_thread_num();
          update_metropolis(myParams, myChains+tempToreplica[i], dsfmt, rank);
        }
        #pragma omp barrier
      }
      // Overrelaxation
      for(int j=0;j<IOVERR;j++){
        #pragma omp parallel
        {
          #pragma omp for
          for(int i=0;i<myParams->n_temp;i++){
            rank=omp_get_thread_num();
            update_overrelax(myParams, myChains+tempToreplica[i], dsfmt, rank);
          }
          #pragma omp barrier
        }
      }
      swap_temperatures(myParams, myChains, dsfmt, rank, tempToreplica, exchAcc);
    }
    // Misure in parallelo
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0;i<myParams->n_temp;i++){
        rank=omp_get_thread_num();
        (myObs+i)->ene_density = get_energy((myChains + tempToreplica[i]), myParams)/myParams->N;
        // Misuro le magnetizzazioni che mi interessano
        get_magns((myChains + tempToreplica[i]), myParams, (myObs+i));
        #ifdef DEBUG
        if(errno == EDOM){
          printf("Funzione get_magn, divisione\n");
          perror("    errno == EDOM");
        }
        if(errno == ERANGE){
          printf("Funzione get_magn, divisione\n");
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
      #pragma omp barrier
    }
    // Scrivo le misure
    // Apro i file delle misure
    for(int i=0;i<myParams->n_temp;i++){
      fprintf(fOut[i], "%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", (myObs+i)->ene_density, (myObs+i)->mSub, (myObs+i)->mA2, (myObs+i)->mT1_par[0], (myObs+i)->mT1_par[1], (myObs+i)->mT1_par[2], (myObs+i)->mT2[0], (myObs+i)->mT2[1], (myObs+i)->mT2[2], (myObs+i)->mE[0], (myObs+i)->mE[1], (myObs+i)->mT1_perp[0], (myObs+i)->mT1_perp[1], (myObs+i)->mT1_perp[2], (myObs+i)->nematic_op[0], (myObs+i)->nematic_op[1],(myObs+i)->bond_nematic[0],(myObs+i)->bond_nematic[1]);
    }
    // Rinormalizzo gli spin ogni 10 000 misure
    if((j+1)%INORM==0){
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          renormalize_spin((myChains +tempToreplica[i]), myParams);
        }
        #pragma omp barrier
      }
    }
    // Salvo le configurazioni a T finita
    if((myParams->iSave!=0)&&((j+1)%myParams->iSave==0)){
      save_Configs(myChains, myParams, myFilenames, tempToreplica);
    }
    // Backup
    if((j+1)%myParams->iBackup==0){
      // Rimuovo il file contentente i nomi e lo riscrivo
      flag=remove(myParams->filenames_file);
      if(flag!=0){
        perror("Errore nel rimuovere filenames in get_measure\n");
        exit(EXIT_FAILURE);
      }
      write_Status(myParams, myChains, myFilenames, tempToreplica);
    }
  }
  // Chiudo i file delle misure
  for(int i=0;i<myParams->n_temp;i++){
    fclose(fOut[i]);
  }
}
