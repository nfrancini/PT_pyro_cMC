// Includo l'header
#include "../head/header.h"

// Geometry procedure that creates the nearest neighbor structure
// and fixes the periodi boundary conditions.
void geometry(param_t *myParams){
  int m1, m2, m3, mu;
  // Lexicographic index
  // iSite = m1 + m2*L + m3*L^2 + mu*L^3
  // with mi in [0,L-1] and mu=0,1,2,3.
  // The geometry is encoded in nn[iSite][dir], where
  // the array nn takes as first index iSite and as second index the
  // nn direction. The value saved is then the Lexicographic index for
  // the corresponding nearest neighbor
  // I need a convenction for the directions of links from iSite to nearest neighbor.
  // The convention is that nn[iSite][0] correspond to the smallest nearest neighbor sublattice
  // and nn[iSite][5] to the largest nearest neighbor sublattice.
  // The the values are saved on "same" tetrahedra and "different" tetrahedra
  // which means that the nearest neighbor tuple with same (m1,m2,m3) values are saved
  // for first.
  // For example, taking an iSite with mu=1, I have
  //  nn[iSite][0]--->(m1,m2,m3;0)      nn with mu=0
  //  nn[iSite][1]--->(m1+1,m2,m3;0)
  //  nn[iSite][2]--->(m1,m2,m3;2)      nn with mu=2
  //  nn[iSite][3]--->(m1+1,m2-1,m3;2)
  //  nn[iSite][4]--->(m1,m2,m3;3)      nn with mu=3
  //  nn[iSite][5]--->(m1+1,m2,m3-1;3)
  // Of course the nearest neighbor tuple is first treated with PBC and then
  // converted to the super index.
  for(int iSite=0;iSite<myParams->N;iSite++){
    get_coordinates(myParams,iSite,&m1,&m2,&m3,&mu);
    // 0-th sublattice
    if(mu==0){
      // mu=1 nearest neighbors
      nn[iSite][0] = get_index(myParams,m1,m2,m3,1); // same tetrahedra (m1,m2,m3;1)
      if(m1==0){
        // diff tetrahedra (m1-1,m2,m3;1) with PBC m1=0 ---> m1-1=L-1
        nn[iSite][1] = get_index(myParams,myParams->size-1,m2,m3,1);
      }
      else{
        nn[iSite][1] = get_index(myParams,m1-1,m2,m3,1);
      }
      // mu=2 nearest neighbors
      nn[iSite][2] = get_index(myParams,m1,m2,m3,2);  // same tetrahedra (m1,m2,m3;2)
      if(m2==0){
        // diff tetrahedra (m1,m2-1,m3;2) with PBC m2=0 ---> m2-1=L-1
        nn[iSite][3] = get_index(myParams,m1,myParams->size-1,m3,2);
      }
      else{
        nn[iSite][3] = get_index(myParams,m1,m2-1,m3,2);
      }
      // mu=3 nearest neighbors
      nn[iSite][4] = get_index(myParams,m1,m2,m3,3);  // same tetrahedra (m1,m2,m3;3)
      if(m3==0){
        // diff tetrahedra (m1,m2,m3-1;3) with PBC m3=0 ---> m3-1=L-1
        nn[iSite][5] = get_index(myParams,m1,m2,myParams->size-1,3);
      }
      else{
        nn[iSite][5] = get_index(myParams,m1,m2,m3-1,3);
      }
      // r_ij, i in mu=0
      rr[iSite][0][0]=1.0/sqrt(2);    // same FCC mu=1 NN
      rr[iSite][0][1]=1.0/sqrt(2);
      rr[iSite][0][2]=0;

      rr[iSite][1][0]=-1.0/sqrt(2);    // different FCC mu=1 NN
      rr[iSite][1][1]=-1.0/sqrt(2);
      rr[iSite][1][2]=0;

      rr[iSite][2][0]=1.0/sqrt(2);    // same FCC mu=2 NN
      rr[iSite][2][1]=0;
      rr[iSite][2][2]=1.0/sqrt(2);

      rr[iSite][3][0]=-1.0/sqrt(2);    // different FCC mu=2 NN
      rr[iSite][3][1]=0;
      rr[iSite][3][2]=-1.0/sqrt(2);

      rr[iSite][4][0]=0;               // same FCC mu=3 NN
      rr[iSite][4][1]=1.0/sqrt(2);
      rr[iSite][4][2]=1.0/sqrt(2);

      rr[iSite][5][0]=0;               // different FCC mu=3 NN
      rr[iSite][5][1]=-1.0/sqrt(2);
      rr[iSite][5][2]=-1.0/sqrt(2);
    }
    // 1-st sublattice
    else if(mu==1){
      // mu=0
      nn[iSite][0] = get_index(myParams,m1,m2,m3,0); // same tetrahedra (m1,m2,m3;0)
      if(m1==myParams->size-1){
        // diff tetrahedra (m1+1,m2,m3;0) with PBC m1=L-1 ---> m1+1=0
        nn[iSite][1] = get_index(myParams,0,m2,m3,0);
      }
      else{
        nn[iSite][1] = get_index(myParams,m1+1,m2,m3,0);
      }
      // mu=2
      nn[iSite][2] = get_index(myParams,m1,m2,m3,2);  // same tetrahedra (m1,m2,m3;2)
      if(m1==myParams->size-1){
        if(m2==0){
          nn[iSite][3] = get_index(myParams,0,myParams->size-1,m3,2);
        }
        else{
          nn[iSite][3] = get_index(myParams,0,m2-1,m3,2);
        }
      }
      else{
        if(m2==0){
          nn[iSite][3] = get_index(myParams,m1+1,myParams->size-1,m3,2);
        }
        else{
          nn[iSite][3] = get_index(myParams,m1+1,m2-1,m3,2);
        }
      }
      // mu = 3
      nn[iSite][4] = get_index(myParams,m1,m2,m3,3);
      if(m1==myParams->size-1){
        if(m3==0){
          nn[iSite][5] = get_index(myParams,0,m2,myParams->size-1,3);
        }
        else{
          nn[iSite][5] = get_index(myParams,0,m2,m3-1,3);
        }
      }
      else{
        if(m3==0){
          nn[iSite][5] = get_index(myParams,m1+1,m2,myParams->size-1,3);
        }
        else{
          nn[iSite][5] = get_index(myParams,m1+1,m2,m3-1,3);
        }
      }
      // r_ij, i in mu=1
      rr[iSite][0][0]=-1.0/sqrt(2);    // same FCC mu=0 NN
      rr[iSite][0][1]=-1.0/sqrt(2);
      rr[iSite][0][2]=0;

      rr[iSite][1][0]=1.0/sqrt(2);    // different FCC mu=0 NN
      rr[iSite][1][1]=1.0/sqrt(2);
      rr[iSite][1][2]=0;

      rr[iSite][2][0]=0;              // same FCC mu=2 NN
      rr[iSite][2][1]=-1.0/sqrt(2);
      rr[iSite][2][2]=1.0/sqrt(2);

      rr[iSite][3][0]=0;              // different FCC mu=2 NN
      rr[iSite][3][1]=1.0/sqrt(2);
      rr[iSite][3][2]=-1.0/sqrt(2);

      rr[iSite][4][0]=-1.0/sqrt(2);   // same FCC mu=3 NN
      rr[iSite][4][1]=0;
      rr[iSite][4][2]=1.0/sqrt(2);

      rr[iSite][5][0]=1.0/sqrt(2);    // different FCC mu=3 NN
      rr[iSite][5][1]=0;
      rr[iSite][5][2]=-1.0/sqrt(2);
    }
    // 2-nd sublattice
    else if(mu==2){
      // mu=0
      nn[iSite][0] = get_index(myParams,m1,m2,m3,0);  // same tetrahedra
      if(m2==myParams->size-1){
        nn[iSite][1] = get_index(myParams,m1,0,m3,0);
      }
      else{
        nn[iSite][1] = get_index(myParams,m1,m2+1,m3,0);
      }
      // mu=1
      nn[iSite][2] = get_index(myParams,m1,m2,m3,1);
      if(m1==0){
        if(m2==myParams->size-1){
          nn[iSite][3] = get_index(myParams,myParams->size-1,0,m3,1);
        }
        else{
          nn[iSite][3] = get_index(myParams,myParams->size-1,m2+1,m3,1);
        }
      }
      else{
        if(m2==myParams->size-1){
          nn[iSite][3] = get_index(myParams,m1-1,0,m3,1);
        }
        else{
          nn[iSite][3] = get_index(myParams,m1-1,m2+1,m3,1);
        }
      }
      // mu=3
      nn[iSite][4] = get_index(myParams,m1,m2,m3,3);
      if(m2==myParams->size-1){
        if(m3==0){
          nn[iSite][5] = get_index(myParams,m1,0,myParams->size-1,3);
        }
        else{
          nn[iSite][5] = get_index(myParams,m1,0,m3-1,3);
        }
      }
      else{
        if(m3==0){
          nn[iSite][5] = get_index(myParams,m1,m2+1,myParams->size-1,3);
        }
        else{
          nn[iSite][5] = get_index(myParams,m1,m2+1,m3-1,3);
        }
      }
      // r_ij, i in mu=2
      rr[iSite][0][0]=-1.0/sqrt(2);    // same FCC mu=0 NN
      rr[iSite][0][1]=0;
      rr[iSite][0][2]=-1.0/sqrt(2);

      rr[iSite][1][0]=1.0/sqrt(2);    // different FCC mu=0 NN
      rr[iSite][1][1]=0;
      rr[iSite][1][2]=1.0/sqrt(2);

      rr[iSite][2][0]=0;              // same FCC mu=1 NN
      rr[iSite][2][1]=1.0/sqrt(2);
      rr[iSite][2][2]=-1.0/sqrt(2);

      rr[iSite][3][0]=0;              // different FCC mu=1 NN
      rr[iSite][3][1]=-1.0/sqrt(2);
      rr[iSite][3][2]=1.0/sqrt(2);

      rr[iSite][4][0]=-1.0/sqrt(2);   // same FCC mu=3 NN
      rr[iSite][4][1]=1.0/sqrt(2);
      rr[iSite][4][2]=0;

      rr[iSite][5][0]=1.0/sqrt(2);    // different FCC mu=3 NN
      rr[iSite][5][1]=-1.0/sqrt(2);
      rr[iSite][5][2]=0;
    }
    else if(mu==3){
      // mu=0
      nn[iSite][0] = get_index(myParams,m1,m2,m3,0);
      if(m3==myParams->size-1){
        nn[iSite][1] = get_index(myParams,m1,m2,0,0);
      }
      else{
        nn[iSite][1] = get_index(myParams,m1,m2,m3+1,0);
      }
      // mu=1
      nn[iSite][2] = get_index(myParams,m1,m2,m3,1);
      if(m1==0){
        if(m3==myParams->size-1){
          nn[iSite][3] = get_index(myParams,myParams->size-1,m2,0,1);
        }
        else{
          nn[iSite][3] = get_index(myParams,myParams->size-1,m2,m3+1,1);
        }
      }
      else{
        if(m3==myParams->size-1){
          nn[iSite][3] = get_index(myParams,m1-1,m2,0,1);
        }
        else{
          nn[iSite][3] = get_index(myParams,m1-1,m2,m3+1,1);
        }
      }
      // mu=2
      nn[iSite][4] = get_index(myParams,m1,m2,m3,2);
      if(m2==0){
        if(m3==myParams->size-1){
          nn[iSite][5] = get_index(myParams,m1,myParams->size-1,0,2);
        }
        else{
          nn[iSite][5] = get_index(myParams,m1,myParams->size-1,m3+1,2);
        }
      }
      else{
        if(m3==myParams->size-1){
          nn[iSite][5] = get_index(myParams,m1,m2-1,0,2);
        }
        else{
          nn[iSite][5] = get_index(myParams,m1,m2-1,m3+1,2);
        }
      }
      // r_ij, i in mu=3
      rr[iSite][0][0]=0;              // same FCC mu=0 NN
      rr[iSite][0][1]=-1.0/sqrt(2);
      rr[iSite][0][2]=-1.0/sqrt(2);

      rr[iSite][1][0]=0;              // different FCC mu=0 NN
      rr[iSite][1][1]=1.0/sqrt(2);
      rr[iSite][1][2]=1.0/sqrt(2);

      rr[iSite][2][0]=1.0/sqrt(2);    // same FCC mu=1 NN
      rr[iSite][2][1]=0;
      rr[iSite][2][2]=-1.0/sqrt(2);

      rr[iSite][3][0]=-1.0/sqrt(2);   // different FCC mu=1 NN
      rr[iSite][3][1]=0;
      rr[iSite][3][2]=1.0/sqrt(2);

      rr[iSite][4][0]=1.0/sqrt(2);   // same FCC mu=2 NN
      rr[iSite][4][1]=-1.0/sqrt(2);
      rr[iSite][4][2]=0;

      rr[iSite][5][0]=-1.0/sqrt(2);    // different FCC mu=2 NN
      rr[iSite][5][1]=1.0/sqrt(2);
      rr[iSite][5][2]=0;
    }
    else{
      printf("mu =/= {0,1,2,3}. Errore geometry \n");
      exit(EXIT_FAILURE);
    }
  }
}

// Implementazione inizializzazione array ausiliari
void initialize_auxiliary_array(param_t *myParams, int *tempToreplica, int *exchAcc){
  for (int i = 0; i < myParams->n_temp; i++) {
    tempToreplica[i]=i;
  }
  for (int i = 0; i < myParams->n_temp-1; i++) {
    exchAcc[i]=0;
  }
}

// Initialize local axes
void initialize_local_axes(){
  for(int mu=0;mu<NUM_SUBLATTICES;mu++){
    // componente x
    if((mu==0)||(mu==3)){
      local_x[mu][0]=-1.0/sqrt(6);
      local_y[mu][0]=1.0/sqrt(2);
      local_z[mu][0]=1.0/sqrt(3);
    }
    else if((mu==1)||(mu==2)){
      local_x[mu][0]=1.0/sqrt(6);
      local_y[mu][0]=-1.0/sqrt(2);
      local_z[mu][0]=-1.0/sqrt(3);
    }
    // componente y
    if((mu==0)||(mu==2)){
      local_x[mu][1]=-1.0/sqrt(6);
      local_y[mu][1]=-1.0/sqrt(2);
      local_z[mu][1]=1.0/sqrt(3);
    }
    else if((mu==1)||(mu==3)){
      local_x[mu][1]=1.0/sqrt(6);
      local_y[mu][1]=1.0/sqrt(2);
      local_z[mu][1]=-1.0/sqrt(3);
    }
    // componente z
    if((mu==0)||(mu==1)){
      local_x[mu][2]=2.0/sqrt(6);
      local_y[mu][2]=0;
      local_z[mu][2]=1.0/sqrt(3);
    }
    else if((mu==2)||(mu==3)){
      local_x[mu][2]=-2.0/sqrt(6);
      local_y[mu][2]=0;
      local_z[mu][2]=-1.0/sqrt(3);
    }
  }
}

// Initialize DM vectors
void initialize_DM_vectors(){
  DM_vectors[0][1][0]=-1;
  DM_vectors[0][1][1]=1;
  DM_vectors[0][1][2]=0;

  DM_vectors[0][2][0]=1;
  DM_vectors[0][2][1]=0;
  DM_vectors[0][2][2]=-1;

  DM_vectors[0][3][0]=0;
  DM_vectors[0][3][1]=-1;
  DM_vectors[0][3][2]=1;

  DM_vectors[1][2][0]=0;
  DM_vectors[1][2][1]=1;
  DM_vectors[1][2][2]=1;

  DM_vectors[1][3][0]=-1;
  DM_vectors[1][3][1]=0;
  DM_vectors[1][3][2]=-1;

  DM_vectors[2][3][0]=1;
  DM_vectors[2][3][1]=1;
  DM_vectors[2][3][2]=0;

  for(int mu=0;mu<NUM_SUBLATTICES;mu++){
    for(int nu=0;nu<NUM_SUBLATTICES;nu++){
      if(mu==nu){
        for(int j=0;j<NUM_COMPS;j++){
          DM_vectors[mu][nu][j]=0;
        }
      }
      else if(mu>nu){
        for(int j=0;j<NUM_COMPS;j++){
          DM_vectors[mu][nu][j]=-DM_vectors[nu][mu][j];
        }
      }
    }
  }
}

// Initialize v vectors
void initialize_v_vectors(){
  // v0
  for(int i=0;i<NUM_COMPS;i++){
    for(int j=0;j<NUM_COMPS;j++){
      if(i<j){
        for(int k=0;k<NUM_COMPS;k++){
          if(k==i){
            vv[0][i][j][k]=1.0/sqrt(2);
            vv[0][j][i][k]=1.0/sqrt(2);
          }
          else if(k==j){
            vv[0][i][j][k]=1.0/sqrt(2);
            vv[0][j][i][k]=1.0/sqrt(2);
          }
          else{
            vv[0][i][j][k]=0;
            vv[0][j][i][k]=0;
          }
        }
      }
      else if(i==j){
        for(int mu=0;mu<NUM_SUBLATTICES;mu++){
          for(int k=0;k<NUM_COMPS;k++){
            vv[mu][i][j][k]=0.0;
          }
        }
      }
    }
  }
  // v1
  for(int i=0;i<NUM_COMPS-1;i++){
    for(int k=0;k<NUM_COMPS;k++){
      if(k==i){
        vv[1][i][2][2]=1.0/sqrt(2);
        vv[1][2][i][2]=1.0/sqrt(2);
        vv[1][i][2][k]=-vv[1][2][i][2];
        vv[1][2][i][k]=-vv[1][2][i][2];
      }
      vv[1][i][i][k]=0.0;
    }
  }
  for(int k=0;k<NUM_COMPS;k++){
    vv[1][0][1][k]=-vv[0][0][1][k];
    vv[1][1][0][k]=-vv[0][0][1][k];
  }
  // v2
  for(int k=0;k<NUM_COMPS;k++){
    vv[2][0][1][k]=-pow(-1,k)*1.0/sqrt(2);
    vv[2][1][0][k]=-pow(-1,k)*1.0/sqrt(2);
    if(k==NUM_COMPS-1){
      vv[2][0][1][k]=0.0;
      vv[2][1][0][k]=0.0;
    }
    vv[2][0][2][k]=-vv[0][0][2][k];
    vv[2][2][0][k]=-vv[0][0][2][k];
    vv[2][1][2][k]=-vv[1][1][2][k];
    vv[2][2][1][k]=-vv[1][1][2][k];
    for(int i=0;i<NUM_COMPS;i++){
      vv[2][i][i][k]=0.0;
    }
  }
  // v3
  for(int k=0;k<NUM_COMPS;k++){
    vv[3][0][1][k]=-vv[2][0][1][k];
    vv[3][1][0][k]=-vv[2][0][1][k];
    vv[3][0][2][k]=-vv[1][0][2][k];
    vv[3][2][0][k]=-vv[1][0][2][k];
    vv[3][1][2][k]=-vv[0][1][2][k];
    vv[3][2][1][k]=-vv[0][1][2][k];
    for(int i=0;i<NUM_COMPS;i++){
      vv[3][i][i][k]=0.0;
    }
  }
}

// Libreria per implementare le funzioni di inizalizzazione
// Le funzioni di gestione input/output sono trattate in un altro file
dsfmt_t * initialize_system(chain_t **myChains, param_t *myParams, obs_t **myObs, filenames_t **myFilenames, int **tempToreplica, int **exchAcc, dsfmt_t *dsfmt, char const *finput, chain_t *zero_T_replica){
  // Se DEBUG è definito allora resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  // 1. Leggo il file di input
  read_from_input_Param(myParams, finput);
  // 2. Alloco la memoria
  nn = malloc(sizeof(int[myParams->N][NUM_NEAREST_NEIGHBOR]));
  rr = malloc(sizeof(double[myParams->N][NUM_NEAREST_NEIGHBOR][NUM_COMPS]));
  zero_T_replica->lattice = malloc(sizeof(double[myParams->N][NUM_COMPS]));
  *tempToreplica = (int *)malloc((myParams->n_temp)*sizeof(int));
  *exchAcc = (int *)malloc((myParams->n_temp-1)*sizeof(int));
  *myObs = (obs_t *)malloc((myParams->n_temp)*sizeof(obs_t));
  *myChains = (chain_t *)malloc((myParams->n_temp)*sizeof(chain_t));
  *myFilenames = (filenames_t *)malloc((myParams->n_temp)*sizeof(filenames_t));
  for(int i=0;i<myParams->n_temp;i++){
    (*myChains +i)->lattice = malloc(sizeof(double[myParams->N][NUM_COMPS]));
    (*myObs +i)->ene_density=0.0;
    (*myObs +i)->mA2=0.0;
    (*myObs +i)->mSub=0.0;
    for(int j=0;j<NUM_COMPS;j++){
      (*myObs +i)->mT1_par[j]=0.0;
      (*myObs +i)->mT2[j]=0.0;
      (*myObs +i)->mT1_perp[j]=0.0;
    }
    for(int j=0;j<NUM_COMPS-1;j++){
      (*myObs +i)->mE[j]=0.0;
    }
  }
  initialize_auxiliary_array(myParams, *tempToreplica, *exchAcc);
  // 3. PBC
  geometry(myParams);
  // 4. Generatore random
  int rank_max = omp_get_max_threads();                     // Numero massimo di threads
  dsfmt = (dsfmt_t *)malloc(rank_max * sizeof(dsfmt_t));    // Alloco dinamicamente la memoria
  dsfmt_init_gen_rand(&dsfmt[0], myParams->dSFMT_seed);     // Inizializzo a partire dall'indice di thread 0
  for(int i=1;i<rank_max;i++){
    dsfmt[i]=dsfmt[i-1];
    dSFMT_jump(&dsfmt[i], jump);
  }
  // 5. Inizializzo gli assi z locali ed i vettori DM
  initialize_local_axes();
  initialize_DM_vectors();
  initialize_v_vectors();
  // 6. Inizializzo a seconda di iStart
  if(myParams->iStart==2){
    // Leggo da file precedente
    read_from_filenames(*myFilenames, myParams);
    read_Status(myParams, *myChains, *myFilenames, *tempToreplica);
  }
  else if(myParams->iStart==1){           // Partenza a freddo (ordinata)
    for(int i=0;i<myParams->n_temp;i++){
      // Azzero acc, count e definisco le temperature di partenza equispaziate
      (*myChains + i)->acc=0;
      (*myChains + i)->count=0;
      // (*myChains + i)->T=myParams->t_max - (myParams->t_max -myParams->t_min)*i/(myParams->n_temp -1);
      (*myChains + i)->T=myParams->t_max * pow(myParams->t_min/myParams->t_max, (double)i/(myParams->n_temp-1));
      (*myChains + i)->acc_Parameter = 0.5;
      for(int j=0;j<myParams->N;j++){
        int mu=(j/((int)pow(myParams->size,3)));
        for(int k=0;k<NUM_COMPS;k++){     // Ordino AIAO
          (*myChains + i)->lattice[j][k]=local_z[mu][k];
        }
      }
    }
  }
  else if(myParams->iStart==0){           // Partenza a caldo (disordinata)
    for(int i=0;i<myParams->n_temp;i++){
      // Azzero acc, count e definisco le temperature di partenza equispaziate
      (*myChains + i)->acc=0;
      (*myChains + i)->count=0;
      (*myChains + i)->T=myParams->t_max - (myParams->t_max -myParams->t_min)*i/(myParams->n_temp -1);
      (*myChains + i)->acc_Parameter = 0.5;
      for(int j=0;j<myParams->N;j++){
        for(int k=0;k<NUM_COMPS;k++){     // Componenti random tra [-1;1) che dopo normalizzo
          (*myChains + i)->lattice[j][k]=2*dsfmt_genrand_close_open(&dsfmt[omp_get_num_threads()])-1.0;
        }
      }
    }
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0;i<myParams->n_temp;i++){  // Rinormalizzo gli array random
        renormalize_spin((*myChains +i), myParams);
      }
      #pragma omp barrier
    }
  }
  return dsfmt;
}

// Procedura che controlla le accettanze Metropolis
void check_metro_acc(chain_t *single_chain, param_t *myParams){
  double acc;
  double step = 0.005;
  acc=(double)single_chain->acc/single_chain->count;
  if(acc>acc_max){
    single_chain->acc_Parameter=single_chain->acc_Parameter + step;
    if(single_chain->acc_Parameter>1.0){
      single_chain->acc_Parameter=0.95;
    }
  }
  else if(acc<acc_min){
    single_chain->acc_Parameter=single_chain->acc_Parameter - step;
    if(single_chain->acc_Parameter<0.01){
      single_chain->acc_Parameter=0.3;
    }
  }
  single_chain->acc=0;
  single_chain->count=0;
}

// Procedura che aggiusta le temperature
// In questa procedura aggiusto anche i parametri di accettanza
void set_temperatures(chain_t *myChains, param_t *myParams, dsfmt_t *dsfmt, int *tempToreplica, int *exchAcc){
  // Array necessari alla procedura
  double g[myParams->n_temp], auxEnergies[myParams->n_temp];
  double beta_prime;
  int rank;
  // L'obietttivo finale è avere più o meno
  // un'accettanza piatta su tutto il sistema
  for(int t=0;t<myParams->iSet;t++){
    // Azzero il contatore degli swap
    for(int i=0;i<myParams->n_temp-1;i++){
      exchAcc[i]=0;
    }
    // Resetto gli array ausiliari
    for(int i=0;i<myParams->n_temp;i++){
      g[i]=0;
      auxEnergies[i]=0;
    }
    // Termalizzazione
    // Prima monte carlo indipendenti in parallelo
    // Poi swap
    for(int j=0;j<myParams->iPT;j++){
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          rank=omp_get_thread_num();
          update_metropolis(myParams, (myChains + tempToreplica[i]), dsfmt, rank);
        }
        #pragma omp barrier
      }
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          check_metro_acc((myChains + tempToreplica[i]), myParams);
        }
        #pragma omp barrier
      }
      // Scambio le temperature
      swap_temperatures(myParams, myChains, dsfmt, rank, tempToreplica, exchAcc);
      // Calcolo le energie medie
      // In parallelo
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          auxEnergies[i]+=get_energy((myChains + tempToreplica[i]), myParams)/myParams->iPT;
        }
        #pragma omp barrier
      }
    }
    // Costruisco il vettore g
    for(int i=1;i<myParams->n_temp;i++){
      if(auxEnergies[i-1]==auxEnergies[i+1]){
        // Energie uguali
        g[i] = 1.0/(myChains+tempToreplica[i])->T;
      }
      else{
        g[i] = (auxEnergies[i-1]/(myChains+tempToreplica[i-1])->T - auxEnergies[i+1]/(myChains+tempToreplica[i+1])->T - auxEnergies[i]*(1.0/(myChains+tempToreplica[i-1])->T - 1.0/(myChains+tempToreplica[i+1])->T)) /(auxEnergies[i-1]-auxEnergies[i+1]);
      }
    }
    if((t%2)==0){
      // Cambio le temperature con indice pari tenendo le dispari fisse
      for(int i=2;i<myParams->n_temp-1;i=i+2){
        beta_prime = (1.0/((myChains+tempToreplica[i])->T) + g[i])/2.0;
        if(((1.0/beta_prime)<myParams->t_min)||((1.0/beta_prime)>myParams->t_max)){
          // DOES NOTHING
        }
        else{
          myChains[tempToreplica[i]].T = 1.0/beta_prime;
        }
      }
    }
    else if((t%2)==1){
      for(int i=1;i<myParams->n_temp-1;i=i+2){
        beta_prime = (1.0/((myChains + tempToreplica[i])->T) + g[i])/2.0;
        if(((1.0/beta_prime)<myParams->t_min)||((1.0/beta_prime)>myParams->t_max)){
          // DOES NOTHING
        }
        else{
          myChains[tempToreplica[i]].T = 1.0/beta_prime;
        }
      }
    }
    else{
      printf("Errore in set_temperatures\n");
      exit(EXIT_FAILURE);
    }
  }
  // CHEK IF EVERY TEMPERATURE IS IN GIVEN RANGE
  for(int i=1;i<myParams->n_temp-1;i++){
    if(((myChains+tempToreplica[i])->T>myParams->t_max)||((myChains+tempToreplica[i])->T<myParams->t_min)){
      printf("Errore in set_temperatures Tmin=%lf, Tmax=%lf , T=%lf\n", myParams->t_min, myParams->t_max, (myChains+tempToreplica[i])->T);
      exit(EXIT_FAILURE);
    }
  }
}

// Procedura per le accettanze
// Da chiamare una volta fissate le temperature
void set_acceptances(chain_t *myChains, param_t *myParams, dsfmt_t *dsfmt, int *tempToreplica){
  int rank;
  for(int j=0;j<myParams->iTerm;j++){
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0;i<myParams->n_temp;i++){
        rank=omp_get_thread_num();
        update_metropolis(myParams, (myChains + tempToreplica[i]), dsfmt, rank);
      }
      #pragma omp barrier
    }
    #pragma omp parallel
    {
      #pragma omp for
      for(int k=0;k<myParams->n_temp;k++){
        check_metro_acc((myChains + tempToreplica[k]), myParams);
      }
      #pragma omp barrier
    }
  }
}

// Procedura che crea i nomi dei file su cui scrivo
void set_filenames(chain_t *myChains, param_t *myParams, filenames_t *myFilenames, int *tempToreplica){
  for(int i=0;i<myParams->n_temp;i++){
    // Creo le stringhe che per i nomi dei file
    sprintf((myFilenames+i)->data_file, "data/data_%d_%lf.dat", myParams->size, (myChains+tempToreplica[i])->T);
    sprintf((myFilenames+i)->backup_file, "backup/backup_%d_%lf.bin", myParams->size, (myChains+tempToreplica[i])->T);
    sprintf((myFilenames+i)->aux_file, "aux/aux_%d_%lf.bin", myParams->size, (myChains+tempToreplica[i])->T);
    if(myParams->iSave!=0){
      sprintf((myFilenames+i)->conf_file, "conf/conf_%d_%lf.txt", myParams->size, (myChains+tempToreplica[i])->T);
    }
  }
  // Scrvo sui file dei dati e delle configurazioni l'header
  for(int i=0;i<myParams->n_temp;i++){
    FILE *fptr1;
    FILE *fptr2;
    fptr1=fopen((myFilenames+i)->data_file, "w");
    if(fptr1==NULL){
      perror("Errore in apertura dei file in set_filenames");
      exit(EXIT_FAILURE);
    }
    else{
      fprintf(fptr1, "%d\t%d\t%d\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", myParams->size, myParams->N, NUM_COMPS, (myChains+tempToreplica[i])->T, myParams->heis, myParams->ising, myParams->PD, myParams->DM, myParams->single_ion);
      fprintf(fptr1, "Ene_dens\tmAIAO\tmx\tmy\tmz\n\n");
    }
    fclose(fptr1);
    if(myParams->iSave!=0){
      fptr2=fopen((myFilenames+i)->conf_file, "w");
      if(fptr2==NULL){
        perror("Errore in apertura dei file in set_filenames");
        exit(EXIT_FAILURE);
      }
      else {
        fprintf(fptr2, "%d\t%d\t%d\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", myParams->size, myParams->N, NUM_COMPS, (myChains+tempToreplica[i])->T, myParams->heis, myParams->ising, myParams->PD, myParams->DM, myParams->single_ion);
        fprintf(fptr2, "iSite\tSx\tSy\tSz\n\n");
      }
      fclose(fptr2);
    }
  }
}

// Procedura per la prima termalizzazione
void burn_in(chain_t *myChains, param_t *myParams, dsfmt_t *dsfmt, int *tempToreplica, int *exchAcc){
  int rank;
  // Resetto le accettanze
  for(int i=0;i<myParams->n_temp-1;i++){
    exchAcc[i]=0;
  }
  // Cicli a vuoto
  for(int j=0;j<myParams->iTerm;j++){
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0;i<myParams->n_temp;i++){
        rank=omp_get_thread_num();
        update_metropolis(myParams, (myChains + tempToreplica[i]), dsfmt, rank);
      }
      #pragma omp barrier
    }
    swap_temperatures(myParams, myChains, dsfmt, rank, tempToreplica, exchAcc);
  }
}

// Procedura per la deallocazione della memoria
void deallocation(chain_t **myChains, param_t *myParams, obs_t **myObs, filenames_t **myFilenames, int **tempToreplica, int **exchAcc, dsfmt_t *dsfmt, chain_t *zero_T_replica){
  free(nn);
  free(rr);
  free(*tempToreplica);
  free(*myObs);
  free(*exchAcc);
  free(dsfmt);
  free(*myFilenames);
  for(int i=0;i<myParams->n_temp;i++){
    free((*myChains +i)->lattice);
  }
  free(*myChains);
  free(zero_T_replica->lattice);
}
