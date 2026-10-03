// Includo l'header
#include "../head/header.h"

// Procedura per eliminare spazi e commenti dal file di input
void remove_white_line_and_comments(FILE *input){
  int temp_i;

  temp_i=getc(input);
  if(temp_i=='\n' || temp_i==' ' || temp_i=='\043') // SCAN PER SPAZI VUOTI E COMMENTI
    {
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\n' || temp_i==' ') // LINEE VUOTE
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i=='\n' || temp_i==' ');
      }
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\043')  // COMMENTI, 043 È IL CODICE ASCII PER #
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i!='\n');
      }
    else
      {
      ungetc(temp_i, input);
      }

    remove_white_line_and_comments(input);
    }
  else
    {
    ungetc(temp_i, input);
    }
  }

// Procedura che legge da file i parametri di sistema
void read_from_input_Param(param_t *Par, char const *finput){
  int flag, temp_i;
  int end=1;
  FILE *fInput;
  char str[100];

  fInput = fopen(finput, "r");
  if (fInput == NULL) {
    perror("Errore in apertura in read_from_input_Param");
    exit(EXIT_FAILURE);
  }
  else{
    while(end==1){
      remove_white_line_and_comments(fInput);

      flag = fscanf(fInput, "%s", str);
      if (flag != 1) {
        perror("Errore di lettura della prima stringa");
        exit(EXIT_FAILURE);
      }
      if(strncmp(str, "size", 4) == 0){
        flag = fscanf(fInput, "%d", &(Par->size));
        if (flag != 1){
          perror("Errore di lettura nella taglia");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iStart", 6) == 0){
        flag = fscanf(fInput, "%d", &(Par->iStart));
        if (flag != 1){
          perror("Errore di lettura nella temperatura");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "T_min", 5) == 0){
        flag = fscanf(fInput, "%lf", &(Par->t_min));
        if (flag != 1){
          perror("Errore di lettura nella temperatura");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "T_max", 5) == 0){
        flag = fscanf(fInput, "%lf", &(Par->t_max));
        if (flag != 1){
          perror("Errore di lettura in alpha");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "heis", 4) == 0){
        flag = fscanf(fInput, "%lf", &(Par->heis));
        if (flag != 1){
          perror("Errore di lettura in eps");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "ising", 5) == 0){
        flag = fscanf(fInput, "%lf", &(Par->ising));
        if (flag != 1){
          perror("Errore di lettura in eps");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "pseudodipolar", 13) == 0){
        flag = fscanf(fInput, "%lf", &(Par->PD));
        if (flag != 1){
          perror("Errore di lettura in eps");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "DM", 2) == 0){
        flag = fscanf(fInput, "%lf", &(Par->DM));
        if (flag != 1){
          perror("Errore di lettura in eps");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "single_ion", 10) == 0){
        flag = fscanf(fInput, "%lf", &(Par->single_ion));
        if (flag != 1){
          perror("Errore di lettura in eps");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iSet", 4) == 0){
        flag = fscanf(fInput, "%d", &(Par->iSet));
        if (flag != 1){
          perror("Errore di lettura in iTerm");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iPT", 3) == 0){
        flag = fscanf(fInput, "%d", &(Par->iPT));
        if (flag != 1){
          perror("Errore di lettura in iDec");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iTerm", 5) == 0){
        flag = fscanf(fInput, "%d", &(Par->iTerm));
        if (flag != 1){
          perror("Errore di lettura in iMis");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iMis", 4) == 0){
        flag = fscanf(fInput, "%d", &(Par->iMis));
        if (flag != 1){
          perror("Errore di lettura in iOverr");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iDec", 4) == 0){
        flag = fscanf(fInput, "%d", &(Par->iDec));
        if (flag != 1){
          perror("Errore di lettura in iStart");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iMinimize", 9) == 0){
        flag = fscanf(fInput, "%d", &(Par->iMinimize));
        if (flag != 1){
          perror("Errore di lettura in iStart");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iSave", 5) == 0){
        flag = fscanf(fInput, "%d", &(Par->iSave));
        if (flag != 1){
          perror("Errore di lettura in iSave");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "iBackup", 7) == 0){
        flag = fscanf(fInput, "%d", &(Par->iBackup));
        if (flag != 1){
          perror("Errore di lettura in iBackup");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "filenames_file", 14) == 0){
        flag = fscanf(fInput, "%s", Par->filenames_file);
        if (flag != 1){
          perror("Errore di lettura in filenames_file");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "log_file", 8) == 0){
        flag = fscanf(fInput, "%s", Par->log_file);
        if (flag != 1){
          perror("Errore di lettura in filenames_file");
          exit(EXIT_FAILURE);
        }
      }
      else if(strncmp(str, "random_seed", 11) == 0){
        flag = fscanf(fInput, "%u", &(Par->dSFMT_seed));
        if (flag != 1){
          perror("Errore di lettura in random_seed");
          exit(EXIT_FAILURE);
        }
      }
      else{
        fprintf(stderr, "Error: unrecognized option %s in the file %s (%s, %d)\n", str, finput, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }

      remove_white_line_and_comments(fInput);

      // CONTROLLO SE LA LINEA LETTA È L'ULTIMA
      temp_i = getc(fInput);
      if(temp_i == EOF){
        end=0;
      }
      else{
        ungetc(temp_i, fInput);
      }
    }

    fclose(fInput);
    Par->N=NUM_SUBLATTICES*(int)(pow(Par->size, 3));
    Par->n_temp = (int)(2*(sqrt((double)Par->N))+1);    // Numero di temperature
    if(Par->n_temp<=20){
      Par->n_temp=20;
    }
  }
}

// Procedura che legge i nomi dei file già salvati
void read_from_filenames(filenames_t *myFilenames, param_t *myParams){
  FILE *fFilenames;
  int flag;
  // Apro il file dove sono scritti i
  // file da simulazione precedente
  fFilenames=fopen(myParams->filenames_file, "r");
  if (fFilenames == NULL) {
    perror("Errore in apertura in read_from_filenames");
    exit(EXIT_FAILURE);
  }
  else{
    // Salvo i nomi dei file
    for(int i=0;i<myParams->n_temp;i++){
      flag = fscanf(fFilenames, "%s", (myFilenames+i)->data_file);
      if (flag != 1) {
        perror("Errore di lettura della prima stringa");
        exit(EXIT_FAILURE);
      }
      flag = fscanf(fFilenames, "%s", (myFilenames+i)->backup_file);
      if (flag != 1) {
        perror("Errore di lettura della prima stringa");
        exit(EXIT_FAILURE);
      }
      flag = fscanf(fFilenames, "%s", (myFilenames+i)->aux_file);
      if (flag != 1) {
        perror("Errore di lettura della prima stringa");
        exit(EXIT_FAILURE);
      }
      if(myParams->iSave!=0){
        flag = fscanf(fFilenames, "%s", (myFilenames+i)->conf_file);
        if (flag != 1) {
          perror("Errore di lettura della prima stringa");
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  fclose(fFilenames);
}

// Procedura che scrive su file lo stato del sistema
void write_Status(param_t *myParams, chain_t *myChains, filenames_t *myFilenames, int *tempToreplica){
  FILE *fConf;
  FILE *fAux;
  FILE *fFilenames;
  // Ho bisogno di una struttura copia
  typedef struct {
    double lattice[myParams->N][NUM_COMPS];
    double T;
    int acc;
    int count;
    double acc_Parameter;
  } cpy_chain_t;

  cpy_chain_t *cpyChains;

  cpyChains = (cpy_chain_t *)malloc((myParams->n_temp)*sizeof(cpy_chain_t));

  // Copio lo Status
  for(int i=0;i<myParams->n_temp;i++){
    (cpyChains+i)->T = myChains[tempToreplica[i]].T;
    (cpyChains+i)->acc = myChains[tempToreplica[i]].acc;
    (cpyChains+i)->count = myChains[tempToreplica[i]].count;
    (cpyChains+i)->acc_Parameter = myChains[tempToreplica[i]].acc_Parameter;
    for(int j=0;j<myParams->N;j++){
      for(int k=0;k<NUM_COMPS;k++){
        (cpyChains+i)->lattice[j][k] = myChains[tempToreplica[i]].lattice[j][k];
      }
    }
  }
  // Salvo la configurazione
  // Quando la leggo è bene riazzerare acc e count
  for(int i=0;i<myParams->n_temp;i++){
    fConf = fopen(myFilenames[i].backup_file, "wb");
    fAux = fopen(myFilenames[i].aux_file, "wb");
    fFilenames = fopen(myParams->filenames_file, "a");
    if((fConf == NULL) || (fAux == NULL) || (fFilenames == NULL)){
      perror("Errore in apertura write_Status");
      exit(EXIT_FAILURE);
    }
    else {
      fwrite(&cpyChains[i], sizeof(cpy_chain_t), 1, fConf);
      fwrite(&tempToreplica[i], sizeof(int), 1, fAux);
      fprintf(fFilenames, "%s\t%s\t%s\t%s\n", myFilenames[i].data_file, myFilenames[i].backup_file, myFilenames[i].aux_file, myFilenames[i].conf_file);
    }
    fclose(fConf);
    fclose(fAux);
    fclose(fFilenames);
  }
  free(cpyChains);
}

// Procedura che legge da file l'ultimo stato salvato del sistema
void read_Status(param_t *myParams, chain_t *myChains, filenames_t *myFilenames, int *tempToreplica){
  FILE *fConf;
  FILE *fAux;
  // Ho bisogno di una struttura copia
  typedef struct {
    double lattice[myParams->N][NUM_COMPS];
    double T;
    int acc;
    int count;
    double acc_Parameter;
  } cpy_chain_t;

  cpy_chain_t *cpyChains;

  cpyChains = (cpy_chain_t *)malloc((myParams->n_temp)*sizeof(cpy_chain_t));

  // Leggo la configurazione da file
  // Dopo averla letta riazzero acc e count
  for(int i=0;i<myParams->n_temp;i++){
    fConf = fopen(myFilenames[i].backup_file, "rb");
    fAux = fopen(myFilenames[i].aux_file, "rb");
    if((fConf == NULL)||(fAux == NULL)){
      perror("Errore in apertura read_Status");
      exit(EXIT_FAILURE);
    }
    else {
      int k = fread(&tempToreplica[i], sizeof(int), 1, fAux);
      int j = fread(&cpyChains[i], sizeof(cpy_chain_t), 1, fConf);
      if((j!=1)||(k!=1)){
        perror("errore in fread nella procedura read_Status");
        exit(EXIT_FAILURE);
      }
    }
    fclose(fConf);
    fclose(fAux);
  }

  // La configurazione copiata viene inserita nella memoria dedicata
  for(int i=0;i<myParams->n_temp;i++){
    myChains[tempToreplica[i]].T = (cpyChains+i)->T;
    myChains[tempToreplica[i]].acc = 0;
    myChains[tempToreplica[i]].count = 0;
    myChains[tempToreplica[i]].acc_Parameter = (cpyChains+i)->acc_Parameter;
    for(int j=0;j<myParams->N;j++){
      for(int k=0;k<NUM_COMPS;k++){
        (myChains + tempToreplica[i])->lattice[j][k] = (cpyChains+i)->lattice[j][k];
      }
    }
  }
  free(cpyChains);
}

// Procedura che salva le configurazioni di spin
void save_Configs(chain_t *myChains, param_t *myParams, filenames_t *myFilenames, int *tempToreplica){
  FILE *fConf;
  // Salvo le configurazioni
  for(int i=0;i<myParams->n_temp;i++){
    fConf=fopen(myFilenames[i].conf_file, "a");
    if((fConf == NULL)){
      perror("Errore in apertura save_Configs");
      exit(EXIT_FAILURE);
    }
    else{
      for(int j=0;j<myParams->N;j++){
        fprintf(fConf, "%d\t%.13lf\t%.13lf\t%.13lf\n", j, myChains[tempToreplica[i]].lattice[j][0], myChains[tempToreplica[i]].lattice[j][1], myChains[tempToreplica[i]].lattice[j][2]);
      }
    }
    fclose(fConf);
  }
}

// Procedura che salva la configurazione a temperatura più bassa
void save_zeroT_config(param_t *myParams, chain_t *single_chain){
  FILE *fptr;
  FILE *fptr2;
  char conf_file[STD_STRING_LENGTH];
  char measure_file[STD_STRING_LENGTH];
  double ene_dens=0;
  double projection=0;
  double magn[NUM_COMPS];
  for(int i=0;i<NUM_COMPS;i++){
    magn[i]=0.0;
  }
  // Creo il nome del file
  sprintf(conf_file, "zero_T_conf/conf_L%d.txt", myParams->size);
  sprintf(measure_file, "zero_T_conf/data_L%d.txt", myParams->size);
  // Apro il file
  fptr = fopen(conf_file, "w");
  fptr2= fopen(measure_file, "w");
  if((fptr==NULL) || (fptr2==NULL)){
    perror("Errore in apertura in save_zeroT_config");
    exit(EXIT_FAILURE);
  }
  else{
    fprintf(fptr, "Sito\tSx\tSz\tSz\n\n");
    fprintf(fptr2, "Ene_dens\tmagn_modulus\n");
    for(int i=0;i<myParams->N;i++){
      for(int j=0;j<NUM_COMPS;j++){
        magn[j]+=single_chain->lattice[i][j]/myParams->N;
      }
      fprintf(fptr, "%d\t%.13lf\t%.13lf\t%.13lf\n", i, single_chain->lattice[i][0], single_chain->lattice[i][1], single_chain->lattice[i][2]);
    }
    ene_dens=get_energy(single_chain, myParams)/myParams->N;
    fprintf(fptr2, "%lf\t%lf\n", ene_dens, sqrt(spin_product(magn, magn)));
  }
  fclose(fptr);
  fclose(fptr2);
}

// Procedura per scrivere il file di log
void write_Log(param_t *myParams, chain_t *myChains, int *tempToreplica, int *exchAcc){
  FILE *fptr;
  int rank_max=omp_get_max_threads();

  fptr=fopen(myParams->log_file, "a");
  if(fptr==NULL){
    perror("Errore in apertura write_Log");
    exit(EXIT_FAILURE);
  }
  // Testa del log
  fprintf(fptr, "+----------------------------------------------------------+\n");
  fprintf(fptr, "| Simulation details for pyro_kramers w Parallel Tempering |\n");
  fprintf(fptr, "+----------------------------------------------------------+\n\n");
  // Debug attivo
  #ifdef DEBUG
  fprintf(fptr, "DEBUG attivo\n\n");
  #endif
  // Caratteristiche geometriche
  fprintf(fptr, "PYROCHLORE LATTICE WITH FCC BRAVAIS LATTICE AND 4 SUBLATTICES.\n");
  fprintf(fptr, "NUMBER OF FCC POINTS PER SIDE %d\n", myParams->size);
  fprintf(fptr, "TOTAL SITES NUMBER %d\n", myParams->N);
  fprintf(fptr, "SPIN COMPONENTS %d\n\n", NUM_COMPS);
  // Partenza
  fprintf(fptr, "STARTING INDEX %d (0 disordered, 1 ordered, 2 from previous file)\n\n", myParams->iStart);
  // Minimizzazione iterativa effettuata o meno
  if(myParams->iMinimize==1){
    fprintf(fptr, "INDEX iMinimize %d, ITERATIVE MINIMIZATION DONE\n\n", myParams->iMinimize);
  }
  else{
    fprintf(fptr, "MINIMIZATION NOT DONE\n\n");
  }
  // Salvataggio di configurazioni per spin static structure factor
  if(myParams->iSave!=0){
    fprintf(fptr, "INDEX iSave %d, %d SPIN CONFIGURATIONS SAVED\n\n", myParams->iSave, myParams->iMis/myParams->iSave);
  }
  else{
    fprintf(fptr, "iSave 0, CONFIGURATIONS NOT SAVED\n\n");
  }
  // Accoppiamenti
  fprintf(fptr, "HEISENBERG COUPLING %lf\n", myParams->heis);
  fprintf(fptr, "ISING COUPLING %lf\n", myParams->ising);
  fprintf(fptr, "PSEUDODIPOLA COUPLING %lf\n", myParams->PD);
  fprintf(fptr, "DM COUPLING %lf\n", myParams->DM);
  fprintf(fptr, "SINGLE ION %lf\n", myParams->single_ion);

  // Indici vari
  fprintf(fptr, "NUMBER OF ITERATIONS TO FIX TEMPERATURES %d\n", myParams->iSet);
  fprintf(fptr, "NUMBER OF SWAPS FOR EACH ITERATION %d\n\n", myParams->iPT);
  // Indice di burn_in
  fprintf(fptr, "BURN IN LOOPS %d\n\n", myParams->iTerm);
  // Numero di misure
  fprintf(fptr, "NUMBER OF MEASURES %d\n", myParams->iMis);
  fprintf(fptr, "NUMBER OF DECORRELATION LOOPS BETWEEN MEASURES %d\n\n", myParams->iDec);
  // Numero di misure tra un backup e l'altro
  fprintf(fptr, "NUMBER OF MEASURES BETWEEN TWO BACKUPS %d\n\n", myParams->iBackup);
  // Temperatura minima, massima e numero di temperature
  fprintf(fptr, "MINIMUM TEMPERATURE %lf\n", myParams->t_min);
  fprintf(fptr, "MAXIMUM TEMPERATURE %lf\n", myParams->t_max);
  fprintf(fptr, "NUMBER OF TEMPERATURES %d\n\n", myParams->n_temp);
  // Accettanze
  fprintf(fptr, "T\t\tmetro_acc\t\tPTacc\t\tmetro_eps\n");
  // Scrivo le caratteristiche per ogni temperatura
  for(int i=0;i<myParams->n_temp-1;i++){
    fprintf(fptr, "%lf\t%lf\t%lf\t%lf\n", (myChains+tempToreplica[i])->T, (double)(myChains+tempToreplica[i])->acc/(myChains+tempToreplica[i])->count, (double)exchAcc[i]/(myParams->iMis * myParams->iDec), (myChains+tempToreplica[i])->acc_Parameter);
  }
  fprintf(fptr, "%lf\t%lf\t--------\t%lf\n\n", (myChains+tempToreplica[myParams->n_temp-1])->T, (double)(myChains+tempToreplica[myParams->n_temp-1])->acc/(myChains+tempToreplica[myParams->n_temp-1])->count, (myChains+tempToreplica[myParams->n_temp-1])->acc_Parameter);
  // Numero di thread usati
  fprintf(fptr, "USED CORE(1 thread 1 core) %d\n\n", rank_max);
  // Calcolo il tempo di esecuzione
  end = clock();
  double exe_time = ((double)end - start)/CLOCKS_PER_SEC;  // Total time
  exe_time=exe_time/rank_max;
  int minutes = (int)(exe_time/60);
  double seconds = exe_time - minutes*60;
  int hours = (int)(minutes/60);
  minutes = minutes - hours*60;
  fprintf(fptr, "EXECUTION TIME PER CPU %d h %d m %lf s\n\n", hours, minutes, seconds);
  // Chiudo il file
  fprintf(fptr, "CORRECT EXECUTION\n\n");
  fclose(fptr);
}
