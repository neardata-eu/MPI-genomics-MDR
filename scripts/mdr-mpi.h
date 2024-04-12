/**************************************
Copyright 2023 Andrés Benavides A and Gonzalo Gómez-Sánchez

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

MPI DATA CONNECTOR 
Connector to process all the samples from the listoffiles.txt file
It is based on using MPI to distribute the load accross machines.

VERSION 1.0 
***************************************/


#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>

/**************************************
 **************PATHS*******************
***************************************/
#define ROOT_PATH "MPI-genomics-MDR"
#define DATA_SET  "/data/simple_10x30x100"
#define DATAPATH ROOT_PATH DATA_SET "/input"
#define OUTPUTPATH ROOT_PATH DATA_SET "/output-mpi-c"
#define LOGPATH ROOT_PATH DATA_SET "/logs"
#define LISTOFFILESPATH DATAPATH "/listoffiles.txt"
#define SAMPLEPATH DATAPATH "/samples/"
#define Npatients  100
#define nsamples 30         // samples per files
#define files2proc  10      // Default value for nfiles
#define CC_RATIO 0.5
#define filter_imp 0.9
#define CV_sets 5

#define MAX_FILENAME_LENGTH 1024

/**************************************
 **************Functions***************
***************************************/

// Read file with labels and keep only cases/controls
void read_labels(const char* labelspath, int* labels) {

    FILE* csv_file = fopen(labelspath, "r");
 
    char line[10*Npatients]; //10 char per patient
    char* token;
    int line_count = 0;
    // Skip first lines
    token=fgets(line, sizeof(line), csv_file);
    token=fgets(line, sizeof(line), csv_file);
    while (fgets(line, sizeof(line), csv_file)) {
        token = strtok(line, "\t");
        labels[line_count] = atoi(&token[strlen(token) - 2]);
        line_count++;
    }

    fclose(csv_file);
    //printf("%s labels read.\n", line);
}

// Generate two set case/ctrl 0->1 and 1->0
int convert_labels(int* cases, int* control) {
    int total_cases=0;
    for (int i = 0; i < Npatients; i++) {
        control[i] = (cases[i] == 0) ? 1 : 0;
        total_cases+=cases[i];
    }
    return total_cases;
}

// Read list of files and convert into a list
void file_to_list(const char* filepath, char* files) {

    FILE* files_list = fopen(filepath, "r");
 
    char line[MAX_FILENAME_LENGTH]; 

    while (fgets(line, sizeof(line), files_list)) {
        strcat(files, line);
    }

    fclose(files_list);
   
}

// Read list of files and convert into list
void get_comb_files(const char* filepath, char* files, int first, int last ) {

    FILE* files_list = fopen(filepath, "r");
    int num_lines=0;
    char line[MAX_FILENAME_LENGTH]; 
    while(num_lines<first)
    {
       char *token = fgets(line, sizeof(line), files_list);
       (void) token;
       num_lines++;
    }

    while (fgets(line, sizeof(line), files_list)) {
        if(num_lines>last-1){
            break;
        }
        strcat(files, line);
        num_lines++;     
    }

    fclose(files_list);
   
}

// Definition of the filter for imputation. We believe the value only if > 0.9
int filter_imputation(float x) {
    return (x > filter_imp) ? 1 : 0;
}

// Transform to key + values
void get_keyval(char *line, char key [nsamples][50], int value[nsamples][Npatients], int iteration){
    
    int i=0, j=0, k=0; int row=0;
    while (line[i]!='\0'){
        if(line[i]==' '){ //change this line to use other separator
            row++;
        }
        else{
            if(row==0){
                key[iteration][k++]=line[i];
            }
            else if(row==1){
                key[iteration][k]='-';
            }
            else if(row==2){
                key[iteration][++k]=line[i];
            }
            else if(row==5){
                break;
            }
        }
        i++;
    }
    key[iteration][++k]='\0'; row=0;
    char data[10]; int v=0; int variante=0;
    while (line[i]!='\0'){
        if(line[i]!=' '){
            data[v++]=line[i];
        }
        else{
            data[v++]='\0';          
            v=0;
            double value_doble = strtod(&data[0], NULL);
            variante +=filter_imputation(value_doble)*(row);
            row++;  
            if(row==3){
                value[iteration][j++]=variante;
                variante=0;
                row=0;
            } 
        }
        i++;
    }
    double value_doble = strtod(&data[0], NULL);
    variante +=filter_imputation(value_doble)*(row);
    value[iteration][j]=variante;
}

// Read sample information and save it into a unordered_map
void read_sample(const char* samplepath, char keys[nsamples][50], int values[nsamples][Npatients]) {

    gzFile file = gzopen(samplepath, "rb");
    char line[MAX_FILENAME_LENGTH+10*Npatients]; 
    int line_count=0;
    gzgets(file, line, sizeof(line)); //skip first line
    while (gzgets(file, line, sizeof(line)) != NULL ) {
        if(line[0] == '#') //skip headers
    	    continue;
        get_keyval(line, keys, values, line_count++);
    }
    gzclose(file);
}


// Transform to key + values
void transformPatients(int patients1 [Npatients], int patients2[Npatients],int* ptcode) {
   
    for (int i=0; i<Npatients;i++) {
        int pt1=(patients1[i]+1)*3;
        int pt2=patients2[i];
        ptcode[i]=(pt1-pt2);       
    }

}

// Apply risk vector to classify the patients
void applyRisk(int* patients,int* risk, int* prediction){
    //casevalues = np.where(risk == 1)
    int casevalues[10];
    int j=0;
    for (int i = 1; i < 10; ++i) {
        if(risk[i]==1){
            casevalues[j]=i;
            j++;
        }
    }
    //prediction[patients==n] = 1
    for(int k=0; k<j;k++){  
        for (int i=0; i<Npatients;i++){ 
            prediction[i]=0;          
            if(patients[i]==casevalues[k]){
                prediction[i]=1;
            }
        }
    }
}

//Count number of cases and number of controls for each column and return the high risk combinations.
//Then, use the high risk predictor to obtain the classification and prediction error.
void getRiskArray(int* patients,float* testError,
                  int trainset [CV_sets] [Npatients], int testset [CV_sets] [Npatients],
                  int* npcases,int* npcontrols,double ccratio){
    
    for (int i = 0; i < CV_sets; ++i) {
        
        int prediction [Npatients];
        int sumCases[10]; 
        int sumControls[10];
        int risk[10];
        int nTrain =0; int nTest =0;
       
        //init sums in zero
        for(int k=0;k<10;++k){
            sumCases[k]=0;
            sumControls[k]=0;
            risk[k]=0;
        }

        // 1 - Get the sets for the iteration
        for (int j=0; j<Npatients;++j){            
            nTrain+=trainset [i] [j];
            nTest+=testset [i] [j];
            prediction[j]=0;
        }

        //2 - Sum the cases controls from training set and the controls from testing set
        for (int j=0; j<Npatients;++j) {
            int caseCounts=patients[j] * npcases[j] * trainset[i][j];
            int controlCounts=patients[j] * npcontrols[j] * trainset[i][j];
            sumCases[caseCounts]+=1;
            sumControls[controlCounts]+=1;
        }

        //3 - Get risk array
        double risk_val=0;
        for (int j = 0; j < 10; ++j) {
            if (sumControls[j] != 0) {
                risk_val = (double)(sumCases[j] / sumControls[j]);
            } else {
                risk_val = 0;
            }
            //4 - Transform to high risk = 1, low risk = 0
            if (risk_val >=  ccratio)
                risk[j] = 1;
        }
        
        // 5 - Classify training set
        applyRisk(patients, risk, &prediction[0]);

        // 6 - Get clasification error
        int cvTestError[Npatients];
        for (int j = 0; j < Npatients; ++j) {
            cvTestError[j] = (prediction[j] + npcases[0]) % 2; //TODO: npcases[j] instead npcases[0]
            cvTestError[j] = (1 - cvTestError[j]) * testset[i][j];
        }

        if (nTest> 0){
            int total=0;
            for(int k=0;k<Npatients;++k){
                total+=cvTestError[k];
            }
            testError[i]=(double)(total)/(double)(nTest);
        }          
        else{
            testError[i]=1; //#TODO: using 0 or 1
        }
    }
}

// Apply MDR to every SNP-SNP combination reading from a dict
void applyMdrDict (int row1 [Npatients], int row2[Npatients], float* testError,
                   int trainset [CV_sets] [Npatients], int testset [CV_sets] [Npatients],
                   int* npcases,int* npcontrols,double ccratio){

    int patients[Npatients];
    transformPatients(row1,row2,&patients[0]); 
    getRiskArray(&patients[0],testError,trainset,testset, npcases, npcontrols, ccratio);
}

// Save mdrerror to output directory ASCII format
void save_output(char keys1[50], char keys2[50],float* testError,char* outfile) {

    gzFile file = gzopen(outfile, "ab");
    gzputs(file, keys1);
    gzputs(file," ");
    gzputs(file, keys2);
    for (int k=0; k<CV_sets;++k){
        char doubleString[20];
        snprintf(doubleString, sizeof(doubleString), " %f", testError[k]); 
        gzputs(file, doubleString);
    }
    gzputs(file,"\n");
    gzclose(file);
    
    /*printf("%s %s ",keys1,keys2);
    for (int k=0; k<CV_sets;++k){
        printf(" %f",testError[k]);
    }
    printf("\n");*/
}

// Save mdrerror to output directory Bynary format
void save_output_bin(int i, int j,float* testError,char* outfile) {

    FILE *file = fopen(outfile, "ab");
    fwrite(&i, sizeof(int), 1, file);
    fwrite(&j, sizeof(int), 1, file);
    for (int k=0; k<CV_sets;++k){
        fwrite(&testError[k], sizeof(float), 1, file);
    }
    fclose(file);
}
