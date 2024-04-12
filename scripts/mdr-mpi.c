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

#include "mdr-mpi.h"

int main(int argc, char *argv[]) {
   
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the rank of the current process and the total number of processes
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check if there are at least 2 processes
    if (size < 2) {
        printf("This program requires at least 2 processes");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    
    // Default value for nfiles
    int nfiles = files2proc; 
    if (argc > 1) {
        nfiles = atoi(argv[1]); 
    }

    // Read patients information
	int labels [Npatients];
	int npcontrols [Npatients];
	int* npcases= &labels[0];

    read_labels(DATAPATH"/labels.sample",&labels[0]);
    int total_cases=convert_labels(npcases,&npcontrols[0]);

	// Get cases/controls ratio. 
    // We will use the number as the high risk/low risk separator
	int total_ctrls=Npatients-total_cases;
	double ccratio = (double)total_cases / (double)total_ctrls;

	// Create training and test set
    double block_size = (Npatients/5.0);
	
	int trainset [CV_sets] [Npatients];
	int testset [CV_sets] [Npatients];

	// Initialize all elements of the array to 1/0
    for (int i = 0; i < CV_sets; ++i) {
	    for (int j = 0; j < Npatients; j++){
		    trainset[i][j] = 1;
		    testset[i][j] = 0;
	    }
    }
    //Set training and test set
 	for (int i = 0; i < CV_sets; ++i) {
		for (int j = (i * block_size); j < (i * block_size + block_size); ++j) 
		{
			trainset[i][j] = 0;
			testset [i][j] = 1;
		}
	}

    // Get files to be processed in the node 
    // TODO: try other distribution int bloque=(nfiles+size-1)/size;
    int node_files = nfiles / size;
    int first = rank * node_files;
    int last = (rank + 1) * node_files;
    if (rank == size - 1){
        last = nfiles;
    }
    int combxnode=last-first;
    printf(" -- Total number of files to be processed %d.\n", nfiles);
    printf(" -- %d combination of files for the node %d:\n", combxnode,rank);
    printf(" --\n");

    // Read files to be processed - TODO send these files to each processor readed
    char files	[files2proc*MAX_FILENAME_LENGTH]; // TODO files2proc or combxnode
    files[0] = '\0';
    get_comb_files(LISTOFFILESPATH,&files[0],first,last);

    // Compute all combinations
    int total_pairs = 0;
    int i=0;
    while(files[i]!= '\0'){
        char f1[MAX_FILENAME_LENGTH];
        char f2[MAX_FILENAME_LENGTH];
        int j=0; int k=0;
        while (files[i]!='-'){
            f1[j++]=files[i++];
        }
        f1[j]='\0'; i++;
        while (files[i]!='\n'){
            f2[k++]=files[i++];
        }
        f2[k]='\0'; i++;
        printf("> Loading files %s and %s at node %d\n",f1,f2,rank);

        // Preparing an empty file to Save results to file
        char outfile[MAX_FILENAME_LENGTH]=OUTPUTPATH"/MDR_";
        strcat(outfile,f1);
        strcat(outfile,"-");
        strcat(outfile,f2);
        strcat(outfile,".bin");
        
        FILE *MDR_outfile = fopen(&outfile[0], "wb");
        fclose(MDR_outfile);

        // Read files and save the data in a key-value arrays
        char keys1[nsamples][50]; //50 char per chrm name
        char keys2[nsamples][50]; //50 char per chrm name
        int values1[nsamples][Npatients]; 
        int values2[nsamples][Npatients]; 
        
        char inputFile1[MAX_FILENAME_LENGTH]=SAMPLEPATH;
        char inputFile2[MAX_FILENAME_LENGTH]=SAMPLEPATH;
        strcat(inputFile1,f1);
        strcat(inputFile1,".gz");
        strcat(inputFile2,f2);
        strcat(inputFile2,".gz");
        read_sample(&inputFile1[0],keys1,values1);
        read_sample(&inputFile2[0],keys2,values2);
        
        // Get all the combinations
        float testError[CV_sets];
        for (int i=0;i<nsamples;i++){
            for (int j=0;j<nsamples;j++) //use j = i+1 to avoid mxm comparations
            {
                applyMdrDict(values1[i],values2[j],&testError[0],
                             trainset, testset,
                             npcases,&npcontrols[0],ccratio);
                //save_output(keys1[i],keys2[j],&testError[0],&outfile[0]); //Save as ASCII
                save_output_bin(i,j,&testError[0],&outfile[0]); 
                total_pairs++;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);   
    printf("total_pairs computed per node %d: %d\n",rank,total_pairs);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
