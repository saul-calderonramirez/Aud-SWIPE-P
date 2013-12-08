/*
 * Name: AudSWIPE-P
 * Authors: Prof. Arturo Camacho
 * 			Bsc. Saul Calderon
 * 			Bsc. Gabriel Alvarado
 * General Description: Parallel C implementation of the AudSWIPE algorithm, originally
 * implemented in MATLAB
 * Archive Description: Testing scripts generator
 */

#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include "AudSWIPEP.h"
#define MI 100
#define MA 500
//#define ST .0001
#define DT .01
#define PARALELISM 1
#define ARGC 12 //plus /0


#define TEST_DIR "../../Respaldos_Muestras/muestras/testDefinitivo/"

static int
one (struct dirent *unused)
{
  return 1;
}

/*
 * Writes a shell command in a file, MPI tests
 *
 */
void writeShCommand(FILE* shTest, char* route, char* archive, double paralelism, char* testFile, char* ptrResultsDir, char* FREQ_RANGE){
	fprintf(shTest, "echo Processing archive: %s \n", archive);
	fprintf(shTest, "mpiexec aswipep -i %s  -o %s/%s_Result.txt -r %s -p %f  -z %s\n", route, ptrResultsDir, archive, FREQ_RANGE, paralelism, testFile);
}

/*
 * Writes a shell command in a file, No MPI tests
 *
 */
void writeShCommandNonMPI(FILE* shTest, char* route, char* archive, char* ptrResultsDir, char* FREQ_RANGE, char* testFile ){
	fprintf(shTest, "echo Processing archive: %s \n", archive);
	if(testFile != NULL){
		fprintf(shTest, "aswipep -i %s  -o %s/%s_Result.txt -r %s  -z %s\n", route, ptrResultsDir, archive, FREQ_RANGE, testFile);
	}
	else{
		fprintf(shTest, "aswipep -i %s  -o %s/%s_Result.txt -r %s\n", route, ptrResultsDir, archive, FREQ_RANGE);
	}
}

/*
 * Runs tests in every subrdirectory, looking for wav files
 * @param mainDirectory, main directory name
 * @param test, test file to store results
 * @param ptrParalism, degree of paralelism.
 * */
void runTestInSubdirectories(char* mainDirectory, char* test, double paralelism, FILE* shTest, int mpi, char *ptrResults, char* ptrFreqRange){
	struct dirent ** instrumentos;
	struct dirent ** instrumento;
	int n_inst;
	int n_muest;
	int i,j;
	n_inst = scandir(mainDirectory, &instrumentos, (int(*)(const struct dirent*))one, alphasort);
	if(n_inst >= 0){
		for(i = 2; i < n_inst; ++i){
			char inst[100];
			char inst_tmp[100];
			strcpy(inst_tmp, mainDirectory);
			strcat(inst_tmp, instrumentos[i]->d_name);
			n_muest = scandir(inst_tmp, &instrumento, (int(*)(const struct dirent*))one, alphasort);
			strcat(inst_tmp, "/");
			strcpy(inst,inst_tmp);
			for(j = 2; j < n_muest; ++j){
				printf("\n Processing file: %s\n", instrumento[j]->d_name);
				strcpy(inst_tmp, inst);
				strcat(inst_tmp, instrumento[j]->d_name);

				if(mpi == 1){
					writeShCommand(shTest, inst_tmp, instrumento[j]->d_name, paralelism, test, ptrResults, ptrFreqRange );
				}
				else{
					writeShCommandNonMPI(shTest, inst_tmp, instrumento[j]->d_name, ptrResults, ptrFreqRange, test );
				}

			}
		}
	}
}

/*
 * Runs tests in one main directory, looking for wav files
 * @param mainDirectory, main directory name
 * @param test, test file to store results
 * @param ptrParalism, degree of paralelism.
 * */
void runTestInDirectory(char* mainDirectory, char* test, double paralelism, FILE* shTest, int mpi, char *ptrResults, char* ptrFreqRange){
	int n_muest;

	int j;

	struct dirent ** muestras;
	n_muest = scandir( mainDirectory, &muestras, (int(*)(const struct dirent*))one, alphasort);
	for(j = 2; j < n_muest; ++j){
		char m[100];
		strcpy(m, mainDirectory);
		strcat(m, muestras[j]->d_name);
		if(mpi == 1){
			writeShCommand(shTest, m, muestras[j]->d_name, paralelism, test, ptrResults, ptrFreqRange );
		}
		else{
			writeShCommandNonMPI(shTest, m, muestras[j]->d_name, ptrResults, ptrFreqRange, test );
		}
	}
}

/*
 * Runs autocorrelation tests
 * */
void runACtests(){
	FILE* shTest = fopen ("runTestsGabriel.sh", "w");
	char* ptrMainDir = TEST_DIR;
	char* ptrResultsDir = "ResultsGabriel";
	fprintf(shTest, "cd ../OMP_Gabriel\n");
	fprintf(shTest, "sudo make clean\n");
	fprintf(shTest, "sudo make install\n");
	runTestInSubdirectories(ptrMainDir, "testGabriel_P100.csv", 1, shTest, 0, ptrResultsDir, "50:2000");
	fclose(shTest);
}


/*
 * Runs paralelism tests
 * */
void runParalelismTests(){
	FILE* shTest = fopen ("runParalelismTests.sh", "w");
	char* ptrMainDir = TEST_DIR;
	char* ptrResultsDir = "ResultsParalelism";
	fprintf(shTest, "sudo make clean\n");
	fprintf(shTest, "sudo make install\n");
	runTestInSubdirectories(ptrMainDir, "test_P0.csv", 0, shTest, 1, ptrResultsDir, "50:2000");
	runTestInSubdirectories(ptrMainDir,  "test_P50.csv", 0.5, shTest, 1, ptrResultsDir, "50:2000");
	runTestInSubdirectories(ptrMainDir, "test_P100.csv", 1, shTest, 1, ptrResultsDir, "50:2000");
	fclose(shTest);
}

/*
 * Generates paralelism tests
 * */
vector generateParalelismLevels(int fs, int min, int max){
	int maxvent = round(log2(8*fs/min));
	int	minvent = round(log2(8*fs/max));
	int numberWs = maxvent - minvent + 1;
	vector paralelismVector = zerov(numberWs);
	int i;
	for(i = 0; i < numberWs; ++i){
		paralelismVector.v[i] = ((double)(i+1)/numberWs) + 0.001 ;
		if(paralelismVector.v[i] > 1) paralelismVector.v[i] = 1;
	}
	return paralelismVector;
}

/*
 * Generates the scripts
 * */
void scriptInstrumentTest(FILE* shTest, char *ptrInstrument, int min, int max, int fs,  int mpi, int test){
		char ptrMainDir[200] = "../../Respaldos_Muestras/muestras/testDefinitivo/";
		char ptrResultsDir[200] = "Results";
		char ptrTestCSV[200] = "Test_";
		char *ptrPtrTest = NULL;
		char ptrFreqRange[30] = "";
		char ptrMin[10];
		char ptrMax[10];
		sprintf ( ptrMin, "%d", min );
		sprintf ( ptrMax, "%d", max );
		strcat(ptrFreqRange, ptrMin);
		strcat(ptrFreqRange, ":");
		strcat(ptrFreqRange, ptrMax);
		strcat(ptrMainDir, ptrInstrument);
		strcat(ptrMainDir, "/");
		strcat(ptrResultsDir, ptrInstrument);
		strcat(ptrTestCSV, ptrInstrument);
		strcat(ptrTestCSV, ".csv");
		fprintf(shTest, "mkdir %s\n", ptrResultsDir);
		//generates the paralelism levels to try (all the quantity of processes possible)
		if(test == 1)ptrPtrTest = ptrTestCSV;
		if(mpi == 1){
			vector paralelismVector = generateParalelismLevels(fs, min, max);
			int i;
			for(i = 0; i < paralelismVector.x; ++i){
				runTestInDirectory(ptrMainDir,  ptrPtrTest, paralelismVector.v[i], shTest, 1, ptrResultsDir, ptrFreqRange);
			}
		}
		else{
			runTestInDirectory(ptrMainDir,  ptrPtrTest, 0, shTest, 0, ptrResultsDir,  ptrFreqRange);
		}
}

/*
 * Generates the parallelism tests
 * */
void runInstrumentsTimeTestsParallel(){
	FILE* shTest = fopen ("runInstrumentsTimeTestsParallel.sh", "w");
	fprintf(shTest, "sudo make clean\n");
	fprintf(shTest, "sudo make install\n");
	fprintf(shTest, "mkdir TimeResults\n" );
	int mpi = 1;
	int timeFile = 1;
	scriptInstrumentTest( shTest, "Basses", 30, 500, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Tuba", 37, 350, 10000, mpi, timeFile);
	scriptInstrumentTest( shTest, "Violas", 131, 2100, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Violin", 196, 3520, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Trumpet", 150, 1175, 10000, mpi, timeFile);
	scriptInstrumentTest( shTest, "Trombones", 73, 700, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Piccolo", 523, 3951, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Flute", 250, 2349, 10000, mpi, timeFile);
	scriptInstrumentTest( shTest, "Oboe", 250, 1760, 10000, mpi, timeFile);
	scriptInstrumentTest( shTest, "Clarinet", 131, 1760, 44100, mpi, timeFile);
	printf("\ntest Script done\n");
	fclose(shTest);
}

/*
 * Generates the sequential tests
 * */
void runInstrumentsTimeTestsSequential(){
	FILE* shTest = fopen ("runInstrumentsTimeTestsNoMPI.sh", "w");
	fprintf(shTest, "sudo make clean\n");
	fprintf(shTest, "sudo make install\n");
	fprintf(shTest, "mkdir TimeResults\n" );
	int mpi = 0;
	int timeFile = 1;
	scriptInstrumentTest( shTest, "Basses", 30, 500, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Tuba", 37, 350, 10000, mpi, timeFile);
	scriptInstrumentTest( shTest, "Violas", 131, 2100, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Violin", 196, 3520, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Trumpet", 150, 1175, 10000, mpi, timeFile);
	scriptInstrumentTest( shTest, "Trombones", 73, 700, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Piccolo", 523, 3951, 44100, mpi, timeFile);
	scriptInstrumentTest( shTest, "Flute", 250, 2349, 10000, mpi, timeFile);
	scriptInstrumentTest( shTest, "Oboe", 250, 1760, 10000, mpi, timeFile);
	scriptInstrumentTest( shTest, "Clarinet", 131, 1760, 44100, mpi, timeFile);
	fclose(shTest);
}

/*
 * Main method
 * @param -p, the degree of paralelism 0 or 1
 * */
/*
int main(int argc, char* argv[]) {
	runInstrumentsTimeTestsSequential();
		exit(EXIT_SUCCESS);
}*/




