#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "swipe.h"

#define MI 100
#define MA 500
#define ST .0001
#define DT .001

static int
one (struct dirent *unused)
{
  return 1;
}

int main(int argc, char* argv[]) {

	struct dirent ** instrumentos;
	struct dirent ** instrumento;
	int n_inst;
	int n_muest;
	int i,j;
	DIR * sb;
	DIR * rl;
	FILE * time_test;
	time_test = fopen ("time_test.txt", "w");
	
	fprintf(time_test, "File\ttime 1\ttime 2\ttime 3\ttotal time\n");
	/*
	
	n_inst = scandir("muestras/sonatina", &instrumentos, (int(*)(const struct dirent*))one, alphasort);
	//struct dirent * direntp;
	if(n_inst >= 0){
		
		for(i = 2; i < n_inst; ++i){
			printf("%s\n", instrumentos[i]->d_name);
			char inst[100];
			char inst_tmp[100];			
			strcpy(inst_tmp,"muestras/sonatina/");			
			strcat(inst_tmp, instrumentos[i]->d_name);
			n_muest = scandir(inst_tmp, &instrumento, (int(*)(const struct dirent*))one, alphasort);
			
			strcat(inst_tmp, "/");
			strcpy(inst,inst_tmp);
			for(j = 2; j < n_muest; ++j){
				char testDir[100];
				strcpy(testDir,"results/");
			
				printf("\t%s\n", instrumento[j]->d_name);
				strcpy(inst_tmp, inst);
				strcat(inst_tmp, instrumento[j]->d_name);
				
				vector p = swipe(inst_tmp,MI , MA, ST, DT, instrumento[j]->d_name, time_test);
				strcat(testDir, instrumento[j]->d_name);
				strcat(testDir,".txt");
				out(p.v, p.x, testDir);
				freev(p);
			
			} 
			
		} 
		
	}*/
	
	//instrumentos = opendir("../muestras/instrumentos");
	n_inst = scandir("muestras/instrumentos", &instrumentos, (int(*)(const struct dirent*))one, alphasort);
	//struct dirent * direntp;
	if(n_inst >= 0){
		
		for(i = 2; i < n_inst; ++i){
			printf("%s\n", instrumentos[i]->d_name);
			char inst[100];
			char inst_tmp[100];
			
			strcpy(inst_tmp,"muestras/instrumentos/");
			
			strcat(inst_tmp, instrumentos[i]->d_name);
			n_muest = scandir(inst_tmp, &instrumento, (int(*)(const struct dirent*))one, alphasort);
			
			strcat(inst_tmp, "/");
			strcpy(inst,inst_tmp);
			for(j = 2; j < n_muest; ++j){
				char testDir[100];
				strcpy(testDir,"results/");
				printf("\t%s\n", instrumento[j]->d_name);
				strcpy(inst_tmp, inst);
				strcat(inst_tmp, instrumento[j]->d_name);
				vector p = swipe(inst_tmp,MI , MA, ST, DT, instrumento[j]->d_name, time_test);
				strcat(testDir, instrumento[j]->d_name);
				strcat(testDir,".txt");
				out(p.v, p.x, testDir);
				freev(p);
			
			} 
			
		} 
		
	}
	
	
	/*
	struct dirent ** muestras;
	n_muest = scandir("../muestras/rl", &muestras, (int(*)(const struct dirent*))one, alphasort);
	for(j = 2; j < n_muest; ++j){
		char m[100];
		printf("\t%s\n", muestras[j]->d_name);
		strcpy(m, "../muestras/rl/");
		strcat(m, muestras[j]->d_name);
		vector p = swipe(m,MI , MA, ST, DT, muestras[j]->d_name, time_test);
		freev(p);
			
	}
	
	
	n_muest = scandir("../muestras/sb", &muestras, (int(*)(const struct dirent*))one, alphasort);
	for(j = 2; j < n_muest; ++j){
		char m[100];
		printf("\t%s\n", muestras[j]->d_name);
		strcpy(m, "../muestras/sb/");
		strcat(m, muestras[j]->d_name);
		vector p = swipe(m,MI , MA, ST, DT, muestras[j]->d_name, time_test);
		freev(p);
			
	}
	
	
	n_muest = scandir("../muestras/keele", &muestras, (int(*)(const struct dirent*))one, alphasort);
	for(j = 2; j < n_muest; ++j){
		char m[100];
		printf("\t%s\n", muestras[j]->d_name);
		strcpy(m, "../muestras/keele/");
		strcat(m, muestras[j]->d_name);
		vector p = swipe(m,MI , MA, ST, DT, muestras[j]->d_name, time_test);
		freev(p);
			
	}+*/
	fclose(time_test);
}

