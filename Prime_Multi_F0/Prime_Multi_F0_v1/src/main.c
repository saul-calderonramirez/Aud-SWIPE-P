#include "AudSWIPEP.h"
/*
*Main method, interfacing with user arguments
*@param argc, number of arguments
*@param argv
*@return program result
*/



void provaCalcS(){
	printf("halloo, probando probando!!\n");

	vector pc = zerov(2);
	vector pcPrimes = zerov(2);
	matrix S = zerom(2, 3);
	pc.v[0] = 50; pc.v[1] = 200;
	//S = [ 0 0.0003495 0.002; 0.0018 0.0338 0.0207];

	S.m[0][0] = 0; S.m[0][1] = 0.0003495; S.m[0][2] = 0.002;
	S.m[1][0] = 0.0018; S.m[1][1] = 0.0338; S.m[1][2] = 0.0207;
	int numPrimes = (int)(pc.v[pc.x - 1] / pc.v[0]);
	intvector numsPrimesReason = primes(numPrimes);

	int i, j;

	for(i = 0; i < numsPrimesReason.x; ++i){
		printf("i: %d\n", numsPrimesReason.v[i]);
		int primeNum = numsPrimesReason.v[i];
		for(j = 0; j < pcPrimes.x; ++j){
			pcPrimes.v[j] = pc.v[j] * primeNum;
		}
		printf("Curr Vector primes:\n");
		printv(pcPrimes);
		matrix Sn = interp1Mat(pc, pcPrimes, S);
	}

}

int main(int argc, char* argv[]) {
		//executeAudSWIPEP(argc, argv);
		vector xO = zerov(2);
		xO.v[0] = 50; xO.v[1] = 200;
		vector yO = zerov(2);
		yO.v[0] = 0; yO.v[1] = 0.0018;
		vector xN = zerov(2);
		xN.v[0] = 100; xN.v[1] = 400;
		vector yN = zerov(2);
		interp1(xO, yO, xN, yN.v);
		printf("yN:\n");
		printv(yN);
		exit(EXIT_SUCCESS);

}
