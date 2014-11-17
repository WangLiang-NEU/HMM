#ifndef _HMM_
#define _HMM_
/************************************************************************/
/* data structure                                                       */
/* @Author Lyon-NEU                                                     */
/* @date  2014.11.16                                                    */
/************************************************************************/
#include <stdio.h>
typedef struct{
	int N;       //number of states 隐藏状态数目
	int M;       //number of observation symbols 观察序列状态数
	double **A;  //A[1..N][1..N]. a[i][j] is the transition probability of going from state i at time t to state j at time t+1
	double **B;  //B[1..N][1..M]. b[j][k] is the probability of observing symbol k in state j
	double *pi;  //pi[1..N] pi[i] is the initial state distribution
}HMM;

void readHMM(FILE* fp,HMM* hmm);
void printHMM(FILE* fp,HMM* hmm);
void copyHMM(HMM* h1,HMM*h2);
void freeHMM(HMM* hmm);

void Forward(HMM *phmm, int T, int *O, double **alpha, double *pprob);
void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob);
#endif