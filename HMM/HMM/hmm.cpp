#include "hmm.h"
/************************************************************************/
/* implementation                                                       */
/************************************************************************/

/*
	*phmm:已知的HMM模型
	T:序列长度
	*O:观察序列
	**alpha 局部概率
	*pprob 最终的观察概率
*/
void Forward(HMM *phmm, int T, int *O, double **alpha, double *pprob)
{
	int i,j;  //状态索引
	int t;    //time index
	double sum;  //temp value when calc local prob

	//1.初始化，计算t=1时刻的所有局部概率
	for (int i=1;i<=phmm->N;i++)
	{
		alpha[1][i]=phmm->pi[i]*phmm->B[i][O[1]];
	}
	//2.归纳，递归计算t>1时的所有局部概率
	for (t=1;t<T;t++)
	{
		for (j=1;j<=phmm->N;j++)
		{
			sum=0.0;
			for (i=1;i<=phmm->N;i++)
			{
				sum+=alpha[t][i]*(phmm->A[i][j]);
			}
			alpha[t+1][j]=sum*(phmm->B[j][O[t+1]]);
		}
	}

	//3.终止，观察序列的概率等于T时刻所有局域概率之和
	*pprob=0.0;
	for (i=1;i<=phmm->N;i++)
	{
		*pprob+=alpha[T][i];
	}
}
/*
	维特比算法
	*phmm: 已知的HMM模型
	T: 观察序列长度
	*O：观察序列
	**delta: 局部概率
	**psi：记录每一步的最优状态索引
	*q:最优路径
	*pprob：最优解
*/
void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob)
{
	int i,j;   /*状态索引*/
	int t;     //时间索引

	int maxvalind;
	double maxval, val;

	/*1. Initialization*/
	for (i=1;i<=phmm->N;i++)
	{
		delta[1][i]=phmm->pi[i]*(phmm->B[i][O[1]]);
		psi[1][i]=0;
	}

	/*2. Recursion，递归的求解每一个单元的值*/
	for (t=2;t<=T;t++)
	{
		for (j=1;j<=phmm->N;j++)
		{
			maxval=0.0;
			maxvalind=1;
			for (i=1;i<=phmm->N;i++)
			{
				val=delta[t-1][i]*(phmm->A[i][j]);
				if (val>maxval)
				{
					maxval=val;
					maxvalind=i;
				}
			}
			delta[t][j]=maxval*(phmm->B[j][O[t]]);
			psi[t][j]=maxvalind;	
		}
	}
	/*3.Termination 从最后一步的状态中选择一个概率最大的*/
	*pprob=0.0;
	q[T]=1;
	for (i=1;i<=phmm->M;i++)
	{
		if (delta[T][i]>*pprob)
		{
			*pprob=delta[T][i];
			q[T]=i;
		}
	}

	/*4. Path backtracking,回溯求解*/
	for (t=T-1;t>0;t--)
	{
		q[t]=psi[t+1][q[t+1]];
	}
}