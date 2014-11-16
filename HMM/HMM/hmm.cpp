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