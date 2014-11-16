#include "hmm.h"
/************************************************************************/
/* implementation                                                       */
/************************************************************************/

/*
	*phmm:��֪��HMMģ��
	T:���г���
	*O:�۲�����
	**alpha �ֲ�����
	*pprob ���յĹ۲����
*/
void Forward(HMM *phmm, int T, int *O, double **alpha, double *pprob)
{
	int i,j;  //״̬����
	int t;    //time index
	double sum;  //temp value when calc local prob

	//1.��ʼ��������t=1ʱ�̵����оֲ�����
	for (int i=1;i<=phmm->N;i++)
	{
		alpha[1][i]=phmm->pi[i]*phmm->B[i][O[1]];
	}
	//2.���ɣ��ݹ����t>1ʱ�����оֲ�����
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

	//3.��ֹ���۲����еĸ��ʵ���Tʱ�����о������֮��
	*pprob=0.0;
	for (i=1;i<=phmm->N;i++)
	{
		*pprob+=alpha[T][i];
	}
}