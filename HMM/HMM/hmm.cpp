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
/*
	ά�ر��㷨
	*phmm: ��֪��HMMģ��
	T: �۲����г���
	*O���۲�����
	**delta: �ֲ�����
	**psi����¼ÿһ��������״̬����
	*q:����·��
	*pprob�����Ž�
*/
void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob)
{
	int i,j;   /*״̬����*/
	int t;     //ʱ������

	int maxvalind;
	double maxval, val;

	/*1. Initialization*/
	for (i=1;i<=phmm->N;i++)
	{
		delta[1][i]=phmm->pi[i]*(phmm->B[i][O[1]]);
		psi[1][i]=0;
	}

	/*2. Recursion���ݹ�����ÿһ����Ԫ��ֵ*/
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
	/*3.Termination �����һ����״̬��ѡ��һ����������*/
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

	/*4. Path backtracking,�������*/
	for (t=T-1;t>0;t--)
	{
		q[t]=psi[t+1][q[t+1]];
	}
}