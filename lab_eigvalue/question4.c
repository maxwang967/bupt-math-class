#include <stdio.h>
#include <math.h>

#define ERR	1.0e-12	//误差限
#define N	11		//矩阵行列数
#define L	1.0e5	//最大迭代次数

double A[N][N]={0};

void Init_A()		//初始化矩阵
{
	for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++){
                A[i][j] = 0;
            }
        }
        for (int i = 1; i < N; i ++){
            A[i][i-1] = 1;
        }
        for (int i = 0; i < N; i ++){
            A[i][10] = -1;
        }
}
/*
void Display_A()
{
	int i,j;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
			printf("%8.4f",A[i][j]);
		printf("\n");
	}
}
*/
int Sgn(double a)
{
	if(a>=0)
		return 1;
	else
		return -1;
}

void On_To_The_Triangle()	//矩阵拟上三角化
{
	int i,j,r,flag=0;
	double cr,dr,hr,tr,temp;
	double ur[N],pr[N],qr[N],wr[N];
	for(r=1;r<=N-2;r++)
	{
		flag=0;
		for(i=r+2;i<=N;i++)
			if(A[i-1][r-1]!=0)
			{
				flag=1;
				break;
			}
		if(0==flag)
			continue;
		dr=0;
		for(i=r+1;i<=N;i++)
			dr+=A[i-1][r-1]*A[i-1][r-1];
		dr=sqrt(dr);
		if(0==A[r][r-1])
			cr=dr;
		else cr=-Sgn(A[r][r-1])*dr;
		hr=cr*cr-cr*A[r][r-1];
		for(i=1;i<=r;i++)
			ur[i-1]=0;
		ur[r]=A[r][r-1]-cr;
		for(i=r+2;i<=N;i++)
			ur[i-1]=A[i-1][r-1];
		for(i=1;i<=N;i++)
		{
			temp=0;
			for(j=1;j<=N;j++)
				temp+=A[j-1][i-1]*ur[j-1];
			pr[i-1]=temp/hr;
		}
		for(i=1;i<=N;i++)
		{
			temp=0;
			for(j=1;j<=N;j++)
				temp+=A[i-1][j-1]*ur[j-1];
			qr[i-1]=temp/hr;
		}
		temp=0;
		for(i=1;i<=N;i++)
		{
			temp+=pr[i-1]*ur[i-1];
			tr=temp/hr;
		}
		for(i=1;i<=N;i++)
		{
			wr[i-1]=qr[i-1]-tr*ur[i-1];
		}
		for(i=1;i<=N;i++)
			for(j=1;j<=N;j++)
				A[i-1][j-1]=A[i-1][j-1]-wr[i-1]*ur[j-1]-ur[i-1]*pr[j-1];
	}
}

void Get_Roots(double eigenvalue[][2],int m,double ss,double tt)	//求一元二次方程的根
{
	double discriminant=ss*ss-4*tt;	//
	if(discriminant<0)
	{
		*(*(eigenvalue+m-2))=0.5*ss;
		*(*(eigenvalue+m-2)+1)=0.5*sqrt(-discriminant);
		*(*(eigenvalue+m-1))=0.5*ss;
		*(*(eigenvalue+m-1)+1)=-0.5*sqrt(-discriminant);
	}
	else
	{
		*(*(eigenvalue+m-2))=0.5*(ss+sqrt(discriminant));
		*(*(eigenvalue+m-2)+1)=0;
		*(*(eigenvalue+m-1))=0.5*(ss-sqrt(discriminant));
		*(*(eigenvalue+m-1)+1)=0;
	}
}

void Get_Mk(double mk[][N],int m,double ss,double tt)	//获取Mk,用于带双步位移的QR分解
{
	int i,j,k;
	for(i=0;i<m;i++)
		for(j=0;j<m;j++)
			*(*(mk+i)+j)=0;
	for(i=0;i<m;i++)
		for(j=0;j<m;j++)
		{
			for(k=0;k<m;k++)
				*(*(mk+i)+j)+=A[i][k]*A[k][j];
			*(*(mk+i)+j)-=ss*A[i][j];
			if(j==i)
				*(*(mk+i)+j)+=tt;
		}
}

void QR_Reslove(double mk[][N],int m)	//QR分解
{
	int i,j,r,flag=0;
	double cr,dr,hr,tr,temp;
	double ur[N],vr[N],pr[N],qr[N],wr[N];
	double B[N][N],C[N][N];
	for(i=0;i<m;i++)
		for(j=0;j<m;j++)
		{
			B[i][j]=*(*(mk+i)+j);
			C[i][j]=A[i][j];
		}
	for(r=1;r<=m-1;r++)
	{
		flag=0;
		for(i=r+1;i<=m;i++)
			if(B[i-1][r-1]!=0)
			{
				flag=1;
				break;
			}
		if(0==flag)
			continue;
		dr=0;
		for(i=r;i<=m;i++)
			dr+=B[i-1][r-1]*B[i-1][r-1];
		dr=sqrt(dr);
		if(0==B[r-1][r-1])
			cr=dr;
		else cr=-Sgn(B[r-1][r-1])*dr;
		hr=cr*cr-cr*B[r-1][r-1];
		for(i=1;i<r;i++)
			ur[i-1]=0;
		ur[r-1]=B[r-1][r-1]-cr;
		for(i=r+1;i<=m;i++)
			ur[i-1]=B[i-1][r-1];
		for(i=1;i<=m;i++)
		{
			temp=0;
			for(j=1;j<=m;j++)
				temp+=B[j-1][i-1]*ur[j-1];
			vr[i-1]=temp/hr;
		}
		for(i=0;i<m;i++)
			for(j=0;j<m;j++)
				B[i][j]-=ur[i]*vr[j];
		for(i=1;i<=m;i++)
		{
			temp=0;
			for(j=1;j<=m;j++)
				temp+=C[j-1][i-1]*ur[j-1];
			pr[i-1]=temp/hr;
		}
		for(i=1;i<=m;i++)
		{
			temp=0;
			for(j=1;j<=m;j++)
				temp+=C[i-1][j-1]*ur[j-1];
			qr[i-1]=temp/hr;
		}
		temp=0;
		for(i=1;i<=m;i++)
		{
			temp+=pr[i-1]*ur[i-1];
			tr=temp/hr;
		}
		for(i=1;i<=m;i++)
		{
			wr[i-1]=qr[i-1]-tr*ur[i-1];
		}
		for(i=1;i<=m;i++)
			for(j=1;j<=m;j++)
				C[i-1][j-1]=C[i-1][j-1]-wr[i-1]*ur[j-1]-ur[i-1]*pr[j-1];
	}
	for(i=0;i<m;i++)
		for(j=0;j<m;j++)
			A[i][j]=C[i][j];
}

void Display_Eigenvalue(double value[][2])	//显示特征值
{
	int i;
	for(i=0;i<N;i++)
	{
		printf("λ%d=%8.4f",i+1,*(*(value+i)));
		if(*(*(value+i)+1)>0)
			printf("+%8.4f",*(*(value+i)+1));
		else if(*(*(value+i)+1)<0)
			printf("%8.4f",*(*(value+i)+1));
		printf("\n");
	}
	printf("\n");
}

int QR_With_Double_Step_Displacement(double eigenvalue[][2])	//带双步位移QR分解求特征值
{
	int i,j,k=1,m=N;
	double s,t;
	double Mk[N][N];

	for(i=0;i<N;i++)
		for(j=0;j<2;j++)
			eigenvalue[i][j]=0;

	do
	{
		k++;
		if(m==1)
		{
			eigenvalue[m-1][0]=A[m-1][m-1];
			m--;
			continue;
		}
		else if(m==2)
		{
			s=A[m-2][m-2]+A[m-1][m-1];
			t=A[m-2][m-2]*A[m-1][m-1]-A[m-1][m-2]*A[m-1][m-2];
			Get_Roots(eigenvalue,m,s,t);	//求一元二次方程的根
			m=0;
			continue;
		}
		else if(m==0)
			return 0;
		else if(fabs(A[m-1][m-2])<=ERR)
		{
			eigenvalue[m-1][0]=A[m-1][m-1];
			m--;
			continue;
		}
		else
		{
			s=A[m-2][m-2]+A[m-1][m-1];
			t=A[m-2][m-2]*A[m-1][m-1]-A[m-1][m-2]*A[m-2][m-1];
			if(fabs(A[m-2][m-3])<=ERR)
			{
				Get_Roots(eigenvalue,m,s,t);	//求一元二次方程的根
				m-=2;
				continue;
			}
			else
			{
				Get_Mk(Mk,m,s,t);		//获取Mk,用于带双步位移的QR分解
				QR_Reslove(Mk,m);		//QR分解
			}
		}
	}
	while(k<L);
	printf("ERROR\n");					//迭代过程不成功，报错
}

main()
{
	double eigenvalue[N][2];

	Init_A();									//初始化矩阵

	On_To_The_Triangle();						//矩阵上三角化

	QR_With_Double_Step_Displacement(eigenvalue);	//带双步位移QR分解求特征值

	Display_Eigenvalue(eigenvalue);				//显示特征值

	return 0;
}
