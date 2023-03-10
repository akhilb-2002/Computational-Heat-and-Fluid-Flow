#include<bits/stdc++.h>
#include<stdio.h>

using namespace std;
typedef long long ll;

// Function to solve the tridiagonal matrix equation
//Solving a 2D grid of cells, with discrete equations of the form, ap*phi[p]=an*phi[n]+as*phi[s]+ae*phi[e]+aw*phi[w]+s, for each cell
//where p is cell P, n is cell N, s is cell S, e is cell E, w is cell W
//Imagine ap,an,as,ae,aw,s to be arrays of size (N+1)*(M+1), in the order given in the main function

void TDMAySWEEPx (double ap[],double an[],double as[],double ae[],double aw[],double s[],long long N,long long M)
{
   double aP[N+1][M+1]={0};
   double aN[N+1][M+1]={0};
   double aS[N+1][M+1]={0};
   double aE[N+1][M+1]={0};
   double aW[N+1][M+1]={0};
   double S[N+1][M+1]={0};

   for(int i=1;i<=N;i++)
   {
    for(int j=1;j<=M;j++)
    {
        aP[i][j]=ap[(j-1)*N+i];
        aN[i][j]=an[(j-1)*N+i];
        aS[i][j]=as[(j-1)*N+i];
        aE[i][j]=ae[(j-1)*N+i];
        aW[i][j]=aw[(j-1)*N+i];
        S[i][j]=s[(j-1)*N+i];
    }
   }

   int noofsweeps= 0; //no of sweeps to be performed

   double phi[N+2][M+2]={0};

   double error=0;

   
   do{
    double phiold[N+2][M+2]={0};

    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=M;j++)
        {
            phiold[i][j]=phi[i][j];
        }
    }

    error=0;

    for(int j=1;j<=M;j++)
    {
        double d[N+1]={0};

        for (int i=1;i<=N;i++)
        {
            d[i]= aW[i][j]*phi[i][j-1]+aE[i][j]*phi[i][j+1]+S[i][j];
        }

        double P[N+1]={0};
        double Q[N+1]={0};

        P[1]=aN[1][j]/aP[1][j];
        Q[1]=d[1]/aP[1][j];

        for (int k=2;k<=N;k++)
        {
            P[k]=aN[k][j]/(aP[k][j]-aS[k][j]*P[k-1]);
            Q[k]=(d[k]+aS[k][j]*Q[k-1])/(aP[k][j]-aS[k][j]*P[k-1]);
        }

        phi[N][j]=Q[N];

        for(int k=N;k>=2;k--)
        {
            phi[k-1][j]=P[k-1]*phi[k][j]+Q[k-1];
        }
       
    }

    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=M;j++)
        {
            error=error+pow((phi[i][j]-phiold[i][j]),2);
        }
    }

    noofsweeps++;
   }while(error>0.001);

    cout<<"After "<<noofsweeps<<" sweeps, the solution is:"<<endl;

    for(int j=1;j<=M;j++)
    {
         for(int i=1;i<=N;i++)
         {
            cout<<"phi["<<i<<"]["<<j<<"]="<<phi[i][j]<<" ";
         }
         cout<<endl;
    }

    cout<<"**Notation : phi[i][j] is the value of phi at cell at ith row and jth column**"<<endl;
}

void TDMAxSWEEPy (double ap[],double an[],double as[],double ae[],double aw[],double s[],long long N,long long M)
{
   double aP[N+1][M+1]={0};
   double aN[N+1][M+1]={0};
   double aS[N+1][M+1]={0};
   double aE[N+1][M+1]={0};
   double aW[N+1][M+1]={0};
   double S[N+1][M+1]={0};

   for(int i=1;i<=N;i++)
   {
    for(int j=1;j<=M;j++)
    {
        aP[i][j]=ap[(j-1)*N+i];
        aN[i][j]=an[(j-1)*N+i];
        aS[i][j]=as[(j-1)*N+i];
        aE[i][j]=ae[(j-1)*N+i];
        aW[i][j]=aw[(j-1)*N+i];
        S[i][j]=s[(j-1)*N+i];
    }
   }

   int noofsweeps= 0; //no of sweeps to be performed

   double phi[N+2][M+2]={0};

   double error=0;

   do{
    double phiold[N+2][M+2]={0};

    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=M;j++)
        {
            phiold[i][j]=phi[i][j];
        }
    }

    error=0;

    for(int i=1;i<=N;i++)
    {
        double d[M+1]={0};

        for (int j=1;j<=M;j++)
        {
            d[j]= aN[i][j]*phi[i+1][j]+aS[i][j]*phi[i-1][j]+S[i][j];
        }

        double P[M+1]={0};
        double Q[M+1]={0};

        P[1]=aE[i][1]/aP[i][1];
        Q[1]=d[1]/aP[i][1];

        for (int k=2;k<=M;k++)
        {
            P[k]=aE[i][k]/(aP[i][k]-aW[i][k]*P[k-1]);
            Q[k]=(d[k]+aW[i][k]*Q[k-1])/(aP[i][k]-aW[i][k]*P[k-1]);
        }

        phi[i][M]=Q[M];

        for(int k=M;k>=2;k--)
        {
            phi[i][k-1]=P[k-1]*phi[i][k]+Q[k-1];
        }
       
    }

    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=M;j++)
        {
            error=error+pow((phi[i][j]-phiold[i][j]),2);
        }
    }

    noofsweeps++;
   }while(error>0.01);

    cout<<"After "<<noofsweeps<<" sweeps, the solution is:"<<endl;

    for(int i=1;i<=N;i++)
    {
         for(int j=1;j<=M;j++)
         {
            cout<<"phi["<<i<<"]["<<j<<"]="<<phi[i][j]<<" ";
         }
         cout<<endl;
    }

    cout<<"**Notation : phi[i][j] is the value of phi at cell at ith row and jth column**"<<endl;
}

int main()
{
    long long int N;
    long long int M;
    cout<<"Number of rows:";
    cin>>N;
    cout<<"Number of columns:";
    cin>>M;

    double ap[(N+1)*(M+1)]={0};
    double an[(N+1)*(M+1)]={0};
    double as[(N+1)*(M+1)]={0};
    double ae[(N+1)*(M+1)]={0};
    double aw[(N+1)*(M+1)]={0};
    double s[(N+1)*(M+1)]={0};

    for(int j=1;j<=M;j++)
    {
        for(int i=1;i<=N;i++)
        {
            cout<<"ap["<<i<<"]["<<j<<"]:";
            cin>>ap[(j-1)*N+i];
            cout<<"an["<<i<<"]["<<j<<"]:";
            cin>>an[(j-1)*N+i];
            cout<<"as["<<i<<"]["<<j<<"]:";
            cin>>as[(j-1)*N+i];
            cout<<"ae["<<i<<"]["<<j<<"]:";
            cin>>ae[(j-1)*N+i];
            cout<<"aw["<<i<<"]["<<j<<"]:";
            cin>>aw[(j-1)*N+i];
            cout<<"s["<<i<<"]["<<j<<"]:";
            cin>>s[(j-1)*N+i];
        }
    }
    //Test code
    /*N=4;
    M=3;
    double ap[]={0,20,30,30,40,30,40,40,50,20,30,30,40,0,0,0,0,0,0,0};
    double an[]={0,10,10,10,0,10,10,10,0,10,10,10,0,0,0,0,0,0,0,0};
    double as[]={0,0,10,10,10,0,10,10,10,0,10,10,10,0,0,0,0,0,0,0};
    double aw[]={0,0,0,0,0,10,10,10,10,10,10,10,10,0,0,0,0,0,0,0};
    double ae[]={0,10,10,10,10,10,10,10,10,0,0,0,0,0,0,0,0,0,0,0};
    double s[]={0,500,500,500,2500,0,0,0,2000,0,0,0,2000,0,0,0,0,0,0,0};*/


    TDMAxSWEEPy(ap,an,as,ae,aw,s,N,M);
    TDMAySWEEPx(ap,an,as,ae,aw,s,N,M);

    return 0;
}


