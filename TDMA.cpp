#include<bits/stdc++.h>

using namespace std;
typedef long long ll;

// Function to solve the tridiagonal matrix equation
// Input in the form ai* phi(i)= bi* phi(i+1) + ci* phi(i-1) + di
//where i is cell P, i+1 is cell E, i-1 is cell W

void TDMA(double a[], double b[], double c[], double d[], long long N) // size of each array a,b,c,d is N+1, to avoid confusion,a[0]=b[0]..=0
{
    double P[N+1]={0};
    double Q[N+1]={0};

    P[1]=b[1]/a[1]; 
    Q[1]=d[1]/a[1];

    for(int i=2;i<N+1;i++)
    {
        P[i]=b[i]/(a[i]-c[i]*P[i-1]);
        Q[i]=(c[i]*Q[i-1]+d[i])/(a[i]-c[i]*P[i-1]);
    }

    double phi[N+1]={0};
    phi[N]=Q[N];

    for(int i=N;i>=2;i--)
    {
        phi[i-1]=P[i-1]*phi[i]+Q[i-1];
    }

    for(int i=1;i<=N;i++)
    {
        cout<<"phi["<<i<<"] = "<<phi[i]<<endl;
    } //phi[0]=0, is a dummy value, i=1 to N have corresponding phi values
}

int main()
{
    long long int N;
    cout<<"Number of cells :";
    cin>>N;

    double a[N+1]={0};
    double b[N+1]={0};
    double c[N+1]={0};
    double d[N+1]={0};

    for(int i=1;i<=N;i++)
    {
        cout<<"Enter a["<<i<<"] of cell "<<i<<": ";
        cin>>a[i];
        cout<<"Enter b["<<i<<"] of cell "<<i<<": ";
        cin>>b[i];
        cout<<"Enter c["<<i<<"] of cell "<<i<<": ";
        cin>>c[i];
        cout<<"Enter d["<<i<<"] of cell "<<i<<": ";
        cin>>d[i];
    }

    double phi[N+1]={0};

    TDMA(a,b,c,d,N);

    return 0;
}


