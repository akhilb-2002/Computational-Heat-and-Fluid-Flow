#include<bits/stdc++.h>
#include<stdio.h>

using namespace std;
typedef long long ll;


//Solving a 2D grid of cells, with discrete equations of the form, ap*phi[p]=an*phi[n]+as*phi[s]+ae*phi[e]+aw*phi[w]+s, for each cell
//where p is cell P, n is adjacent cell radially outward, s is adjacent cell radially inward,
// e is adjacent cell in the angular clockwise direction, w is adjacent cell in the angular anti-clockwise direction
//Imagine ap,an,as,ae,aw,s to be arrays of size (N+1)*(M+1), in the order given in the main function

void PipewithinnerandouterradiussolvedbyangularTDMAsweepRadial (double Rout,double Rin, long long divtheta, long long int divrad)
{
    //divtheta is the number of divisions in the angular direction
    //divrad is the number of divisions in the radial direction

    double dtheta=2*3.14/double(divtheta);
    double dr=(Rout-Rin)/double(divrad);

    double gamma=0.001; double Sc=0.01;

    double ap[divrad+1][divtheta+1]={0};
    double an[divrad+1][divtheta+1]={0};
    double as[divrad+1][divtheta+1]={0};
    double ae[divrad+1][divtheta+1]={0};
    double aw[divrad+1][divtheta+1]={0};
    double s[divrad+1][divtheta+1]={0};

    //Boundary conditions
    for(int i=1;i<=divrad;i++)
    {
        for(int j=1;j<=divtheta;j++)
        {
            //Mid-region of the pipe
            if(i!=divrad && i!=1)
           {
            double Rp=Rin+(i-0.5)*dr;
            double dV=Rp*dr*dtheta;

            ae[i][j]=(gamma*dr)/(Rp*dtheta);
            aw[i][j]=(gamma*dr)/(Rp*dtheta);
            an[i][j]=(gamma*(Rin + i*dr)*dtheta)/dr;
            as[i][j]=(gamma*(Rin+(i-1)*dr)*dtheta)/dr;

            ap[i][j]=an[i][j]+as[i][j]+ae[i][j]+aw[i][j];
            s[i][j]=Sc*dV;
           }

           //Inner boundary of the pipe
           else if(i==1)
           {
            double Rp=Rin+(i-0.5)*dr;
            double dV=Rp*dr*dtheta;

            ae[i][j]=(gamma*dr)/(Rp*dtheta);
            aw[i][j]=(gamma*dr)/(Rp*dtheta);
            an[i][j]=(gamma*(Rin + i*dr)*dtheta)/dr;
            as[i][j]=0;

            ap[i][j]=ae[i][j]+aw[i][j]+ an[i][j] + 2*(gamma*(Rin + (i-1)*dr)*dtheta)/dr;
            s[i][j]=Sc*dV;
           }

           //Outer boundary of the pipe
           else if(i==divrad)
           {
            double Rp=Rin+(i-0.5)*dr;
            double dV=Rp*dr*dtheta;

            ae[i][j]=(gamma*dr)/(Rp*dtheta);
            aw[i][j]=(gamma*dr)/(Rp*dtheta);
            an[i][j]=0; 
            as[i][j]=(gamma*(Rin+(i-1)*dr)*dtheta)/dr;

            ap[i][j]=ae[i][j]+aw[i][j]+2*(gamma*(Rin+i*dr)*dtheta)/dr+as[i][j];
            s[i][j]=Sc*dV;
           }
        }
    }

   int noofsweeps= 0; //no of sweeps to be performed

   double phi[divrad+2][divtheta+2]={0};

   double error=0;
   
   do{
    double phiold[divrad+2][divtheta+2]={0};

    for(int i=1;i<=divrad;i++)
    {
        for(int j=1;j<=divtheta;j++)
        {
            phiold[i][j]=phi[i][j];
        }
    }

    error=0;

    // Since matrix is non-tridiagonal here, we use the Sherman-Morrison formula to solve the matrix equation
    for(int i=1;i<=divrad;i++)
    {
        double D[divtheta+1]={0};
        double C[divtheta+1]={0};
        double B[divtheta+1]={0};
        double A[divtheta+1]={0};

        for (int j=1;j<=divtheta;j++)
        {
            D[j]= an[i][j]*phi[i+1][j]+as[i][j]*phi[i-1][j]+s[i][j];
        }

        for(int k=1;k<=divtheta;k++)
        {
            C[k]=aw[i][k];
            B[k]=ae[i][k];
            A[k]=ap[i][k];
        }

        double u[divtheta+1]={0};
        double v[divtheta+1]={0};
        double y[divtheta+1]={0};
        double q[divtheta+1]={0};

        u[1]=gamma;
	    u[divtheta]=-B[divtheta];
	
	    v[1]=1;
	    v[divtheta]=-C[1]/gamma;
	
	    C[1]=0;
	    A[1]=A[1]-gamma;
	    A[divtheta]=A[divtheta]-(B[divtheta]*C[divtheta])/gamma;
	    B[divtheta]=0;
	
	    for(int k=2;k<=divtheta;k++){
		    double w=-C[k]/A[k-1];	
		    A[k]=A[k]+w*B[k-1];
		    D[k]=D[k]-w*D[k-1];
		    u[k]=u[k]-w*u[k-1];
	    }
	    y[divtheta]=D[divtheta]/A[divtheta];
	    q[divtheta]=u[divtheta]/A[divtheta];
	
	    for(int k=divtheta-1;k>=1;k--){
		    y[k]=(D[k]+B[k]*y[k+1])/A[k];
		    q[k]=(u[k]+B[k]*q[k+1])/A[k];
	    }
	    double Value=(v[1]*y[1]+v[divtheta]*y[divtheta])/(1+v[divtheta]*q[divtheta]+v[divtheta]*q[divtheta]);

        // We literally solve, phi=y- (q*transpose(v)*y/(1+transpose(v)*q))

	    for(int k=1;k<=divtheta;k++){
		    phi[i][k]=y[k]-q[k]*Value;
        }
    }

    for(int i=1;i<=divrad;i++)
    {
        for(int j=1;j<=divtheta;j++)
        {
            error=error+abs(phi[i][j]-phiold[i][j]);
        }
    }

    noofsweeps++;
   }while(error>1e-9);

    cout<<"After "<<noofsweeps<<" sweeps, the solution is:"<<endl;

    for(int i=1;i<=divrad;i++)
    {
        for(int j=1;j<=divtheta;j++)
        {
            cout<<phi[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"The notation is as follows:"<<endl;
    cout<<"phi[i][j] is the value of phi at the point (i,j)"<<endl;
    cout<<"i is the radial coordinate and j is the angular coordinate"<<endl;

}

int main()
{
    double Rout,Rin;
    long long int divtheta,divrad;

    cout<<"Enter the inner radius of the pipe: ";
    cin>>Rin;
    cout<<"Enter the outer radius of the pipe: ";
    cin>>Rout;
    cout<<"Enter the number of divisions in the radial direction: ";
    cin>>divrad;
    cout<<"Enter the number of divisions in the angular direction: ";
    cin>>divtheta;


    PipewithinnerandouterradiussolvedbyangularTDMAsweepRadial(Rin,Rout,divrad,divtheta);
    return 0;
}


