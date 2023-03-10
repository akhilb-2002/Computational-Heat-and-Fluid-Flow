#include<bits/stdc++.h>

using namespace std;
typedef long long ll;

int main()
{
    long long int N;
	cout<<"Number of iterations: ";
	cin>>N;
	
	double phi1[N]={0};
	double phi2[N]={0};
	double phi3[N]={0};
	double phi4[N]={0};
	double phi5[N]={0};
	double error[N]={0};

	phi1[0]=100;phi2[0]=100;phi3[0]=100;phi4[0]=100;phi5[0]=100;

	for(int i=1;i<N;i++)
		{
			phi1[i]=(204 + phi2[i-1])/3;
			phi2[i]=(phi1[i] + phi3[i-1] + 12)/2;
			phi3[i]=(phi2[i] + phi4[i-1] + 20)/2;
			phi4[i]=(phi3[i] + phi5[i-1] + 28)/2;
			phi5[i]=(1236 + phi4[i])/3;

			error[i]=sqrt(pow((phi1[i]-phi1[i-1]),2) + pow((phi2[i]-phi2[i-1]),2) + pow((phi3[i]-phi3[i-1]),2) + pow((phi4[i]-phi4[i-1]),2) + pow((phi5[i]-phi5[i-1]),2));
		}

	for(int i=0;i<N;i++)
		{
			cout<<i<<"   "<<phi1[i]<<"   "<<phi2[i]<<"   "<<phi3[i]<<"   "<<phi4[i]<<"   "<<phi5[i]<<"  "<<error[i];
			cout<<endl;
		}
    return 0;
}


