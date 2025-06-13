#include <stdio.h>
#include <stdlib.h>
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include <string.h>
#include <math.h>
int natoms,step,natom,read_return,nFrame;
float time,prec;
matrix box;
rvec *x;
XDRFILE *xtc;
int a;
float protein_mass[2][3]={0},distAB;


extern float Distance(int n1, int n2) {
	float Dis,Dis2,X,Y,Z;
	X=fabsf(x[n1][0]-x[n2][0]);
	if(X>8) X=16-X;
	Y=fabsf(x[n1][1]-x[n2][1]);
	if(Y>8) Y=16-Y;
	Z=fabsf(x[n1][2]-x[n2][2]);
	if(Z>8) Z=16-Z;

	Dis2=X*X+Y*Y+Z*Z;
	Dis=sqrt(Dis2);
	return Dis;
}
extern float DistanceC(int n1, int n2)//计算 原子到质心的距离 
{
	float Dis,Dis2,X,Y,Z;
	X=fabsf(protein_mass[n1][0]-x[n2][0]);
	if(X>8) X=16-X;
	Y=fabsf(protein_mass[n1][1]-x[n2][1]);
	if(Y>8) Y=16-Y;
	Z=fabsf(protein_mass[n1][2]-x[n2][2]);
	if(Z>8) Z=16-Z;
	
	Dis2=X*X+Y*Y+Z*Z;
	Dis=sqrt(Dis2);
	return Dis;
 }
extern float DistanceOfPointTOLine(int i)
{ 	
	float ab,as,bs,cos_A,sin_A;
	 ab=distAB;// ab=sqrt(pow((a->x - b->x),2.0)+pow((a->y - b->y),2.0)+pow((a->z - b->z),2.0))
	 as=DistanceC(0,i);
	 bs=DistanceC(1,i);
	 cos_A=(pow(as,2.0)+pow(ab,2.0)-pow(bs,2.0))/(2*ab*as);
	 sin_A=sqrt(1-pow(cos_A,2.0));
	return as*sin_A;
}
extern void get_cID( int step )
{
	int line=0;
	char str[101];
	float X,Y,Z;
	FILE *f1 ;  
	if ( (f1 = fopen("mass.txt", "rt")) == NULL )  
	{
        puts("Fail to open file!");
        exit(0);
    }
    while( fgets (str, 101, f1) != NULL ) 
	{
		if(line == step) 
		{
			sscanf(str,"%d%f%f%f%f%f%f",&line,&protein_mass[0][0],&protein_mass[0][1],&protein_mass[0][2],&protein_mass[1][0],&protein_mass[1][1],&protein_mass[1][2]);
			break;
		}
		line++;
	}
	X=fabsf(protein_mass[0][0]-protein_mass[1][0]);
	if(X>8) X=16-X;
	Y=fabsf(protein_mass[0][1]-protein_mass[1][1]);
	if(Y>8) Y=16-Y;
	Z=fabsf(protein_mass[0][2]-protein_mass[1][2]);
	if(Z>8) Z=16-Z;
	distAB=sqrt(X*X+Y*Y+Z*Z);
	fclose(f1);
} 

void num_of_water() { 
	int n,s,Count=0;
	FILE *fWc;
	fWc = fopen("number_of_interfacial_water.dat","a+");
	for(n=4373; n<=539389; n+=4) 
	{
		for (int i=1;i<=2186;i++)
		{
			if(Distance (n-1,i-1)<=2)
			{
				for (int j=2187;j<=4372;j++)
				{
			 		if(Distance (n-1,j-1)<=2)
					{
						if(DistanceOfPointTOLine(n-1)<=2.6)
						{ 
						Count=Count+1;
						goto here;
						}
					}
			 	}
			goto here;	
			}
		}
		here:s=1;
	}
	fprintf(fWc,"%f %d",time,Count);
	fprintf(fWc,"\n");
	fclose(fWc);
}




void main() {
	xtc=xdrfile_open ("md.xtc","r");
	read_xtc_natoms ("md.xtc",&natoms);
	a=sizeof(x[0]);
	x = calloc(natoms, sizeof(x[0]));
	FILE *fz;
	nFrame=0;
	while (1) {
		read_return=read_xtc (xtc,natoms,&step,&time,box,x,&prec);
		if (read_return!=0) break;
		nFrame++;
		get_cID(time);
		num_of_water();
		printf("%d\n",nFrame);
	}
}
