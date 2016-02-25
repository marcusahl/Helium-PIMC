//
//  main.cpp
//  optimizedharm
//
//  Created by Metalmarcus on 18/02/14.
//  Copyright (c) 2014 Metalmarcus. All rights reserved.
//


#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "mtrand.h"


using namespace std;

const int jmax = 6;

int M=0; //number of steps
double dt; //imaginary time step
double B; //Inverse temperature
double temperature[jmax];
const int l=5; //Number of bead levels
const int length = pow(2,l)+1;

const int tau = 10000; //Correlation time
const int tfinal = 6*tau; //Total run time
int tsamp;
const int samp= tfinal-tau; //sampling time
const int swaptime = 100; //swaps per timestep
const int particles=2; //# of particles
const int dim=2; //dimensions


double L=0; //box side length or disc diameter
const int space=100000; // possible coordinates in the box
int ss=space-1;


double r [space]; //grid of space points

double Ttemp;
double Vtemp;
double R_mean_trial [dim];
double R_diff;
double R_trial_diff;

double h[space];
double d=1.241314;

double ri [2]; //attempted new particle spot
double rj [2];

int permutation [particles]; //Keeps track of which line it enters in the PBC
int inversePermutation [particles];

double A; //acceptance ratio
double S; //Second derivative of action
double E; //Energy
double D; //Distance
double RMS; //RootMeanSquare distance

double RMSsamp;
double Dsamp;
double Esamp;
double E2samp;
double Ssamp;
double Tsamp;
double Vsamp;
double Dav[jmax];
double Sav[jmax];
double Eav[jmax];
double E2av[jmax];
double Cv[jmax];
double Tav[jmax];
double Vav[jmax];
double RMSav [jmax];
double ACCEPT [jmax];
double SWAP [jmax];

int acceptance=0;
int swapCount =0;
double dE;
double dEqi[3];
double dT;
double dV;
double dAction;
double Pot [space]; //potential energy table
double T; //Inter-bead potential
double V; //Intra-bead potential

const double PI  =3.141592653589793238463;
//double lamb=0.5; //mass = hbar = 1
double lamb=6.0596; //Å^-2*K
double mHe=931.4943228;

bool has_saved=false;
double saved_norm;
unsigned int seed;
bool classical = false;
bool periodic = false;

//Function prototypes

double getRandom();
unsigned int good_seed();
void bisection (double ***);
void centeroid (double ***);
double boxMuller();
double equipartition(double);
double norm (double, double);
double externPot(double);
void potential();
void worldLines(double ***);
double min(double, double);
void swap(double ***);
int factorial (int);
double averageDistance (double ***);
int sumFactorial (int);
double rootMeanSquare (double ***);
void newswap(double ***);
double pbc(double);


MTRand53 mt(good_seed());

int main(int argc, char* argv [])
{
    M=atoi(argv[1]);
    L=atof(argv[2]);
    
    //double *Temp = new double[M+1];
    double ***R=new double** [particles];
    
    for(int kk=0; kk<particles; kk++)
    {
        R[kk]= new double* [M+1];
        
    }
    
    for(int pp=0; pp<particles; pp++)
    {
        for(int kk=0; kk<M+1; kk++)
        {
            R[pp][kk] = new double[dim];
        }
    }
    
    potential(); //Tabulates the potential values
    
    for (int jj=0; jj<jmax; jj++)
    {
        temperature[jj]=1+0.4*jj;
        B=1/temperature[jj];
        dt=B/M;
        
        worldLines(R);  //Initializing the worldlines as straight lines in imaginary time, equally spaced in the box
        T=0;
        V=0;
        
        for(int pp=0; pp<particles; pp++) //calculates the Energy for the system
        {
            
            for(int j=0; j<M; j++)
            {
                V+=externPot(norm(R[pp][j][0], R[pp][j][1]));
                
                for(int ii=0; ii<dim; ii++)
                {
                    
                    int i_plus=j+1;
                    T+=M*pow((R[pp][i_plus][ii]-R[pp][j][ii]),2)/(4*lamb*B*B);
                }
                
                for (int ppp=pp; ppp<particles; ppp++)
                {
                    if(ppp!=pp)
                    {
                        double r0=R[pp][j][0]-R[ppp][j][0];
                        double r1=R[pp][j][1]-R[ppp][j][1];
                        
                        double r_diff = norm(r0,r1);
                        
                        int j_approx = floor(r_diff*ss/L);
                        int j_plus=j_approx+1;
                        V+=(Pot[j_approx]+Pot[j_plus])/(2*M);
                    }
                    
                }
                
                
            }
            
        }
        
        E=dim*particles/(2*dt)-T+V;
        S=-dim*particles/(2*B*dt)+2*T/(B);
        tsamp=-tau;
        D=averageDistance(R);
        RMS=rootMeanSquare(R);
        Esamp=0;
        E2samp=0;
        Ssamp=0;
        acceptance=0;
        swapCount=0;
        Tsamp=0;
        Vsamp=0;
        Dsamp=0;
        RMSsamp=0;
        
       // cout<<T<<" "<<V<<" "<<E<<endl;
        
        for(int t=0; t<tfinal; t++)
        {
            
            //cout<<t<<endl;
            
            for(int pp=0; pp<particles; pp++) //Bisection moves
            {
                for(int kk=0; kk<M; kk++)
                {
                    bisection(R);
                }
            }

            if(tsamp>=0)
                
            {
                if (classical ==false)
                {
                    int dummyint=0;
                    //swap(R);
                    while(dummyint<swaptime)
                    {
                        newswap(R);
                        dummyint++;
                    }
                }

                    
            }
                
            

            if(classical==true)
            {
                centeroid(R);
            }
            
            
            if (tsamp>0)
            {
                
                D=averageDistance(R);
                RMS=rootMeanSquare(R);
                RMSsamp+=RMS;
                Dsamp+=D;
                Esamp+=E;
                E2samp+=(E*E);
                Ssamp+=S;
                Tsamp+=T;
                Vsamp+=V;
                
                /*if (RMS!=RMS || D>L/1.412156237095048801688724)
                {
                    for (int kk=0; kk<M; kk++)
                    {
                        cout<<"R[0]["<<kk<<"]=("<<R[0][kk][0]<<", "<<R[0][kk][1]<<")\tR[1]["<<kk<<"]=("<<R[1][kk][0]<<", "<<R[1][kk][1]<<")"<<endl;
                    }
                    cout<<"-----------\nt="<<t<<"\nD="<<D<<"\nRMS="<<RMS<<endl;
                    return -1;
                }*/
            }
            
            /*double Ttest=0;
            double Vtest=0;
            for(int pp=0; pp<particles; pp++) //calculates the Energy for the system
            {
                
                for(int j=0; j<M; j++)
                {
                    Vtest+=externPot(norm(R[pp][j][0], R[pp][j][1]));
                    
                    for(int ii=0; ii<dim; ii++)
                    {
                        int i_plus=j+1;
                        if (periodic)
                        {
                            double Rnew = (R[pp][i_plus][ii]-R[pp][j][ii])-L*round((R[pp][i_plus][ii]-R[pp][j][ii])/L);
                            Ttest+=M*(pow(Rnew, 2))/(4*lamb*B*B);
                        }
                        else
                        {
                            Ttest+=M*pow((R[pp][i_plus][ii]-R[pp][j][ii]),2)/(4*lamb*B*B);
                        }
                    }
                    
                    for (int ppp=pp; ppp<particles; ppp++)
                    {
                        if(ppp!=pp)
                        {
                            double r0=R[pp][j][0]-R[ppp][j][0];
                            double r1=R[pp][j][1]-R[ppp][j][1];
                            
                            double r_diff = norm(r0,r1);
                            
                            int j_approx = floor(r_diff*ss/L);
                            int j_plus=j_approx+1;
                            Vtest+=(Pot[j_approx]+Pot[j_plus])/(2*M);
                        }
                        
                    }
                    
                    
                }
                
            }
            
            double Etest=dim*particles/(2*dt)-Ttest+Vtest;
            
            if (abs(Etest-E)>1 || E!=E)
            {
                cout<<"Equilibrium energy="<<dim*particles/(2*dt)<<endl;
                cout<<"Variabel\t|\tExakt \t | Inkrementellt\t|"<<endl;
                cout<<"T\t|\t"<<Ttest<<"\t | \t"<<T<<endl;
                cout<<"V\t|\t"<<Vtest<<"\t | \t"<<V<<endl;
                cout<<"E\t|\t"<<Etest<<"\t | \t"<<E<<"\nAccepted moves="<<acceptance<<"\nAccepted Swaps="<<swapCount<<endl;
                cout<<"-----------------"<<endl;
                for (int kk=0; kk<M+1; kk++)
                {
                    cout<<"R[0]["<<kk<<"]=("<<R[0][kk][0]<<", "<<R[0][kk][1]<<" ), |R[0]|="<<norm(R[0][kk][0],R[0][kk][1])
                    <<"\tR[0]["<<kk<<"]=("<<R[1][kk][0]<<", "<<R[1][kk][1]<<"), |R[1]|="<<norm(R[1][kk][0],R[1][kk][1])<<endl;
                }
                cout<<"-----------------"<<endl;
                for(int pp=0; pp<particles; pp++)
                {
                    cout<<"Worldine "<<pp<<"---->"<<permutation[pp]<<endl;
                }
                return 9;
            }*/
            
            tsamp++;
        }
        
        
        
        RMSav[jj]=RMSsamp/samp;
        Dav[jj]=Dsamp/samp;
        Tav[jj]=Tsamp/samp;
        Vav[jj]=Vsamp/samp;
        Eav[jj]=Esamp/samp;
        E2av[jj]=E2samp/samp;
        Sav[jj]=Ssamp/samp;
        Cv[jj]=B*B*(E2av[jj]-Sav[jj]-Eav[jj]*Eav[jj]);
        ACCEPT[jj]=(double)acceptance/(double)(M*particles*tfinal);
        SWAP[jj]=(double)swapCount/(double)(tsamp*swaptime);
        //cout<<"Acceptance:"<<acceptance<<"/"<<M*particles*tfinal<<"="<<ACCEPT[jj]<<endl;
        //cout<<"Swap ratio:"<<swapCount<<"/"<<tsamp*swaptime<<"="<<SWAP[jj]<<"\n"<<endl;
    }
   
    
    for(int k=0; k<jmax; k++)
    {
        //cout<<temperature[k]<<" | "<<Eav[k]/particles<<" | "<<E2av[k]/particles<<" | "<<Cv[k]/particles<<" | "<<Dav[k]<<" | "<<RMSav[k]<<" | "<<ACCEPT[k]<<" | "<<SWAP[k]<<"\n";
        cout<<temperature[k]<<" "<<Eav[k]/particles<<" "<<E2av[k]/particles<<" "<<Cv[k]/particles<<" "<<Dav[k]<<" "<<RMSav[k]<<" "<<ACCEPT[k]<<" "<<SWAP[k]<<" "<<L-10<<" "<<particles<<"\n";

    }
    
   /* for(int k=0; k<particles; k++)
    {
        cout<<k<<"----->"<<permutation[k]<<endl;
        cout<<"Timeslice:"<<0<<"\t"<<"("<<R[k][0][0]<<", "<<R[k][0][1]<<
        ")\n"<<"Timeslice:"<<M<<"\t("<<R[k][M][0]<<", "<<R[k][M][1]<<")\n"<<endl;


        
    }*/
    
    
    for(int pp=0; pp<particles; pp++)
    {
        for(int kk=0; kk<M+1; kk++)
        {
            delete[] R[pp][kk];
        }
    }
    
    for(int kk=0; kk<particles; kk++)
    {
        delete[] R[kk];
        
    }
    
    delete [] R;
    //delete [] Temp;
    
    return 0;
}


//Functions


inline double averageDistance(double ***R) //Average distance between beads
{
    double dist=0;
    for(int pp=0; pp<particles; pp++)
    {
        for(int kk=0; kk<M; kk++)
        {
            double Rx=R[pp][kk][0];
            double Ry=R[pp][kk][1];
            for (int ppp=pp+1; ppp<particles; ppp++)
            {
                double xComp=Rx-R[ppp][kk][0];
                double yComp=Ry-R[ppp][kk][1];
                dist+=norm(xComp,yComp)/(M*sumFactorial(particles-1));
            }
        }
    }
    return dist;
}

void bisection(double *** R) //Bisection method for moving world-line segments
{
    double **Rtemp = new double *[length];  //Saving the old coordinates
    double **Rtrial = new double *[length]; //Trial array for new coordinates

    int startBead = floor(getRandom()*M);   //time slice of first bead
    int p= floor(getRandom()*particles); //particle
    bool included [length]; //Keeps track if the bead has been updated
    bool hasV [length]; //Keeps track if the bead's potential energy change has been calculated
    int stepsize=(length-1)/2;
    
    for (int pp=0; pp<length; pp++)
    {
        Rtemp[pp]= new double[dim];
        Rtrial[pp]=new double[dim];
    }
    
    for (int ii=0; ii<dim; ii++)
    {
        for (int kkk=0; kkk<length; kkk++)
        {
            int bead = startBead+kkk;
            if (bead<M)
            {
                Rtemp[kkk][ii]=R[p][bead][ii];
                Rtrial[kkk][ii]=R[p][bead][ii];
            }
            else
            {
                int perBead=bead-floor(bead/M)*M;
                Rtemp[kkk][ii]=R[permutation[p]][perBead][ii];
                Rtrial[kkk][ii]=R[permutation[p]][perBead][ii];
            }
            
        }
    }
    
    for (int ii=0; ii<length; ii++)
    {
        included[ii]=false;
        hasV[ii]=false;
    }
    included[0]=included[length-1]=true;
    hasV[0]=hasV[length-1]=true;
    
    dT=0;
    dV=0;
    bool run = true;
    dAction=0;
    Vtemp=0;
    
    while (run)
    {
        for (int pp=0; pp<length; pp+=stepsize)
        {
            if (included[pp]==false)
            {
                for (int ii=0; ii<dim; ii++)
                {

                    double eta=boxMuller()*sqrt(lamb*stepsize*dt);
                    if (periodic)
                    {
                        double R1=Rtemp[pp][ii]-Rtrial[pp-stepsize][ii];
                        double R2=Rtrial[pp+stepsize][ii]-Rtemp[pp][ii];
                        R1-=L*round(R1/L);
                        R2-=L*round(R2/L);
                        
                        Rtrial[pp][ii]=(2*Rtemp[pp][ii]-R1+R2)/2+eta;
                        
                        if (abs(Rtrial[pp][ii])>(L/2))
                        {
                            Rtrial[pp][ii]-=L*round(Rtrial[pp][ii]/L);
                        }

                        
                    }
                    else
                    {
                        Rtrial[pp][ii]=(Rtrial[pp-stepsize][ii]+Rtrial[pp+stepsize][ii])/2+eta;
                    }
 //                   cout<<"R["<<pp<<"]["<<ii<<"]="<<Rtemp[pp][ii]<<"\neta="<<eta<<"\nR'["<<pp<<"]["<<ii<<"]="<<Rtrial[pp][ii]<<"\n-------"<<endl;
                }
                included[pp]=true;
            }
        }
        
        for (int pp=0; pp<length; pp+=stepsize)
        {
            if (included[pp]==true && hasV[pp]==false)
            {
                for (int ppp=0; ppp<particles; ppp++)
                {
                    
                    if ((startBead+pp)<M)
                    {
                        if(ppp!=p)
                        {
                            double r_diff[dim];
                            double r_trial_diff[dim];
                            double rDiff;
                            double rTrialDiff;
                            
                            for (int ii=0; ii<dim; ii++)
                            {
                                int bead = startBead+pp;
                                r_diff[ii] = Rtemp[pp][ii]-R[ppp][bead][ii];
                                r_trial_diff[ii] = Rtrial[pp][ii]-R[ppp][bead][ii];
                            }
                            
                            rDiff=norm(r_diff[0], r_diff[1]);
                            rTrialDiff=norm(r_trial_diff[0], r_trial_diff[1]);
                            
                            int i_approx = floor(rDiff*ss/L);
                            int i_plus = i_approx+1;
                            int j_approx = floor(rTrialDiff*ss/L);
                            int j_plus=j_approx+1;
                            dV+=((Pot[j_approx]+Pot[j_plus])-(Pot[i_approx]+Pot[i_plus]))/(2*M);
                        }
                    }
                    
                    else if((startBead+pp)>=M)
                    {
                        
                        if (ppp!=permutation[p])
                        {
                            double r_diff[dim];
                            double r_trial_diff[dim];
                            double rDiff;
                            double rTrialDiff;
                            
                            for (int ii=0; ii<dim; ii++)
                            {
                                int bead = startBead+pp-floor((startBead+pp)/M)*M;
                                r_diff[ii] = Rtemp[pp][ii]-R[ppp][bead][ii];
                                r_trial_diff[ii] = Rtrial[pp][ii]-R[ppp][bead][ii];
                            }
                            
                            rDiff=norm(r_diff[0], r_diff[1]);
                            rTrialDiff=norm(r_trial_diff[0], r_trial_diff[1]);
                            
                            int i_approx = floor(rDiff*ss/L);
                            int i_plus = i_approx+1;
                            int j_approx = floor(rTrialDiff*ss/L);
                            int j_plus=j_approx+1;
                            dV+=((Pot[j_approx]+Pot[j_plus])-(Pot[i_approx]+Pot[i_plus]))/(2*M);
                        }
                    }
                    
                }
                dV+=externPot(norm(Rtrial[pp][0],Rtrial[pp][1]));
                //dV+=externPot(norm(Rtrial[pp][0],Rtrial[pp][1]))-externPot(norm(Rtemp[pp][0],Rtemp[pp][1]));


                //cout<<dV<<endl;
                hasV[pp]=true;
            }
        }
        double ran= getRandom();
        double dA=-B*stepsize*dV+dAction;
        double A=exp(dA);
        
        if (A<ran)
        {
            run=false;
        }
        
        else
        {
       //     cout<<"\nStepsize ="<<stepsize<<"\ndV="<<dV<<"\nA="<<A<<"\nran="<<ran<<endl;

         //   cout<<Rtemp[stepsize][0]<<", "<<Rtemp[stepsize][1]<<"---->"<<Rtrial[stepsize][0]<<", "<<Rtrial[stepsize][1]<<endl;
            if (stepsize>1)
            {
                dAction = B*stepsize*dV;
                stepsize/=2;
            }
            else
            {
                for(int pp=0; pp<length-1; pp++)
                {
                    for (int ii=0; ii<dim; ii++)
                    {
                        if (periodic)
                        {
                            double Rnew = (Rtrial[pp+1][ii]-Rtrial[pp][ii]);
                            double Rold = (Rtemp[pp+1][ii]-Rtemp[pp][ii]);
                            dT+=M*(pow(pbc(Rnew), 2)
                                   -pow(pbc(Rold),2))/(4*lamb*B*B);
                        }
                        
                        else
                        {
                            dT+=M*(pow((Rtrial[pp+1][ii]-Rtrial[pp][ii]), 2)
                             -pow((Rtemp[pp+1][ii]-Rtemp[pp][ii]),2))/(4*lamb*B*B);
                        }
                    }
                }
                
                E+=-dT+dV;
                S+=2*dT/B;
                T+=dT;
                V+=dV;
                //cout<<dT/M<<"\n"<<dV/M<<endl;
                //cout<<T<<"\n"<<V<<"\n"<<endl;
                acceptance++;
                
                for (int ii=0; ii<dim; ii++)
                {
                    for (int pp=0; pp<length; pp++)
                    {
                        int bead = startBead+pp;
                        if(bead<M)
                        {
                            R[p][bead][ii]=Rtrial[pp][ii];
                        }
                        else
                        {
                            int perBead=bead-floor(bead/M)*M;
                            R[permutation[p]][perBead][ii]=Rtrial[pp][ii];
                        }
                    }
                    if ((startBead+length>M+1))
                    {
                        R[p][M][ii]=R[permutation[p]][0][ii];
                    }
                }
                
                run=false;
                
            }
        }
    }
    
    
    
    for (int pp=0; pp<length; pp++)
    {
        delete[] Rtemp[pp];
        delete[] Rtrial[pp];
    }
    
    delete [] Rtrial;
    delete [] Rtemp;
    return;
    
}

double boxMuller() //Box-Mueller tranform, casts 2 uncorrelated uniformly distributed random numbers on 2 uncorrelated Gaussian numbers
{
    
    double x,y,r, r2;
    
    if(has_saved){
        has_saved=false;
        return saved_norm;
    }
    
    else{
        
        do{
            x=2.0*getRandom()-1;
            y=2.0*getRandom()-1;
            r2=x*x+y*y;
        }while(r2>1.0 || r2==0.0);
        
        r=sqrt(-2.0*log(r2)/r2);
        saved_norm=r*y; /* save second random number */
        has_saved=true;
        return r*x; /* return first random number */
    }
}

void centeroid(double *** R) //Centeroid movements of the world-line
{
    double **R_mean=new double *[particles];
    double **R_trial=new double *[M];
    double rk [dim];
    
    for(int pp=0; pp<particles; pp++)
    {
        R_mean[pp]=new double [dim];
        
    }
    
    for (int kk=0; kk<M; kk++)
    {
        R_trial[kk]=new double [dim];
    }
    
    
    for(int pp=0; pp<particles; pp++) //Centeroid moves
    {
        
        for(int ii=0; ii<dim; ii++)
        {
            R_mean[pp][ii]=0;
            for(int kk=0; kk<M; kk++)
            {
                R_mean[pp][ii]+=R[pp][kk][ii]/M;
            }
        }
    }
    
    for(int n=0; n<particles; n++)
    {
        dV=0;
        int p=floor(getRandom()*particles);
        double eta[dim];
        
        for(int ii=0; ii<dim; ii++)
        {
            eta[ii]=boxMuller()*sqrt((lamb*dt)/M);
            R_mean_trial[ii]=R_mean[p][ii]+eta[ii];
        }
        for(int pp=0; pp<particles; pp++)
        {
            if(pp!=p)
            {
                double rx=R_mean[p][0]-R_mean[pp][0];
                double ry=R_mean[p][1]-R_mean[pp][1];
                double rx_trial=R_mean_trial[0]-R_mean[pp][0];
                double ry_trial=R_mean_trial[1]-R_mean[pp][1];
                
                R_diff=norm(rx, ry);
                R_trial_diff=norm(rx_trial, ry_trial);
                
                int i_approx=floor(R_diff*ss/L);
                int i_plus=i_approx+1;
                int i_trial_approx=floor(R_trial_diff*ss/L);
                int i_trial_plus=i_trial_approx+1;
                
                dV+=(Pot[i_trial_approx]+Pot[i_trial_plus])/2-(Pot[i_approx]+Pot[i_plus])/2;
            }
        }
        
        //Check if any bead is moved outside the wall
        for(int kk=0; kk<M; kk++)
        {
            rk[0]=rk[1]=0;
            for (int ii=0; ii<dim; ii++)
            {
                R_trial[kk][ii]=R[p][kk][ii]+eta[ii];
                rk[ii]=R[p][kk][ii]+eta[ii];
            }
            dV+=externPot(norm(rk[0], rk[1]));
        }
        
        if(dV<0)
        {
            Vtemp=0;
            
            for(int kk=0; kk<M; kk++)
            {
                for(int pp=0;pp<particles;pp++)
                {
                    if(pp!=p)
                    {
                        double rx=R[p][kk][0]-R[pp][kk][0];
                        double ry=R[p][kk][1]-R[pp][kk][1];
                        double rx_trial=R_trial[kk][0]-R[pp][kk][0];
                        double ry_trial=R_trial[kk][1]-R[pp][kk][1];
                        R_diff=norm(rx,ry);
                        R_trial_diff=norm(rx_trial,ry_trial);
                        
                        int i_approx=floor(R_diff*ss/L);
                        int i_plus=i_approx+1;
                        int i_trial_approx=floor(R_trial_diff*ss/L);
                        int i_trial_plus=i_trial_approx+1;
                        
                        Vtemp+=(Pot[i_trial_approx]+Pot[i_trial_plus])/(2*M)-(Pot[i_approx]+Pot[i_plus])/(2*M);
                    }
                }
                for (int ii=0; ii<dim; ii++)
                {
                    R[p][kk][ii]=R_trial[kk][ii];
                    
                    
                    if (kk==0)
                    {
                        R[p][M][ii]=R[p][0][ii];
                    }
                    
                }
                
            }
            V+=Vtemp;
            E+=Vtemp;
            
        }
        
        
        else
        {
            double p_ran=getRandom();
            A=exp(-dt*dV);
            if (A>p_ran)
            {
                Vtemp=0;
                for(int kk=0; kk<M; kk++)
                {
                    
                    for(int pp=0;pp<particles;pp++)
                    {
                        if(pp!=p)
                        {
                            double rx=R[p][kk][0]-R[pp][kk][0];
                            double ry=R[p][kk][1]-R[pp][kk][1];
                            double rx_trial=R_trial[kk][0]-R[pp][kk][0];
                            double ry_trial=R_trial[kk][1]-R[pp][kk][1];
                            R_diff=norm(rx,ry);
                            R_trial_diff=norm(rx_trial,ry_trial);
                            
                            int i_approx=floor(R_diff*ss/L);
                            int i_plus=i_approx+1;
                            int i_trial_approx=floor(R_trial_diff*ss/L);
                            int i_trial_plus=i_trial_approx+1;
                            
                            
                            Vtemp+=(Pot[i_trial_approx]+Pot[i_trial_plus])/(2*M)-(Pot[i_approx]+Pot[i_plus])/(2*M);
                        }
                    }
                    
                    for (int ii=0; ii<dim; ii++)
                    {
                        R[p][kk][ii]=R_trial[kk][ii]-L*round(R_trial[kk][ii]/L);
                        
                        
                        if (kk==0)
                        {
                            R[p][M][ii]=R[p][0][ii];
                        }
                        
                    }
                    
                }
                V+=Vtemp;
                E+=Vtemp;
                
            }
        }
    }
    
    
    for(int kk=0; kk<particles; kk++)
    {
        delete[] R_mean[kk];
        
    }
    
    for (int kk=0; kk<M; kk++)
    {
        delete[] R_trial[kk];
    }
    
    delete[] R_trial;
    delete[] R_mean;
}

double equipartition(double tau) //Calculates the Equipartition part of the kinetic energy as a function of the timestep 'tau'
{
    double eq= particles*M*log(4*PI*tau*lamb)/(2*tau);
    return eq;
}

inline double externPot(double R) //Adds an externalpotential to the box
{
    double potential;
 
    if (periodic)
    {
        return potential =0;
    }
    //Circular domain, hard walls
    if (abs(R)>(0.5*(L-10)))
    {
        return potential=pow(10.0, 1000000);
    }
    else
    {
        return potential=0;
    }
    
    //Harmonic trap
    //return potential=pow(R,2)/2;
    
    //Free particles
    //return potential=0.0;
}

inline int factorial (int x) //factorial of a number
{
    int value=1;
    while (x>1)
    {
        value=x*value;
        x--;
    }
    return value;
}

inline double getRandom() //Creates a random number using MT-twister
{
    double newRandom=mt();
    return newRandom;
}

unsigned int good_seed() //Chooses a good seed for the random generator
{
    unsigned int random_seed, random_seed_a, random_seed_b;
    std::ifstream file ("/dev/urandom", std::ios::binary);
    if (file.is_open())
    {
        char * memblock;
        int size = sizeof(int);
        memblock = new char [size];
        file.read (memblock, size);
        file.close();
        random_seed_a = *reinterpret_cast<int*>(memblock);
        delete[] memblock;
    }
    else
    {
        random_seed_a = 0;
    }
    //random_seed_b = std::time(0);
    random_seed_b = time(0);
    random_seed = random_seed_a xor random_seed_b;
    /*std::cout << "random_seed_a = " << random_seed_a << std::endl;
     std::cout << "random_seed_b = " << random_seed_b << std::endl;
     std::cout << " random_seed =  " << random_seed << std::endl;*/
    return random_seed;
}

double min(double arg1, double arg2)
{
    if (arg1==NAN)
    {
        return 0;
    }
    
    else if (arg2==NAN)
    {
        return 0;
    }
    
    else if (arg1<arg2)
    {
        return arg1;
    }
    
    else
    {
        return arg2;
    }
}

inline double norm(double Rx, double Ry) //Calculates the R2-norm
{
    double NORM;
    double R1[2]={Rx,Ry};
    if (periodic)
    {
        for (int ii=0; ii<dim; ii++)
        {
            double Rnew=pbc(R1[ii]);
            R1[ii]=Rnew;

        }
    }
    
    return NORM = sqrt(abs(R1[0]*R1[0]+R1[1]*R1[1]));
}

double pbc(double Rx)
{
    double newR;
    return newR=Rx-L*round(Rx/L);
}

inline void potential() //Tabulates the potential interaction between particles for different distances
{
    for(int i=0; i<space; i++)
    {
        r[i]=i*(L)/(ss*2.9673);
        if(r[i]<d)
        {
            h[i]=exp(-pow(d/r[i]-1,2));
        }
        else
        {
            h[i]=1;
        }
        Pot[i]=0; //non-interacting particles
        //Pot[i]=pow(r[i],2)/2;//Harmonic interaction potential
        //Pot[i]=2*(pow((L/particles)/(r[i]),12)-2*pow((L/particles)/(r[i]),6)); //Lennard - Jones potential
        
       /*Pot[i]=10.8*(0.54485046e6*exp(-13.353384*r[i])
                     -(1.3732412*pow(r[i],-6)
                       +0.4253785*pow(r[i],-8)
                       +0.1781*pow(r[i],-10))*h[i]);  //Aziz potential*/
        
        
    }
    //Pot[0]=100000000000000000;
}

double rootMeanSquare(double ***R)
{
    double rootMean;
    for(int pp=0; pp<particles; pp++)
    {
        for(int kk=0; kk<M; kk++)
        {
            double Rx=R[pp][kk][0];
            double Ry=R[pp][kk][1];
            for (int ppp=pp+1; ppp<particles; ppp++)
            {
                double xComp=Rx-R[ppp][kk][0];
                double yComp=Ry-R[ppp][kk][1];
                double Rnorm= norm(xComp, yComp);
                rootMean+=Rnorm*Rnorm;
            }
        }
    }
    
    rootMean/=M*sumFactorial(particles-1);
    rootMean=sqrt(rootMean);
    return rootMean;
}

void swap(double ***R) //swap move of two worldlines
{
    int startBead;   //time slice of first bead
    int p = floor(getRandom()*particles);
    bool included [length];
    bool hasV [length];
   // bool run = false;
    bool run = true;
    int stepsize=(length-1)/2;
    int p2; //Other particle

    
    startBead=floor(getRandom()*M);

    do
    {
        p2=floor(getRandom()*particles);
    }while(p2==p);
    
    
    double ***Rtrial = new double **[2];
    double ***Rtemp = new double **[2];
    
    for (int pp=0; pp<2; pp++)
    {
        Rtrial[pp]=new double *[length];
        Rtemp[pp]=new double *[length];
        
    }
    
    for (int pp=0; pp<2; pp++)
    {
        for (int kk=0; kk<length; kk++)
        {
            Rtrial[pp][kk]=new double [dim];
            Rtemp[pp][kk]=new double [dim];
        }
    }
    //Save the initial configuration
    for (int kk=0; kk<length; kk++)
    {
        for (int ii=0; ii<dim; ii++)
        {
            int kkk=startBead+kk;
            
            if (kkk>M)
            {
                int permkkk=kkk-floor(kkk/M)*M;
                Rtemp[0][kk][ii]=R[permutation[p]][permkkk][ii];
                Rtemp[1][kk][ii]=R[permutation[p2]][permkkk][ii];

            }
            
            else
            {
                Rtemp[0][kk][ii]=R[p][kkk][ii];
                Rtemp[1][kk][ii]=R[p2][kkk][ii];
            }
        }
        included[kk]=false;
        hasV[kk]=false;
    }

    //Setup the trial paths for the new bead paths
    for (int ii=0; ii<dim; ii++)
    {
        Rtrial[0][0][ii]=R[p][startBead][ii];
        Rtrial[1][0][ii]=R[p2][startBead][ii];
        
        if (startBead+length>M)
        {
            int endBead=startBead+length-1-floor((startBead+length-1)/M)*M;
            Rtrial[0][length-1][ii]=R[permutation[p2]][(endBead)][ii];
            Rtrial[1][length-1][ii]=R[permutation[p]][(endBead)][ii];
        }
        else
        {
            Rtrial[0][length-1][ii]=R[p2][(startBead+length-1)][ii];
            Rtrial[1][length-1][ii]=R[p][(startBead+length-1)][ii];
        }
        
    }
    included[0]=included[length-1]=hasV[0]=hasV[length-1]=true;
    
    double dAction=0;
    dT=0;
    dV=0;
    while(run)
    {
        //Move the trial bead(s)
        for(int kk=0; kk<length; kk+=stepsize)
        {
            if (included[kk]==false)
            {
                
                for (int pp=0; pp<2; pp++)
                {
            
                    for (int ii=0; ii<dim; ii++)
                    {
                        Rtrial[pp][kk][ii]=(R[pp][kk-stepsize][ii]+R[pp][kk+stepsize][ii])/2+boxMuller()*sqrt(lamb*stepsize*dt);
                    }
                           
                }
                included[kk]=true;
            }
            
                           
        }
        //Calculate action cost of move
        
        for (int kk=0; kk<(length); kk+=stepsize)
        {
            int currentBead=startBead+kk;
            
            if (included[kk]==true && hasV[kk]==false)
            {
                for (int ppp=0; ppp<particles; ppp++)
                {

                    if((currentBead<M) && (ppp!=p))
                    {
                        
                        if(ppp==p2)
                        {
                            double r_diff[dim];
                            double r_trial_diff[dim];
                            double rDiff;
                            double rTrialDiff;
                        
                            for (int ii=0; ii<dim; ii++)
                            {
                                r_diff[ii] = Rtemp[0][kk][ii]-Rtemp[1][kk][ii];
                                r_trial_diff[ii] = Rtrial[0][kk][ii]-Rtrial[1][kk][ii];
                            }
                        
                            rDiff=norm(r_diff[0], r_diff[1]);
                            rTrialDiff=norm(r_trial_diff[0], r_trial_diff[1]);
                        
                            int i_approx = floor(rDiff*ss/L);
                            int i_plus = i_approx+1;
                            int j_approx = floor(rTrialDiff*ss/L);
                            int j_plus=j_approx+1;
                            dV+=((Pot[j_approx]+Pot[j_plus])-(Pot[i_approx]+Pot[i_plus]))/(2*M);
                        }
                    
                        else
                        {
                       
                            for (int pp=0; pp<2; pp++)
                            {
                            
                                double r_diff[dim];
                                double r_trial_diff[dim];
                                double rDiff;
                                double rTrialDiff;
                        
                                for (int ii=0; ii<dim; ii++)
                                {
                                    int bead = (startBead+kk+M)%M;
                                    r_diff[ii] = Rtemp[pp][kk][ii]-R[ppp][bead][ii];
                                    r_trial_diff[ii] = Rtrial[pp][kk][ii]-R[ppp][bead][ii];
                                }
                        
                                rDiff=norm(r_diff[0], r_diff[1]);
                                rTrialDiff=norm(r_trial_diff[0], r_trial_diff[1]);
                        
                                int i_approx = floor(rDiff*ss/L);
                                int i_plus = i_approx+1;
                                int j_approx = floor(rTrialDiff*ss/L);
                                int j_plus=j_approx+1;
                                dV+=((Pot[j_approx]+Pot[j_plus])-(Pot[i_approx]+Pot[i_plus]))/(2*M);
                            }
                            
                        }
                    }
                    else if ((currentBead>(M-1)) && (ppp !=permutation[p]))
                    {
                        if(ppp==permutation[p2])
                        {
                            double r_diff[dim];
                            double r_trial_diff[dim];
                            double rDiff;
                            double rTrialDiff;
                            
                            for (int ii=0; ii<dim; ii++)
                            {
                                r_diff[ii] = Rtemp[0][kk][ii]-Rtemp[1][kk][ii];
                                r_trial_diff[ii] = Rtrial[0][kk][ii]-Rtrial[1][kk][ii];
                            }
                            
                            rDiff=norm(r_diff[0], r_diff[1]);
                            rTrialDiff=norm(r_trial_diff[0], r_trial_diff[1]);
                            
                            int i_approx = floor(rDiff*ss/L);
                            int i_plus = i_approx+1;
                            int j_approx = floor(rTrialDiff*ss/L);
                            int j_plus=j_approx+1;
                            dV+=((Pot[j_approx]+Pot[j_plus])-(Pot[i_approx]+Pot[i_plus]))/(2*M);
                            
                        }
                        
                        else
                        {
                            
                            for (int pp=0; pp<2; pp++)
                            {
                                
                                double r_diff[dim];
                                double r_trial_diff[dim];
                                double rDiff;
                                double rTrialDiff;
                                
                                for (int ii=0; ii<dim; ii++)
                                {
                                    int bead = (startBead+kk+M)%M;
                                    r_diff[ii] = Rtemp[pp][kk][ii]-R[ppp][bead][ii];
                                    r_trial_diff[ii] = Rtrial[pp][kk][ii]-R[ppp][bead][ii];
                                }
                                
                                rDiff=norm(r_diff[0], r_diff[1]);
                                rTrialDiff=norm(r_trial_diff[0], r_trial_diff[1]);
                                
                                int i_approx = floor(rDiff*ss/L);
                                int i_plus = i_approx+1;
                                int j_approx = floor(rTrialDiff*ss/L);
                                int j_plus=j_approx+1;
                                dV+=((Pot[j_approx]+Pot[j_plus])-(Pot[i_approx]+Pot[i_plus]))/(2*M);
                            }
                            
                        }
                        
                    }
                }
                
                for (int pp=0; pp<2; pp++)
                {
                    dV+=externPot(norm(Rtrial[pp][kk][0],Rtrial[pp][kk][1]));
                }
                
                hasV[kk]=true;
            }
            
        }
        //Test the action change
        double ran= getRandom();
        double dA=-stepsize*(dV)+dAction;

        double A=min(exp(dA),1.0);
        //       cout<<"\ndV="<<dV<<"\nA="<<A<<"\nran="<<ran<<endl;
        
        if (A<ran)
        {
            run=false;
        }
        
        else
        {
            //cout<<stepsize<<endl;
            if (stepsize>1)
            {
                dAction = stepsize*dV;
                stepsize/=2;
            }
            else
            {
                for(int pp=1; pp<length; pp++)
                {
                    for (int ii=0; ii<dim; ii++)
                    {
                        dT+=(pow((Rtrial[0][pp][ii]-Rtrial[0][pp-1][ii])/(dt), 2)
                             -pow((Rtemp[0][pp][ii]-Rtemp[0][pp-1][ii])/(dt),2))/(4*lamb*M);
                        dT+=(pow((Rtrial[1][pp][ii]-Rtrial[1][pp-1][ii])/(dt), 2)
                             -pow((Rtemp[1][pp][ii]-Rtemp[1][pp-1][ii])/(dt),2))/(4*lamb*M);
                    }
                }
                
                E+=-dT+dV;
                S+=2*dT/B;
                T+=dT;
                V+=dV;
                //cout<<T<<"\n"<<V<<"\n"<<endl;
                swapCount++;
                
                
                
                //swap all follwowing worldline coordinates with each other
                
                int pQuicksave=permutation[p];
                permutation[p]=permutation[p2];
                permutation[p2]=pQuicksave;
                
                for (int ii=0; ii<dim; ii++)
                {
                    for (int kk=0; kk<length; kk++)
                    {
                        int currentBead=startBead+kk;
                        if (currentBead>(M-1))
                        {
                            int perCurrentBead=currentBead-floor(currentBead/M)*M;
                            R[permutation[p]][perCurrentBead][ii]=Rtrial[0][kk][ii];
                            R[permutation[p2]][perCurrentBead][ii]=Rtrial[1][kk][ii];
                        }
                        
                        else
                        {
                            R[p][currentBead][ii]=Rtrial[0][kk][ii];
                            R[p2][currentBead][ii]=Rtrial[1][kk][ii];
                        }
                        
                    }
                    
                    
                    
                    double Rquicksave;
                    for (int kk=startBead+length; kk<M; kk++)
                    {
                        
                        Rquicksave=R[p][kk][ii];
                        R[p][kk][ii]=R[p2][kk][ii];
                        R[p2][kk][ii]=Rquicksave;
                        
                    }
                    
                   /*Rquicksave=R[permutation[p]][M][ii];
                    R[permutation[p]][M][ii]=R[permutation[p2]][M][ii];
                    R[permutation[p2]][M][ii]=Rquicksave;//*/
                    
                }
              
                
                for (int ii=0; ii<dim; ii++)
                {
                    R[p][M][ii]=R[permutation[p]][0][ii];
                    R[p2][M][ii]=R[permutation[p2]][0][ii];
                }
                
                run=false;
                
            }
        }

        
    }
    
    for (int pp=0; pp<2; pp++)
    {
        for (int kk=0; kk<length; kk++)
        {
            delete[] Rtrial[pp][kk];
            delete[] Rtemp[pp][kk];
            
        }
    }
    
    for (int pp=0; pp<2; pp++)
    {
        delete[] Rtrial[pp];
        delete[] Rtemp[pp];
        
    }
    
    delete[] Rtrial;
    delete[] Rtemp;
    return;
}

int sumFactorial (int x) //factorial of a number
{
    int value=0;
    while (x>0)
    {
        value+=x;
        x--;
    }
    return value;
}

inline void worldLines(double *** R) //Initiates the world-lines
{
    for(int pp=0; pp<particles; pp++)
    {
        for(int i=0;i<M;i++)
        {
            //R[pp][i][0]=2.9673*(pow(-1, floor(pp/2))*(0.5+floor(pp/4)));
            R[pp][i][0]=R[pp][i][1]=0;                             //Harmonic potential start
            //R[pp][i][1]=-L/sqrt(32)+2*pp*L/sqrt(32     //Lennard-Jones start
            //R[pp][i][1]=2.9673*pow(-1,pp)/2;                 //Aziz-potential start
            
           // cout<<R[pp][i][0]<<"\n"<<R[pp][i][1]<<endl;
            
        }
        R[pp][M][0]=R[pp][0][0];
        R[pp][M][1]=R[pp][0][1];
        permutation[pp]=inversePermutation[pp]=pp;
    }
}

void newswap(double ***R) //swap move of two worldlines
{
    //General variables
    int startBead=floor(getRandom()*(M));   //time slice of first bead
    int p = floor(getRandom()*particles); //First particle in permutation cycle
    int stepsize=(length-1)/2;

    //Booleans for while loops
    bool run = true;
    bool permute = true;
    
    //Variables for permutations
    double probability [particles]; //probability intervall for permutations
    double lastProbability=0; //
    double particleTable[particles]; // Table of free-energy density matricies
    int particle_goes_to[particles]; //List over with which each particle switches with
    int permutationTable[particles]; //Keeps track of permutation
    int p2=0; //second particle
    double Q=0; //Normalizing factor in permutation transition
    vector<int> particlesUsed(1,NAN); //Table of included particle in the permutation cycle
    int numberOfParticlesUsed=0;
    bool includedInCycle[particles];
    
    
    //Variables for wiggle move
    bool included [length];
    bool hasV [length];
    
    //initialize tables
    for (int pp=0; pp<particles; pp++)
    {
        particleTable[pp]=0;
        particle_goes_to[pp]=pp;
        permutationTable[pp]=permutation[pp];
        includedInCycle[pp]=false;
    }
    
    //Create table of free energy matricies
    int endBead=startBead+length-M*(int)floor((startBead+length)/M);
    for (int pp=0; pp<particles; pp++)
    {
        
        particleTable[pp]=pow(4*PI*lamb*length*dt, dim/2);

        if (p!=pp)
        {
            double exponent=0;
            for (int ii=0; ii<dim; ii++)
            {
                if (startBead>(M-length))
                {
                    if (periodic)
                    {
                        double R1=R[p][startBead][ii]-R[permutation[pp]][endBead][ii];
                        exponent-=(pow(R1-L*round(R1/L), 2))/(4*lamb*length*dt);
                    }
                    else
                    {
                        exponent-=(pow(R[p][startBead][ii]-R[permutation[pp]][endBead][ii], 2))/(4*lamb*length*dt);
                    }
                }
                else
                {
                    if (periodic)
                    {
                        double R1=R[p][startBead][ii]-R[pp][endBead][ii];
                        exponent-=(pow(R1-L*round(R1/L), 2))/(4*lamb*length*dt);
                    }
                    else
                    {
                        exponent-=(pow(R[p][startBead][ii]-R[pp][endBead][ii], 2))/(4*lamb*length*dt);
                    }
                }
            }
            particleTable[pp]*=exp(exponent);
            
        }
        
        else
        {
            particleTable[pp]=0; //Excludes the identity permutation
        }
        probability[pp]=particleTable[pp]+lastProbability;
        lastProbability=probability[pp];
        
        Q+=particleTable[pp];

    }
    
    //Testing if permutation is accepted
    double testExponent=0;
    
    for (int ii=0; ii<dim; ii++)
    {
        testExponent-=(pow(R[p][startBead][ii]-R[p][endBead][ii], 2))/(4*lamb*length*dt); //Free energy matrix of current configuration for the line
    }
    
    double denom= Q+pow(4*PI*lamb*length*dt, dim/2)*exp(testExponent);
    double C= Q/denom;
    double p_rand=getRandom();
    if (C<p_rand || C!=C)
    {
        return;
    }
    
    else
    {
        //Randomizing a partner for the permutation
        double rand=getRandom();
        for (int pp=0; pp<particles; pp++)
        {
            probability[pp]/=Q;
            if (rand<probability[pp])
            {
                p2=pp;
                particle_goes_to[p]=p2;
                permutationTable[p]=permutation[p2];
                break;
            }
        }
        
        
        particlesUsed[0]=p2;
        numberOfParticlesUsed++;
        includedInCycle[p2]=true;
        //Adding more particles to the permutation cycle

        while (permute)
        {
            lastProbability=0;
            Q=0;
            for (int pp=0; pp<particles; pp++)
            {
                
                particleTable[pp]=pow(4*PI*lamb*length*dt, dim/2);
                bool notAllowedParticle=false;
                
                for (int ppp=0; ppp<numberOfParticlesUsed; ppp++)
                {
                    if (particlesUsed[ppp]==pp)
                    {
                        notAllowedParticle=true;
                    }
                }
                
                if (notAllowedParticle==false)
                {
                    double exponent=0;
                    for (int ii=0; ii<dim; ii++)
                    {
                        if (startBead>(M-length))
                        {
                            if (periodic)
                            {
                                double R1=R[p][startBead][ii]-R[permutation[pp]][endBead][ii];
                                exponent-=(pow(R1-L*round(R1/L), 2))/(4*lamb*length*dt);
                            }
                            else
                            {
                                exponent-=(pow(R[p][startBead][ii]-R[permutation[pp]][endBead][ii], 2))/(4*lamb*length*dt);
                            }
                        }
                        else
                        {
                            if (periodic)
                            {
                                double R1=R[p][startBead][ii]-R[pp][endBead][ii];
                                exponent-=(pow(R1-L*round(R1/L), 2))/(4*lamb*length*dt);
                            }
                            else
                            {
                                exponent-=(pow(R[p][startBead][ii]-R[pp][endBead][ii], 2))/(4*lamb*length*dt);
                            }
                        }
                    }
                    particleTable[pp]*=exp(exponent);
                    
                }
                
                else
                {
                    particleTable[pp]=0;
                }
                probability[pp]=particleTable[pp]+lastProbability;
                lastProbability=probability[pp];
            
                Q+=particleTable[pp];
                
            }
            
            //Acceptence test for permutation
            double testExponent=0;
            
            for (int ppp=0; ppp<numberOfParticlesUsed; ppp++)
            {
                for (int ii=0; ii<dim; ii++)
                {
                    if (startBead>(M-length))
                    {
                        if (periodic)
                        {
                            double R1=R[p2][startBead][ii]-R[permutation[p2]][endBead][ii];
                            testExponent-=(pow(R1-L*round(R1/L), 2))/(4*lamb*length*dt);
                        }
                        
                        testExponent-=(pow(R[p2][startBead][ii]-R[permutation[p2]][endBead][ii], 2))/(4*lamb*length*dt);
                        
                    }
                    else
                    {
                        if (periodic)
                        {
                            double R1=R[p2][startBead][ii]-R[permutation[p2]][endBead][ii];
                            testExponent-=(pow(R1-L*round(R1/L), 2))/(4*lamb*length*dt);
                        }
                        
                        testExponent-=(pow(R[p2][startBead][ii]-R[p2][endBead][ii], 2))/(4*lamb*length*dt);
                    }
                }

            }
            
            

            
            double denom= Q+pow(4*PI*lamb*length*dt, dim/2)*exp(testExponent);
            double C= Q/denom;
            double p_rand=getRandom();
            
            if (C<p_rand || C!=C)
            {
                return;
            }
            
            for (int pp=0; pp<particles; pp++)
            {
                probability[pp]/=Q;
            }

            
            //Next particle in permutation cycle
            double rand=getRandom();
            for (int pp=0; pp<particles; pp++)
            {
                
                if (rand<probability[pp])
                {
                    particle_goes_to[p2]=pp;
                    permutationTable[p2]=permutation[pp];
                    p2=pp;
                    particlesUsed.push_back(p2);
                    includedInCycle[p2]=true;
                    numberOfParticlesUsed++;
                    break;
              
                }
            }
            
            if (p2==p)
            {
                permute=false;
            }
        
        }
        
    }
    
    
    
    double ***Rtrial = new double **[numberOfParticlesUsed];
    double ***Rtemp = new double **[numberOfParticlesUsed];
    
    for (int pp=0; pp<numberOfParticlesUsed; pp++)
    {
        Rtrial[pp]=new double *[length];
        Rtemp[pp]=new double *[length];
        
    }
    
    for (int pp=0; pp<numberOfParticlesUsed; pp++)
    {
        for (int kk=0; kk<length; kk++)
        {
            Rtrial[pp][kk]=new double [dim];
            Rtemp[pp][kk]=new double [dim];
        }
    }
    //Save the initial configuration
    
    for (int kk=0; kk<length; kk++)
    {
        for (int pp=0; pp<numberOfParticlesUsed; pp++)
        {
            for (int ii=0; ii<dim; ii++)
            {
                int kkk=startBead+kk;
                
                if (kkk>M)
                {
                    int permkkk=kkk-floor(kkk/M)*M;
                    Rtemp[pp][kk][ii]=R[permutation[particlesUsed[pp]]][permkkk][ii];
                    
                }
                
                else
                {
                    Rtemp[pp][kk][ii]=R[particlesUsed[pp]][kkk][ii];
                }
            }
        }
        included[kk]=false;
        hasV[kk]=false;
    }
    
    
    //Setup the trial paths for the new bead paths
    
    for (int pp=0; pp<numberOfParticlesUsed; pp++)
    {
        for (int ii=0; ii<dim; ii++)
        {
            Rtrial[pp][0][ii]=R[particlesUsed[pp]][startBead][ii];
            
            if (startBead+length>M)
            {
                int endBead=startBead+length-1-floor((startBead+length-1)/M)*M;
                Rtrial[pp][length-1][ii]=R[permutationTable[particlesUsed[pp]]][(endBead)][ii];
            }
            else
            {
                Rtrial[pp][length-1][ii]=R[particle_goes_to[particlesUsed[pp]]][(startBead+length-1)][ii];
            }
            
        }
    }
    
    
    
    
    included[0]=included[length-1]=hasV[0]=hasV[length-1]=true;
    
    double dAction=0;
    dT=0;
    dV=0;
    while(run)
    {
        //Move the trial bead(s)
        /*
         FIX
                THIS
                        PERIODIC
                                BOUNDARY
         CONDITION
         
         */
        for(int kk=0; kk<length; kk+=stepsize)
        {
            if (included[kk]==false)
            {
                
                for (int pp=0; pp<numberOfParticlesUsed; pp++)
                {
                    
                    for (int ii=0; ii<dim; ii++)
                    {
                        double eta=boxMuller()*sqrt(lamb*stepsize*dt);
                        if (periodic)
                        {
                            
                            double Rdist=Rtrial[pp][kk+stepsize][ii]-Rtrial[pp][kk-stepsize][ii];
                            if (abs(Rdist)>L/2)
                            {
                                Rtrial[pp][kk][ii]=0.5*(Rtrial[pp][kk+stepsize][ii]+Rtrial[pp][kk-stepsize][ii])+0.5*L+eta;
                            }
                            else
                            {
                                Rtrial[pp][kk][ii]=(Rtrial[pp][kk+stepsize][ii]+Rtrial[pp][kk-stepsize][ii])/2+eta;
                            }
                            if (abs(Rtrial[pp][kk][ii])>(L/2))
                            {
                                Rtrial[pp][kk][ii]=pbc(Rtrial[pp][kk][ii]);
                            }
                            
                            
                        }
                        else
                        {
                            Rtrial[pp][kk][ii]=(Rtrial[pp][kk-stepsize][ii]+Rtrial[pp][kk+stepsize][ii])/2+eta;
                        }
                    }
                    
                }
                included[kk]=true;
            }
            
            
        }
        //Calculate action cost of move
        
        for (int kk=0; kk<(length); kk+=stepsize)
        {
            int currentBead=startBead+kk;
            
            if (included[kk]==true && hasV[kk]==false)
            {
                
                for (int ppp=0; ppp<particles; ppp++)
                {
                    
                    if (includedInCycle[ppp]==false)
                    {
                        if (currentBead<M)
                        {
                            for (int pp=0; pp<numberOfParticlesUsed; pp++)
                            {
                                
                                double r_diff[dim];
                                double r_trial_diff[dim];
                                double rDiff;
                                double rTrialDiff;
                                    
                                for (int ii=0; ii<dim; ii++)
                                {
                                    r_diff[ii] = Rtemp[pp][kk][ii]-R[ppp][currentBead][ii];
                                    r_trial_diff[ii] = Rtrial[pp][kk][ii]-R[ppp][currentBead][ii];
                                }
                                    
                                rDiff=norm(r_diff[0], r_diff[1]);
                                rTrialDiff=norm(r_trial_diff[0], r_trial_diff[1]);
                                    
                                int i_approx = floor(rDiff*ss/L);
                                int i_plus = i_approx+1;
                                int j_approx = floor(rTrialDiff*ss/L);
                                int j_plus=j_approx+1;
                                dV+=((Pot[j_approx]+Pot[j_plus])-(Pot[i_approx]+Pot[i_plus]))/(2*M);
                            }

                        }
                        
                        else
                        {
                            int bead=currentBead-round(currentBead/M)*M;
                                
                            for (int pp=0; pp<numberOfParticlesUsed; pp++)
                            {
                                    
                                double r_diff[dim];
                                double r_trial_diff[dim];
                                double rDiff;
                                double rTrialDiff;
                                    
                                for (int ii=0; ii<dim; ii++)
                                {
                                    r_diff[ii] = Rtemp[pp][kk][ii]-R[permutation[ppp]][bead][ii];
                                    r_trial_diff[ii] = Rtrial[pp][kk][ii]-R[permutation[ppp]][bead][ii];
                                }
                                    
                                rDiff=norm(r_diff[0], r_diff[1]);
                                rTrialDiff=norm(r_trial_diff[0], r_trial_diff[1]);
                                    
                                int i_approx = floor(rDiff*ss/L);
                                int i_plus = i_approx+1;
                                int j_approx = floor(rTrialDiff*ss/L);
                                int j_plus=j_approx+1;
                                dV+=((Pot[j_approx]+Pot[j_plus])-(Pot[i_approx]+Pot[i_plus]))/(2*M);
                            }

                        }
                            
                    }
                    
                }
                
                for (int pp=0; pp<numberOfParticlesUsed; pp++)
                {
                    
                    for (int ppp=pp; ppp<numberOfParticlesUsed; ppp++)
                    {
                        if (particlesUsed[ppp]!=particlesUsed[pp])
                        {
        
                            double r_diff[dim];
                            double r_trial_diff[dim];
                            double rDiff;
                            double rTrialDiff;
                            
                            for (int ii=0; ii<dim; ii++)
                            {
                                r_diff[ii] = Rtemp[pp][kk][ii]-Rtemp[ppp][kk][ii];
                                r_trial_diff[ii] = Rtrial[pp][kk][ii]-Rtrial[ppp][kk][ii];
                            }
                                
                            rDiff=norm(r_diff[0], r_diff[1]);
                            rTrialDiff=norm(r_trial_diff[0], r_trial_diff[1]);
                                
                            int i_approx = floor(rDiff*ss/L);
                            int i_plus = i_approx+1;
                            int j_approx = floor(rTrialDiff*ss/L);
                            int j_plus=j_approx+1;
                            dV+=((Pot[j_approx]+Pot[j_plus])-(Pot[i_approx]+Pot[i_plus]))/(2*M);
                        }
                    }
                }
            
            }
        
        
            for (int pp=0; pp<numberOfParticlesUsed; pp++)
            {
                dV+=externPot(norm(Rtrial[pp][kk][0],Rtrial[pp][kk][1]));
            }
                
            hasV[kk]=true;
        }
            
        
            //Test the action change
        double ran= getRandom();
        double dA=-B*stepsize*(dV)+dAction;
        double A=exp(dA);
        //  cout<<"\ndV="<<dV<<"\nA="<<A<<"\nran="<<ran<<endl;
        
        
        if (A<ran)
        {
            run=false;
        }
        
        else
        {
            //cout<<stepsize<<endl;
            if (stepsize>1)
            {
                dAction = B*stepsize*dV;
                stepsize/=2;
            }
            else
            {
                for (int pp=0; pp<numberOfParticlesUsed; pp++)
                {
                    for(int kk=1; kk<length; kk++)
                    {
                        for (int ii=0; ii<dim; ii++)
                        {
                            if (periodic)
                            {
                                double Rnew = (Rtrial[pp][kk][ii]-Rtrial[pp][kk-1][ii]);
                                double Rold = (Rtemp[pp][kk][ii]-Rtemp[pp][kk-1][ii]);
                                
                                dT+=M*(pow(pbc(Rnew), 2)
                                       -pow(pbc(Rold),2))/(4*lamb*B*B);
                            }
                            else
                            {
                                dT+=M*(pow((Rtrial[pp][kk][ii]-Rtrial[pp][kk-1][ii]), 2)
                                   -pow((Rtemp[pp][kk][ii]-Rtemp[pp][kk-1][ii]), 2))/(4*lamb*B*B);
                            }
                        }
                    }
                }
                
                
                E+=-dT+dV;
                S+=2*dT/B;
                T+=dT;
                V+=dV;
                if (E!=E)
                {
                    cout<<T<<"\n"<<V<<"\n"<<endl;

                    
                    cout<<"-----------------"<<endl;
                    for (int kk=0; kk<length; kk++)
                    {
                        cout<<"R[0]["<<startBead+kk<<"]=("<<Rtemp[0][kk][0]<<", "<<Rtemp[0][kk][1]<<" )---->("<<Rtrial[0][kk][0]<<", "<<Rtrial[0][kk][1]<<")"
                        <<"\tR[0]["<<startBead+kk<<"]=("<<Rtemp[1][kk][0]<<", "<<Rtemp[1][kk][1]<<" )---->("<<Rtrial[1][kk][0]<<", "<<Rtrial[1][kk][1]<<")"<<endl;
                    }
                    cout<<"-----------------"<<endl;
                    for(int pp=0; pp<particles; pp++)
                    {
                        cout<<"Worldine "<<pp<<"---->"<<permutation[pp]<<endl;
                    }
                }
                //cout<<T<<"\n"<<V<<"\n"<<endl;
                swapCount++;
                
                
                
                //swap all follwowing worldline coordinates with each other
                
                for (int pp=0; pp<numberOfParticlesUsed; pp++)
                {
                    permutation[particlesUsed[pp]]=permutationTable[particlesUsed[pp]];
                }
                
                
                for (int ii=0; ii<dim; ii++)
                {
                    for (int kk=0; kk<length; kk++)
                    {
                        for (int pp=0; pp<numberOfParticlesUsed; pp++)
                        {
                            int currentBead=startBead+kk;
                            if (currentBead>(M-1))
                            {
                                int perCurrentBead=currentBead-floor(currentBead/M)*M;
                                R[permutation[particlesUsed[pp]]][perCurrentBead][ii]=Rtrial[pp][kk][ii];
                            }
                            
                            else
                            {
                                R[particlesUsed[pp]][currentBead][ii]=Rtrial[pp][kk][ii];
                            }
                            
                        }
                    }
                    
                        
                        
                    double Rquicksave;
                    for (int kk=startBead+length; kk<M; kk++)
                    {
                            
                        Rquicksave=R[particlesUsed[0]][kk][ii];
                        for (int pp=0; pp<numberOfParticlesUsed-1; pp++)
                        {
                            R[particlesUsed[pp]][kk][ii]=R[particle_goes_to[particlesUsed[pp]]][kk][ii];
                        }
                        R[particlesUsed[numberOfParticlesUsed-1]][kk][ii]=Rquicksave;
                            
                    }
                        
                        /*Rquicksave=R[permutation[p]][M][ii];
                         R[permutation[p]][M][ii]=R[permutation[p2]][M][ii];
                         R[permutation[p2]][M][ii]=Rquicksave;//*/
                        
                    
                }
                
                for (int ii=0; ii<dim; ii++)
                {
                    for (int pp=0; pp<numberOfParticlesUsed; pp++)
                    {
                        R[particlesUsed[pp]][M][ii]=R[permutation[particlesUsed[pp]]][0][ii];
                    }
                    
                }
                
                run=false;
                
            }
        }
    }
    
    
    for (int pp=0; pp<numberOfParticlesUsed; pp++)
    {
        for (int kk=0; kk<length; kk++)
        {
            delete[] Rtrial[pp][kk];
            delete[] Rtemp[pp][kk];
            
        }
    }
    
    for (int pp=0; pp<numberOfParticlesUsed; pp++)
    {
        delete[] Rtrial[pp];
        delete[] Rtemp[pp];
        
    }
    
    delete[] Rtrial;
    delete[] Rtemp;
    return;
}









