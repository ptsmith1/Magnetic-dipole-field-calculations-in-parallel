#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

double avg_1D(std::vector<double> &array)
{
    double sum=0;
    int count=0;
    for (int i =0; i<array.size();i++)
    {
        if (array[i]!=0)
        {
            sum = sum + array[i];
        }
        else
        {
            count++;
        }
    }
    return(sum/(array.size()-count));
}

//each call of this function returns the x,y,z component of the total B
//field for a single point within the ellipsoid

double* getB(double m[],int nx,int nxny,int i,double k, double l,int N,int x,int y,int z)
    {
        int j;//tracks the array position of the dipole acting on dipole i
        int x2,y2,z2; //coordinates of acting dipole
        int dx,dy,dz;
        double mag,rdotm,r_mag_cubed,k_div_r_mag_cubed;
        double r_hat[3] ={0};//stores unit vector from j to i
        double* out = new double[3];//stores cumulative bx,by,bz values acting on i
        z= i/(nxny);//calculates x,y,z coordinates starting from (0,0,0) of each dipole i
        y=(i-nxny*z)/nx;
        x=i%nx;
        out[0]=0;
        out[1]=0;
        out[2]=0;
        for (j=0;j<N;j++)
        {
            if (i!=j)//to filter out particles from acting on themselves
            {
                z2= j/(nxny);//x,y,z coordinates of each dipole j
                y2=(j-nxny*z2)/nx;
                x2=j%nx;
                dx=x2-x;
                dy=y2-y;
                dz=z2-z;
                mag = sqrt((dx)*(dx)+(dy)*(dy)+(dz)*(dz)); //magnitude of r
                r_mag_cubed = (mag*l)*(mag*l)*(mag*l);
                k_div_r_mag_cubed=k/r_mag_cubed;
                r_hat[0] = (dx)/mag;//calculate and normalise r
                r_hat[1] = (dy)/mag;
                r_hat[2] = (dz)/mag;
                rdotm = (r_hat[0]*m[0]+r_hat[1]*m[1]+r_hat[2]*m[2]);
                //calculate contribution to B at i due to j
                out[0] = out[0]+k_div_r_mag_cubed*((3*r_hat[0]*rdotm-m[0]));
                out[1] = out[1]+k_div_r_mag_cubed*((3*r_hat[1]*rdotm-m[1]));
                out[2] = out[2]+k_div_r_mag_cubed*((3*r_hat[2]*rdotm-m[2]));
            }
        }
        return(out);
    }

bool inside_ellipsoid(float x,float y,float z,float rx,float ry,float rz)
{
    //returns true/false depending on if a particles coordinates are within the ellipsoid
    bool in = (((x-rx)*(x-rx))/(rx*rx)+((y-ry)*(y-ry))/(ry*ry)+((z-rz)*(z-rz))/(rz*rz))<=1.00;
    return in;
}

void print_xyplane(std::vector<double> const &Bx,int nx,int ny,int layer)
{
    //for a chosen z layer, prints a x-y grid of Bx/By/Bz values
    for (int i =0;i<ny;i++)
    {
        //the outer loop iterates over each row, the inner loop iterates from the start of that
        //row on the specified layer and the end of that row,ie over each column
        for (int j =layer*nx*ny+i*nx;j<layer*nx*ny+i*nx+nx;j++)
        {

                std::cout<<Bx[j]<<" ";
        }
        std::cout<<std::endl;
    }
}

void print_zyplane(std::vector<double> const &Bx,int nz,int ny,int nx,int layer)
{
    //for a chosen z layer, prints a z-y grid of Bx/By/Bz values
    for (int i =0;i<ny;i++)
    {
        //takes significantly longer to print due to non contiguous memory access
        for (int j =layer+i*ny;j<layer+nx*ny*nz+i*ny;j+=nx*ny)
        {

                std::cout<<Bx[j]<<" ";
        }
        std::cout<<std::endl;
    }
}

    
int main (int argc, char *argv[])
{
    auto t1 = Clock::now(); //timer
///////////////////////////////Simulation user input variables//////////////////////////////////////
    
    double u = 9.274009e-24; //magnetic dipole moment=bohr magneton
    float rx= 10e-9; //x,y,z radius of ellipsoid
    float ry= 10e-9;
    float rz= 20e-9;
    const double l=3e-10; //lattice seperation
    double m[3] = {1,0,0}; //magnetic moment direction
    
///////////////////////////////////////Simulation setup///////////////////////////////////////////// 
  
    double k = u * 1e-7;
    int nx=2*rx/l; //lattice dimensions
    int ny=2*ry/l;
    int nz=2*rz/l;
    int N=nx*ny*nz;
    int nxny = nx*ny;
    int z,y,x;
    
///////////////////////////////////////////Simulation////////////////////////////////////////
    
    double B_self[3]={0,0,0};//stores the self generated magnetic field in each axis
    std::vector<double> Bx(N,0);//stores x,y,z component of the magnetic field respectively
    std::vector<double> By(N,0);
    std::vector<double> Bz(N,0);
    double* ptr;//A temporary pointer which stores the output of the result of getB
    
    B_self[0] = (k*8*M_PI)/(3*(l*l*l))*m[0];//calculate B due to self
    B_self[1] = (k*8*M_PI)/(3*(l*l*l))*m[1];
    B_self[2] = (k*8*M_PI)/(3*(l*l*l))*m[2];
    
    //The following structure loops through all dipoles, for each dipole within the ellipsoid it calls the getB function which calculates the magnetic field due to the interaction of all N dipoles. It then combines the self generated and interaction magnetic fields
    for (int i=0;i<N;i++)
    {
        z= i/(nx*ny);
        y=(i-nx*ny*z)/nx;
        x=i%nx;
        if (inside_ellipsoid(x,y,z,((float)(nx)-1)/2,((float)(ny)-1)/2,((float)(nz)-1)/2) ==1)//only calculates B for dipoles within ellipsoid
        {
            ptr=getB(m,nx,nxny,i,k,l,N,x,y,z);
            Bx[i]=ptr[0] + B_self[0];//combines interacting and self generated fields
            By[i]=ptr[1] + B_self[1];
            Bz[i]=ptr[2] + B_self[2];
        }
    }
    
////////////////////////////////Final parameter calculation and data output/////////////////////////
   
    double D[3] ={0,0,0}; //demagnetizing factors
    double B_avg[3]={0,0,0}; //Average magnetic fields
    //find average B in each axis
    B_avg[0] = avg_1D(Bx);
    B_avg[1] = avg_1D(By);
    B_avg[2] = avg_1D(Bz); 
    //calculate demganetizing factor in each axis
    D[0] = 1-(B_avg[0]*(l*l*l))/(k*4*M_PI);
    D[1] = 1-(B_avg[1]*(l*l*l))/(k*4*M_PI);
    D[2] = 1-(B_avg[2]*(l*l*l))/(k*4*M_PI);
    std::cout<<std::setprecision(30);
    std::cout<<"Average Bx,By,Bz is:("<<B_avg[0]<<","<<B_avg[1]<<","<<B_avg[2]<<")"<<std::endl;
    std::cout<<"Average Dx,Dy,Dz is:("<<D[0]<<","<<D[1]<<","<<D[2]<<")"<<std::endl;
    std::cout<<std::setprecision(6);
    for (int i = 0;i<nz;i++)
    {
        if (i==nz/2)//remove this to print all layers
        {
            print_xyplane(Bx,nx,ny,i);//replacing fBx with fBy/fBz prints corresponding component of the field
        }
    }
    
    for (int i = 0;i<nx;i++)
    {
        //print_zyplane(Bx,nz,ny,nx,i);
    }
    auto t2 = Clock::now();
    double timetaken =(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count())*1e-9;
    std::cout <<"Time taken: "<<timetaken<<"  seconds"<<std::endl;
}

