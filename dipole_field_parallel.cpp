#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <mpi.h>
#include <iomanip>

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
    return(sum/(array.size()-count)); //only average over dipoles inside ellipsoid
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

void errcheck (int ierror, const char* msg) //checks for MPI errors and returns the name of the function causing the error
{
  if (ierror != MPI_SUCCESS)
    {
      printf("Error in %s! \n",msg);
    }
}
    
int main (int argc, char *argv[])
{
////////////////////////////////////////MPI Initiation//////////////////////////////////////////////
    int ierror = 0;
    ierror = MPI_Init(&argc, &argv); errcheck(ierror, "MPI_Init");
    double start_time = MPI_Wtime();
    int rank,num_processors;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
    
///////////////////////////////Simulation user input variables//////////////////////////////////////
    
    double u = 9.274009e-24; //magnetic dipole moment=bohr magneton
    float rx=10e-9; //x,y,z radius of ellipsoid, minimum value for rx,ry,rz is 0.6e-9
    float ry=10e-9;
    float rz =20e-9;
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
    std::vector<int> proc_start_pos((num_processors+1),0); //This stores the dipole that each processor starts working on
    
    //The following structure is responsible for load balancing of the processors by dividing the dipoles between processors so that each has the same amount of dipoles that actualy contribute to the magnetic field.
    //Done in serial due to low computation requirements
    if (rank ==0) 
    {
        int inside_count=0;
        int inside_total=0;
        
        for (int i =0;i<N;i++) //This loop adds up the total dipoles inside the ellopse and stores it as inside_total
        {
            z= i/(nx*ny); //calculates the 3D coordinates of each particle
            y=(i-nx*ny*z)/nx;
            x=i%nx;
            
            if (inside_ellipsoid(x,y,z,((float)(nx)-1)/2,((float)(ny)-1)/2,((float)(nz)-1)/2) ==1)
                inside_count++; //checks if point x,y,z is inside ellipsoid
        }
        
        inside_total=inside_count; 
        inside_count =0;
        int proc_start_pos_count=0;
        //This loop assigns each processor its starting dipole by counting through all points in the lattice and adding 1 to inside_count if the particle is inside the ellipsoid. Every time inside_count reaches a multiple of inside_total/num_processors, the current i is saved to the proc_start_pos array. This results in each processor being assigned the same amount of dipoles inside the ellipsoid but the processors which compute the dipoles near the poles of the ellipsoid will iterate over more dipoles, but most will be outside the ellipsoid and so no calculation is required.
        for (int i =0;i<N;i++) 
        {
            z= i/(nx*ny);
            y=(i-nx*ny*z)/nx;
            x=i%nx;
            
            if (inside_ellipsoid(x,y,z,((float)(nx)-1)/2,((float)(ny)-1)/2,((float)(nz)-1)/2) ==1)
                inside_count++;
            
            if (inside_count==(inside_total/num_processors))
            {
                proc_start_pos[proc_start_pos_count+1]=i;
                proc_start_pos_count++;
                inside_count=0;
            }
        }
        std::cout<<"N is "<<N<<" and nx,ny,nz are "<<nx<<","<<ny<<","<<nz<<std::endl;
        std::cout<<"Total non 0 dipoles "<<inside_total<<std::endl;
        proc_start_pos[num_processors]=N;
    }
    
    //Sends each processor the array of starting positions, they only need a single start position (their own) but it is easier to do it this way and the memory requirement is insignificant. This is a synchronous communication
    ierror = MPI_Bcast(&proc_start_pos[0],num_processors+1,MPI_INT,0,MPI_COMM_WORLD);
    proc_start_pos[num_processors+1]=N; errcheck(ierror, "MPI_Bcast");
    for (int i =0; i<num_processors;i++)
    {
        if(i==rank)
        {
            std::cout<<"My rank is  "<<rank<<" and my starting dipole is "<<proc_start_pos[rank]<<std::endl;
        }
    }
    
///////////////////////////////////////////Simulation////////////////////////////////////////
    
    int N_proc=N/num_processors;
    double B_self[3]={0,0,0};//stores the self generated magnetic field in each axis
    std::vector<double> Bx(N,0);//stores x,y,z component of the magnetic field respectively
    std::vector<double> By(N,0);
    std::vector<double> Bz(N,0);
    double* ptr;//A temporary pointer which stores the output of the result of getB
    
    B_self[0] = (k*8*M_PI)/(3*(l*l*l))*m[0];//calculate B due to self
    B_self[1] = (k*8*M_PI)/(3*(l*l*l))*m[1];
    B_self[2] = (k*8*M_PI)/(3*(l*l*l))*m[2];
    
    //The following structure loops through all dipoles assigned to a processor, for each dipole within the ellipsoid it calls the getB function which calculates the magnetic field due to the interaction of all N dipoles. It then combines the total magnetic field at each dipole
    for (int i=proc_start_pos[rank];i<proc_start_pos[rank+1];i++)
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
    
////////////////////////////////////Vector reduction and deallocation///////////////////////////////
    
    std::vector<double> fBx(N,0); //stores the final magnetic field after reduction from all procs
    std::vector<double> fBy(N,0);
    std::vector<double> fBz(N,0);
    
    //sums the value of B at each position over all processors, for each dipole only the output from 1 processor will be non zero
    ierror = MPI_Reduce(&Bx[0], &fBx[0], N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); errcheck(ierror, "MPI_Reduce");
    MPI_Reduce(&By[0], &fBy[0], N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Bz[0], &fBz[0], N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    //deallocates Bx,By,Bz on all processors
    std::vector<double>().swap(Bx);
    std::vector<double>().swap(By);
    std::vector<double>().swap(Bz);
    
////////////////////////////////Final parameter calculation and data output/////////////////////////
   
    double D[3] ={0,0,0}; //demagnetizing factors
    double B_avg[3]={0,0,0}; //Average magnetic fields
    if (rank==0) 
    {
        //find average B in each axis
        B_avg[0] = avg_1D(fBx);
        B_avg[1] = avg_1D(fBy);
        B_avg[2] = avg_1D(fBz); 
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
                print_xyplane(fBx,nx,ny,i);//replacing fBx with fBy/fBz prints corresponding component of the field
            }
        }
        
        for (int i = 0;i<nx;i++)
        {
            //print_zyplane(fBx,nz,ny,nx,i);
        }
    }

    double end_time = MPI_Wtime(); //timer
    std::cout << "Run time: " << end_time -start_time << " seconds" <<std::endl; 
    ierror = MPI_Finalize(); errcheck(ierror, "MPI_Finalize");
}

