#include "master.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "mpi.h"
//#include "test.h"

using namespace std;


int main(int argc, char* argv[])
{
        //test();
        

        //master m(0.99,1,0.005, 0.005, 0.005);
        //cout<<abs(m.x()(0))<<endl;

        MPI_Status status;

        int myid, numprocess;

        int groupn(40);
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocess);

        int everygroup(groupn/numprocess);

        if(myid==0)
        {
                ofstream outfile("result");
                std::vector<double> res(groupn);
                
                for(int i=0; i<everygroup; ++i)
                {
                        double delta(0.000001*(everygroup*myid+i)+0.96497);
                        master M(delta, delta, 0, 0.005, 0.005, 0.005);
                        res.at(i)=abs(M.x()(0));
                        


                }

                

                for(int id=1; id<numprocess; ++id)
                {
                        //std::vector<double> res(everygroup);

                        
                        MPI_Recv(&res[everygroup*id], everygroup, MPI_DOUBLE, id, id, MPI_COMM_WORLD, &status);
                        
                }
                
                for(int i=0; i<groupn; ++i)
                {
                        outfile<<0.96497+0.000001*i<<"\t"<<res.at(i)<<endl;
                }

                
                outfile.close();
                

        }else
        {
                std::vector<double> res(everygroup);
                for(int i=0; i<everygroup; ++i)
                {
                        double delta(0.000001*(everygroup*myid+i)+0.96497);
                        master M(delta, delta, 0, 0.005, 0.005, 0.005);
                        res.at(i)=abs(M.x()(0));

                }

                MPI_Send(&res[0], everygroup, MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);

        }



        MPI_Finalize();
}