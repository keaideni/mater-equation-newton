#ifndef MASTER_H
#define MASTER_H
#include <complex>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class master
{
private:
        Matrix3cd _diff;
        Vector3cd _x, _fun;
        complex<double> i;
public:
        const Matrix3cd& diff()const{return _diff;};
        const Vector3cd& x()const{return _x;};
        const Vector3cd& fun()const{return _fun;};


        void calcuf(const double& deltar, const double& deltaq, const double& j,
                const double& f, const double& kappa, const double& gamma)
        {
                _fun(0)=-i*deltar*_x(0)-i*_x(1)-i*f-kappa/2*_x(0)+i*j*_x(0);
                _fun(1)=-i*(deltaq*_x(1)-2.0*_x(0)*_x(2))-gamma/2*_x(1);
                _fun(2)=-i*(_x(0)*conj(_x(1))-conj(_x(0))*_x(1))-gamma*(_x(2)+1.0/2);

        };

        void calcudiff(const double& deltar, const double& deltaq, const double& j,
                const double& f, const double& kappa, const double& gamma)
        {
                _diff(0,0)=-i*deltar-kappa/2+i*j; _diff(0,1)=-i; _diff(0,2)=0;
                _diff(1,0)=i*2.0*_x(2); _diff(1,1)=-i*deltaq-gamma/2; _diff(1,2)=i*2.0*_x(0);
                _diff(2,0)=-i*conj(_x(1)); _diff(2,1)=i*conj(_x(0)); _diff(2,2)=-gamma;
        };


        master(const double& deltar, const double& deltaq, const double& j,
                const double& f, const double& kappa, const double& gamma):
        i(0,1),
        _x(1.0,1.0,1.0)
        {
                Vector3cd tempx(_x);
                double err(1);int ii(0);


                while(err>0.0000001)
                {
                        calcuf(deltar, deltaq, j, f, kappa, gamma);
                        calcudiff(deltar, deltaq, j, f, kappa, gamma);

                        _x-=_diff.inverse()*_fun;

                        err=abs(_x(0)-tempx(0))+abs(_x(1)-tempx(1))+abs(_x(2)-tempx(2));

                        tempx=_x;

                        //cout<<setw(10)<<++ii<<"  =>  "<<setw(20)<<err<<endl;

                }
                


        };
        ~master(){};
        
};



#endif // MASTER_H