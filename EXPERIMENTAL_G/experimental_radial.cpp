#include <iostream>
#include "xyz.h"

using namespace std;

int main(int argc, char const *argv[])
{
    
    double binsize=atof(argv[3]);
    double borderx=atof(argv[4]);
    double bordery=atof(argv[5]);
    double borderz=atof(argv[6]);

    double micron_pixel_xy=atof(argv[7]);
    double micron_pixel_z=atof(argv[8]);

    cout<<"* Reading file "<<argv[1]<<endl;
    xyz File(argv[1],binsize, micron_pixel_xy, micron_pixel_z) ;

    File.print_conf();
    cout<<"* Removing borders... "<<endl;
    File.remove_borders(borderx, bordery, borderz);

    cout<<"* Computing the ideal gas.."<<endl;
    File.compute_ideal_g();
    cout<<"* Computing the real distribution.."<<endl;
    File.compute_real_g();
    cout<<"* Writing the g(r) to "<<argv[2]<<endl;
    File.write_g_to(argv[2]);
    return 0;
}