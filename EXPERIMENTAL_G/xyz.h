#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "histogram.h"
using namespace std;

class xyz
{

public:
    xyz(string filename, double binw, double micro_pixel_xy,double micro_pixel_z);
    void read_conf(string filename);
    void print_conf();

    void remove_borders(double a, double b, double c);
    void compute_ideal_g();
    void compute_real_g();
    void write_g_to(string file);


    std::vector<double> original_x,original_y,original_z;
    std::vector<double> x,y,z;
    std::vector<char> original_type;

    std::vector<double> Ideal;

    int N;
    double max_x, max_y, max_z,min_x, min_y, min_z;
    double binwidth;

    histogram ideal_hist,real_hist;

    std::vector<double> g;

};