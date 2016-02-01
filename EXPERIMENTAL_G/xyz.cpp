#include "xyz.h"


double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double distance(double x1,double y1,double z1,double x2,double y2,double z2){
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}

xyz::xyz(string filename, double binw, double micron_pixel_xy,double micron_pixel_z){

    binwidth=binw;
    ifstream fin(filename.c_str(), ios::in);
    fin>>N;
    this->original_type.resize(N); this->original_x.resize(N);this->original_y.resize(N);this->original_z.resize(N);
    string line;

    max_x=0, max_y=0,max_z=0;
    
    fin>>line; getline(fin, line); cout<<line; //one word

    for (int i = 0; i < this->N; ++i)
    {
        // cout<<"reading"<<endl;
        fin>>this->original_type[i]>>this->original_x[i]>>this->original_y[i]>>this->original_z[i];

        original_x[i]*=micron_pixel_xy;
        original_y[i]*=micron_pixel_xy;
        original_z[i]*=micron_pixel_z;
        if(original_x[i]>max_x) max_x=original_x[i];
        if(original_y[i]>max_y) max_y=original_y[i];
        if(original_z[i]>max_z) max_z=original_z[i];
    }
    cout<<"Box edges "<<max_x<<" "<<max_y<<" "<<max_z<<endl;
    

fin.close();
}

void xyz::print_conf(){
    
    cout<<"====> The number of particles is "<<N<<endl;
    for (int i = 0; i < x.size(); ++i)
    {
        std::cout<<x[i]<<" "<<y[i]<<" "<<z[i]<<endl;
    }
}

void xyz::remove_borders(double a, double b, double c){

    int removed=0;
    min_x=a, min_y=b,min_z=c;
    for (int i = 0; i < N; ++i)
    {
        if(original_x[i]<a || original_y[i]<b || original_z[i]<c || original_x[i]>max_x-a || original_y[i]>max_y-b || original_z[i]>max_z-c)
        {
            removed++;
        }

        x.push_back(original_x[i]);
        y.push_back(original_y[i]);
        z.push_back(original_z[i]);
    }
    N-=removed;
    
    cout<<"====> Removed "<<removed<<" particles. The remaining particles are "<<N<<endl;

    cout<<"====> The number density is "<< N/((max_x-min_x)*(max_y-min_y)*(max_z-min_z))<<" particle/micron^3."<<endl;
    ideal_hist.setup(0,max_x*0.5, binwidth);
    real_hist.setup(0,max_x*0.5, binwidth);
}

void xyz::compute_ideal_g(){

    std::vector<double> xi(N),yi(N),zi(N);

    for (int i = 0; i < N; ++i)
    {
        xi[i]=fRand(min_x, max_x);
        yi[i]=fRand(min_y, max_y);
        zi[i]=fRand(min_z, max_z);
    }

    double d;

    // compute the distances
    for (int i = 0; i < N-1; ++i)
    {
        cerr<<i<<'\r';
        for (int j = i+1; j < N; ++j)
        {
            d=distance(xi[i],yi[i],zi[i],xi[j],yi[j],zi[j]);
            ideal_hist.record(d);
        }
    }

}

void xyz::compute_real_g(){

    double d;

    // compute the distances
    for (int i = 0; i < N-1; ++i)
    {
        cerr<<i<<'\r';
        for (int j = i+1; j < N; ++j)
        {
            d=distance(x[i],y[i],z[i],x[j],y[j],z[j]);
            real_hist.record(d);
        }
    }

    g.resize(real_hist.counts.size());
    for (int i = 0; i < g.size(); ++i)
    {
        if(ideal_hist.counts[i]!=0)
            g[i]=real_hist.counts[i]/(double)ideal_hist.counts[i];
    }

}

void xyz::write_g_to(string filename){
    ofstream fout(filename.c_str(), ios::out);
    fout<<"#r g"<<endl;
    for (int i = 0; i < g.size(); ++i)
    {
        fout<<binwidth*(i+1)-0.5*binwidth<<" "<<g[i]<<" "<<real_hist.counts[i]<<" "<<ideal_hist.counts[i]<<endl;
    }

}