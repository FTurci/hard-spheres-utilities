#include "histogram.h"


histogram::histogram(){

}

void histogram::setup(double minv, double maxv, double binwidth){
    // cerr<<"** Setting the histogram: binsize = "<< binwidth<<endl;
    min=minv;
    max=maxv;
    inv_binWidth=1./binwidth;
    binCount=(max-min)/binwidth;
    counts.resize(binCount);
    // cerr<<"** Setting the histogram: bincount = "<< binCount<<endl;
}
void histogram::record(double datum) {
    int bin = (int)((datum - min) *inv_binWidth);

    if (bin < 0) {
        lowerOutlierCount++;
    } else if (bin >= binCount) {
        upperOutlierCount++;
    } else {
    
        counts[bin]++;
    }
}