#include <vector>
#include <iostream>

using namespace std;

class histogram {
  public:
    // Pick whichever constructor feels more natural to you
    histogram();

    void setup(double minv, double maxv, double binwidth);
    void record(double datum);

    std::vector<int> counts;
  private:
    double min,max;
    double inv_binWidth;
    int binCount;
    int lowerOutlierCount, upperOutlierCount;
    
};