#ifndef Cluster_h
#define Cluster_h

#include "AnalysisObject.h"

#include <unordered_map>
#include<set>
#include<vector>
#include<string>

class Cluster : public AnalysisObject{
    
    public:

    Cluster();
    ~Cluster() override;

    void Begin(const CommandLine& cl) override;
    void Analyze(unsigned int iClz, unsigned int iSec) override;
    void End() override;

    TH2F* CreatePixelHeight(unsigned int iClz, unsigned int row, unsigned int col, unsigned int iSec);

    float CalculateCoG(TH2F* hPixelheightmap);

    private:

    const unsigned int m_MatrixDim = 5;
    unsigned int m_Pixel_thr;   // to identify a cluster
    
    void Book() override;


    // single pixel cluster distribution and multiplicity considering all sectors
    TH1F* hSingle_pxl_cluster;
    TH1F* hMultiplicity;
    TH2F* hPixelheight_map;

    // corresponding canvases
    TCanvas* cSingle_pxl_cluster;
    TCanvas* cMultiplicity;
    TCanvas* cPixelheight;

    // histo for sectors
    TH1F* hMatrix_sector[m_Nsectors];
    TH1F* hCluster_sector[m_Nsectors];
    TH1F* hMultiplicity_sector[m_Nsectors];
    //TH1F* hRow_pxl_sector[m_Nsectors];
    //TH1F* hCol_pxl_sector[m_Nsectors];

    // canvases for sector
    TCanvas* cMatrix_sector;
    TCanvas* cCluster_sector;
    TCanvas* cMultiplicity_sector;
    //TCanvas* cRow_pxl_sector;
    //TCanvas* cCol_pxl_sector;


    
    std::vector<std::pair<int,int>> seed_coordinate;
    std::vector<int> cluster_seed_coordinate;
    std::vector<std::pair<int,int>> seed_diff_coord;
    std::set<std::pair<int, int>> seen_seed_coord;


    // map for pixel height
    std::vector<float>pixelHeight;

};

struct Clz{
    static unsigned int Count;
    static unsigned int Mult;
    static unsigned int Height;
    static unsigned int MatrixHeight;
    static const unsigned int MatrixSize = 25;
};

#endif