#include "Cluster.h"
#include "Debug.h"
#include "TreeData.h"
#include<set>

// static variable of Clz struct
unsigned int Clz::Count = 0;
unsigned int Clz::Mult = 0;
unsigned int Clz::Height = 0;
unsigned int Clz::MatrixHeight = 0;

Cluster::Cluster():
AnalysisObject(),
hSingle_pxl_cluster(nullptr), hMultiplicity(nullptr),  hPixelheight_map(nullptr),
cSingle_pxl_cluster(nullptr), cMultiplicity(nullptr), cPixelheight(nullptr),
hMatrix_sector{nullptr}, hCluster_sector{nullptr},  
hMultiplicity_sector{nullptr},  
// hRow_pxl_sector{nullptr}, hCol_pxl_sector{nullptr},
cMatrix_sector(nullptr), cCluster_sector(nullptr),  cMultiplicity_sector(nullptr)
// ,cRow_pxl_sector(nullptr), cCol_pxl_sector(nullptr)
{
    
}

Cluster::~Cluster(){

}

void Cluster::Book(){

    // single pixel clusters
    Save( hSingle_pxl_cluster = new TH1F("hSingle_pxl_cluster", "Single pxl cluster distr", 100, 0, 3700) );
    hSingle_pxl_cluster->GetXaxis()->SetTitle("Signal [ADC]");
    Save( CreateCanvas(cSingle_pxl_cluster, hSingle_pxl_cluster) );

    // cluster multiplicity
    Save( hMultiplicity = new TH1F("hMultiplicity", "Cluster multiplicity", m_Nrows, 0, m_Nrows) );
    hMultiplicity->GetXaxis()->SetTitle("Cluster size [# of pixels]");
    Save( CreateCanvas(cMultiplicity, hMultiplicity) );

    Save( hPixelheight_map = new TH2F("hPixelheight_map", "Pixel height map",5, -0.5, 5- 0.5, 5, - 0.5,5 - 0.5) );
    hPixelheight_map->GetXaxis()->SetTitle("Column");
    hPixelheight_map->GetYaxis()->SetTitle("Row");
    Save( CreateCanvas(cPixelheight, hPixelheight_map) );


  
    // creation of histograms for sectors
    for(unsigned int i=0; i < m_Nsectors; ++i){

        Save( hCluster_sector[i] = new TH1F ( ("hCluster_sector_" + std::to_string(i)).c_str(),
            ("Cluster signal sector # " + std::to_string(i)).c_str(), 100, 0, 25000) );   
        hCluster_sector[i]->GetXaxis()->SetTitle("Signal [ADC]");

        Save( hMatrix_sector[i] = new TH1F ( ("hMatrix_sector_" + std::to_string(i)).c_str(),
            ("Matrix sector # " + std::to_string(i)).c_str(), 100, 0, 25000) );
        hMatrix_sector[i]->GetXaxis()->SetTitle("Signal [ADC]");

        Save( hMultiplicity_sector[i] = new TH1F ( ("hMultiplicity_sector_" + std::to_string(i)).c_str(),
            ("Multiplicity sector # " + std::to_string(i)).c_str(), m_Nrows, 0, m_Nrows) );
        hMultiplicity_sector[i]->GetXaxis()->SetTitle("Number of pixels");


        //Save( hRow_pxl_sector[i] = new TH1F ( ("hRow_pxl_sector_" + std::to_string(i)).c_str(),
        //    ("Cluster size along rows sector # " + std::to_string(i)).c_str(), m_MatrixDim, 0, m_MatrixDim) );
        //hRow_pxl_sector[i]->GetXaxis()->SetTitle("Pixel row");

        //Save( hCol_pxl_sector[i] = new TH1F ( ("hCol_pxl_sector_" + std::to_string(i)).c_str(),
        //    ("Cluster size along cols sector # " + std::to_string(i)).c_str(), m_MatrixDim, 0, m_MatrixDim) );
        //hCol_pxl_sector[i]->GetXaxis()->SetTitle("Pixel column");

    }

    Save( CreateCanvasForSector(cCluster_sector, hCluster_sector) );
    Save( CreateCanvasForSector(cMatrix_sector, hMatrix_sector) );
    Save( CreateCanvasForSector(cMultiplicity_sector, hMultiplicity_sector) );

    return;
}

void Cluster::Begin(const CommandLine& cl){

    std::cout << "Starting cluster analysis ... \n";
    std::cout << std::endl;

    // default value for pixel threshold = 4
    m_Pixel_thr = (cl.contains("pixelThr")) ? (cl.getValue("pixelThr")) : 4;
    DBPRINT("Pixel threshold set to : ", m_Pixel_thr);
    Cluster::Book();

    for(int i = 0; i<Clz::MatrixSize; ++i)
        pixelHeight.push_back(0);
    
    return;
}

void Cluster::Analyze(unsigned int iClz, unsigned int iSec){

    TreeData<int>* dataInt = TreeData<int>::GetInstance();
    TreeData<float>* dataFloat = TreeData<float>::GetInstance();

    //Clz::Height = 0;

    int row_seed = dataInt->Get("row_seed")->at(iClz);
    int col_seed = dataInt->Get("col_seed")->at(iClz);
    
    // loop over the matrix (sub-matrix 5x5)
    for(unsigned int k=0; k < Clz::MatrixSize; ++k){
        // set threshold on neighbor pixels
        unsigned int neighbourIndex = Clz::MatrixSize*iClz + k;
        float pulseHeight = dataFloat->Get("ph_nxn")->at(neighbourIndex);
        float noise = dataFloat->Get("noise_nxn")->at(neighbourIndex);
        if( pulseHeight/noise > m_Pixel_thr ){
            // update cluster mass
            Clz::Height += pulseHeight;
            ++Clz::Mult;
        }

            if(iSec==2){//&& row_seed==11 && col_seed==16){
            pixelHeight[k] += pulseHeight;
            hPixelheight_map->SetBinContent(hPixelheight_map->FindFixBin(k/5, k%5), pixelHeight[k]);
        }


     
        //std::cout << "Clz height = " << Clz::Height << std::endl;

        // update matrix mass
        Clz::MatrixHeight += pulseHeight;
    }
    /*
    // save seed coordinates
    if(iSec==2){
        std::pair<int,int> p = std::make_pair(row_seed, col_seed);
        seed_coordinate.push_back(p);

        std::string title_s = "Pixel height map Seed = (" + std::to_string(p.first) + "," + std::to_string(p.second) + std::to_string(iClz);
        const char* title = title_s.c_str();
        std::string name_s = "hPixelheight_map" + std::to_string(p.first) + std::to_string(p.second) + std::to_string(iClz);
        const char* name = name_s.c_str();
        hPixelheight_map->SetTitle(title);
        hPixelheight_map->Add(Cluster::CreatePixelHeight(iClz, p.first, p.second, 2));
                
                //pixelHeight_seed->Add(CreatePixelHeight(cluster_seed_coordinate[i], seed_coordinate[i].first, seed_coordinate[i].second, 2));
    }*/
       
        
    // seed height + all overthreshold pixel height
    hCluster_sector[iSec]->Fill(Clz::Height); 
    
    // sum of all pixel heights in the sub-matrix
    hMatrix_sector[iSec]->Fill(Clz::MatrixHeight);

    // multiplicity distribution: how many pixels are over threshold
    hMultiplicity_sector[iSec]->Fill(Clz::Mult);
    hMultiplicity->Fill(Clz::Mult);

    // cluster with single pixel, i.e. with seed only
    if(Clz::Mult == 1) hSingle_pxl_cluster->Fill(Clz::Height);

    return;
}


TH2F* Cluster::CreatePixelHeight(unsigned int iClz, unsigned int row, unsigned int col, unsigned int iSec){

    TreeData<float>* dataFloat = TreeData<float>::GetInstance();
    std::string name_s = "hPixelheight_map"+(std::to_string(row)+(std::to_string(col)));
    const char* name = name_s.c_str();
    TH2F* hPixelheight_map = new TH2F(name, "Pixel height map", 5, -0.5, 4.5, 5, -0.5, 4.5);
    for(int x=0; x<sqrt(Clz::MatrixSize); ++x){
        for(int y=0; y<sqrt(Clz::MatrixSize); ++y){
            unsigned int k= x+y*5;
            unsigned int neighbourIndex = Clz::MatrixSize*iClz + k;
            float pulseHeight = dataFloat->Get("ph_nxn")->at(neighbourIndex);
            hPixelheight_map->SetBinContent(x,y, hPixelheight_map->GetBinContent(x,y)+ pulseHeight);
        }
    }
    return hPixelheight_map;
}

float Cluster::CalculateCoG(TH2F* hPixel_map){
    float cog_x = 0;
    float cog_y = 0;
    float sum = 0;

    // center of gravity
    for(int x = 0; x < sqrt(Clz::MatrixSize); ++x){
        for(int y = 0; y < sqrt(Clz::MatrixSize); ++y){  
            cog_y = cog_y + y*hPixel_map->GetBinContent(x,y);
            cog_x = cog_x + x*hPixel_map->GetBinContent(x,y);
            sum = sum + hPixel_map->GetBinContent(x,y);
        }
    }

    return (cog_x/sum, cog_y/sum);
}



void Cluster::End(){

    hPixelheight_map->Scale(1/hPixelheight_map->GetBinContent(hPixelheight_map->GetMaximumBin()));
    //hPixelheight_map->Scale(-1);
    //hPixelheight_map->GetZaxis()->SetRangeUser(5e+05, hPixelheight_map->GetBinContent(hPixelheight_map->GetMaximum()));
    std::cout << "Cluster analysis completed. Plotting the results ... \n";
    m_Already_analyzed = true;

    if(m_PlotAll){
        for(auto it = m_PlotPairs.begin(); it != m_PlotPairs.end(); ++it)
            Cluster::Plot(it->first, it->second);
        for(auto it = m_PlotPairsSectors.begin(); it != m_PlotPairsSectors.end(); ++it)
            Cluster::PlotForSector(it->first, it->second);
    }

    else if(m_PlotSec){
        for(auto it = m_PlotPairsSectors.begin(); it != m_PlotPairsSectors.end(); ++it)
            Cluster::PlotForSector(it->first, it->second);
    }

    else{
        for(auto it = m_PlotPairs.begin(); it != m_PlotPairs.end(); ++it)
            Cluster::Plot(it->first, it->second);
    }

    return;
}

