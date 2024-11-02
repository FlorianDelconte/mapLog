#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/PGMWriter.h"

#include "CLI11.hpp"
#include "UnrolledMap.h"

using namespace DGtal;
typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;


void
scaleCloud(std::vector<Z3i::RealPoint> &pts, double scale){
  for(unsigned int i = 0; i <  pts.size(); i++){
    Z3i::RealPoint &pi =  pts[i];
    pi *= scale;
  }
}

int
main(int argc,char **argv)
{
    int MaxdF {4};
    int minRV {-5};
    int intensityPerUnit {3};
    int pad {200};
    
    std::string meshInput="/volWork/these/source/BERTRAND/exemple_bruit/hetre2.off";
    std::string centerlineInput="/volWork/these/source/BERTRAND/exemple_bruit/hetre2_centerline.xyz";
    std::string cylPointsInput="/volWork/these/source/BERTRAND/exemple_bruit/hetre2-cyl";
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    app.description("Allowed options are: ");
    app.add_option("-i,--input", meshInput , "path to mesh (.off)");
    app.add_option("-l,--centerline", centerlineInput , "path to centerline (.xyz)");
    app.add_option("-c,--cylpoints", cylPointsInput , "path to cylPoints (.xyz)");
    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    
    trace.info()<<"read mesh..."<<std::endl;
    DGtal::Mesh<Z3i::RealPoint> mesh(true);
    MeshReader<Z3i::RealPoint>::importOFFFile(meshInput, mesh, false);
    trace.info()<<"read centerline..."<<std::endl;
    std::vector<Z3i::RealPoint> centerline = PointListReader<Z3i::RealPoint>::getPointsFromFile(centerlineInput);
    trace.info()<<"read cylindrical point..."<<std::endl;
    std::vector<Z3i::RealPoint> cylPoints = PointListReader<Z3i::RealPoint>::getPointsFromFile(cylPointsInput);
    /****************************/
    /*COMPUTE MAP DISCRETISATION*/
    /****************************/
    std::vector<double> radiusComp;
    for (const auto& point : cylPoints) {
        radiusComp.push_back(point[0]); // Accès à la première composante du vecteur cylindrique : (radius,angle,height)
    }
    UnrolledMap unrolled_map(radiusComp,MaxdF,minRV,intensityPerUnit, pad);
    trace.info()<<"compute discretisation..."<<std::endl;
    unrolled_map.computeDicretisation(cylPoints);
    unrolled_map.computeNormalizedImageMultiScale();
    unrolled_map.computeGRAYImage();
    Image2dGrayScale relief = unrolled_map.getReliefImageGrayScale();
    PGMWriter<Image2dGrayScale>::exportPGM("output.pgm", relief);
    
}