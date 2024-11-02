//  Main_Unroll.cpp
//
//  Created by Florian Delconte on 09/07/2021.
//  Copyright © 2021 xxx. All rights reserved.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#include <ctime>
#include <iostream>
#include <fstream>
#include <utility>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/writers/STBWriter.h"
#include "DGtal/io/readers/STBReader.h"
#include "DGtal/io/writers/PGMWriter.h"
#include "DGtal/io/readers/PGMReader.h"

#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/shapes/Mesh.h"

#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"


#include "segmentation/Patch.h"
#include "segmentation/DefectSegmentation.h"
#include "segmentation/SegmentationHelper.h"


#include "Centerline/Centerline.h"
#include "Centerline/CenterlineHelper.h"

#include "Common/CylindricalCoordinateSystem.h"
#include "Common/IOHelper.h"



#include "CLI11.hpp"


using namespace DGtal;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;


std::string REPO="/volWork/these/source/examples/TH/FINAL/";
//std::string XYZ_BRUT=REPO+"xyz_brut/";
std::string TRUNK_MESH_CLEAN="/volWork/these/source/examples/TH/ALL/clean1_Log/log/trunk/trunk_s15/mesh/";
//std::string BRANCH=REPO+"branch/";
std::string B_T_C_cloud=REPO+"b_t_c_cloud/";
std::string CENTERLINE=REPO+"centerline/";
//std::string COLOR=REPO+"color/";
//std::string LUMINANCE=REPO+"luminance/";
//std::string NORMALS=REPO+"normals/";
std::string OUTPUT_MAP_RELIEF=REPO+"XYZ_MAP/rm/";
std::string OUTPUT_MAP_HSV=REPO+"mesh_map/hsv/";
std::string OUTPUT_MAP_VG=REPO+"mesh_map/vg/";
std::string OUTPUT_DISCRETISATION=REPO+"discretisation/";

struct IdColor{
  Color operator()( const unsigned int & aValue ) const{
    return DGtal::Color(aValue);
  }
};
//TODO:fuse Face of DGtal;
struct facept{
	Z3i::RealPoint p0;
  Z3i::RealPoint p1;
  Z3i::RealPoint p2;
};

DGtal::Color
colorFromHSB(double h, double saturation, double value){
  double r, g, b;
  DGtal::Color::HSVtoRGB(r, g, b, h,saturation, value);
  return DGtal::Color(r*255.0,g*255.0,b*255.0);

}

void
computeNormal(std::vector<Z3i::RealPoint> &vectorOfNormals,DGtal::Mesh<Z3i::RealPoint> mesh, std::vector<CylindricalPoint> cylPts, std::vector<Z3i::RealPoint> cl, std::vector<Z3i::RealPoint>pcl) {
  trace.info()<<"Start compute normal ..."<<std::endl;
  for (int i = 0; i < mesh.nbFaces(); i++){
    Face currentFace = mesh.getFace(i);
    
    Z3i::RealPoint currentFaceP0 = mesh.getVertex(currentFace.at(0));
    Z3i::RealPoint currentFaceP1 = mesh.getVertex(currentFace.at(1));
    Z3i::RealPoint currentFaceP2 = mesh.getVertex(currentFace.at(2));
    Z3i::RealPoint center = (currentFaceP0+currentFaceP1+currentFaceP2)/3.0;
    Z3i::RealPoint v0 = currentFaceP0 - currentFaceP2;
    Z3i::RealPoint v1 = currentFaceP1 - currentFaceP2;
    Z3i::RealPoint vnormal=v0.crossProduct(v1);
      Z3i::RealPoint vnormalised(vnormal/vnormal.norm());
    //truc
    int segmentID=cylPts.at(i).segmentId;
    //trace.info()<<i<<" "<<segmentID<<std::endl;
    Z3i::RealPoint directionFiber = cl[segmentID + 1] - cl[segmentID];
      directionFiber = directionFiber / directionFiber.norm();
    double dist = directionFiber.dot(pcl[i] - cl[segmentID]);
    Z3i::RealPoint proj = cl[segmentID] + dist*directionFiber;
    Z3i::RealPoint vectorRadial = pcl[i] - proj;
      Z3i::RealPoint vectorRadial_normalized = (vectorRadial/vectorRadial.norm());
    double angle_radial_normal=acos((vnormalised.dot(vectorRadial_normalized))/(vnormalised.norm()*vectorRadial_normalized.norm()));
    double chroma=sin(angle_radial_normal)*vnormal.norm();

    if (vnormalised.dot(vectorRadial_normalized) < 0){
      vnormalised *= -1.0;
    }
    Z3i::RealPoint faceSvect = vectorRadial_normalized.crossProduct(Z3i::RealPoint(0.0,0.0,1.0));
      faceSvect = faceSvect/faceSvect.norm();

    vectorOfNormals[i] = Z3i::RealPoint(Z3i::RealPoint(0.0,0.0,1.0).dot(vnormalised),
                                          faceSvect.dot(vnormalised),
                                          vectorRadial_normalized.dot(vnormalised));
    vectorOfNormals[i] = vectorOfNormals[i]/vectorOfNormals[i].norm();
  }

}

void
computeNormal(std::vector<Z3i::RealPoint> &vectorOfNormals,std::vector<facept> fom, std::vector<CylindricalPoint> cylPts, std::vector<Z3i::RealPoint> cl, std::vector<Z3i::RealPoint>pcl) {
  trace.info()<<"Start compute normal ..."<<std::endl;
  for (int i = 0; i < fom.size(); i++){
    facept currentFace = fom.at(i);
    Z3i::RealPoint currentFaceP0 = currentFace.p0;
    Z3i::RealPoint currentFaceP1 = currentFace.p1;
    Z3i::RealPoint currentFaceP2 = currentFace.p2;
    Z3i::RealPoint center = (currentFaceP0+currentFaceP1+currentFaceP2)/3.0;
    Z3i::RealPoint v0 = currentFaceP0 - currentFaceP2;
    Z3i::RealPoint v1 = currentFaceP1 - currentFaceP2;
    Z3i::RealPoint vnormal=v0.crossProduct(v1);
    Z3i::RealPoint vnormalised(vnormal/vnormal.norm());
    //truc
    int segmentID=cylPts.at(i).segmentId;
    //trace.info()<<i<<" "<<segmentID<<std::endl;
    Z3i::RealPoint directionFiber = cl[segmentID + 1] - cl[segmentID];
      directionFiber = directionFiber / directionFiber.norm();
    double dist = directionFiber.dot(pcl[i] - cl[segmentID]);
    Z3i::RealPoint proj = cl[segmentID] + dist*directionFiber;
    Z3i::RealPoint vectorRadial = pcl[i] - proj;
      Z3i::RealPoint vectorRadial_normalized = (vectorRadial/vectorRadial.norm());
    double angle_radial_normal=acos((vnormalised.dot(vectorRadial_normalized))/(vnormalised.norm()*vectorRadial_normalized.norm()));
    double chroma=sin(angle_radial_normal)*vnormal.norm();

    if (vnormalised.dot(vectorRadial_normalized) < 0){
      vnormalised *= -1.0;
    }
    Z3i::RealPoint faceSvect = vectorRadial_normalized.crossProduct(Z3i::RealPoint(0.0,0.0,1.0));
      faceSvect = faceSvect/faceSvect.norm();

    vectorOfNormals[i] = Z3i::RealPoint(Z3i::RealPoint(0.0,0.0,1.0).dot(vnormalised),
                                          faceSvect.dot(vnormalised),
                                          vectorRadial_normalized.dot(vnormalised));
    vectorOfNormals[i] = vectorOfNormals[i]/vectorOfNormals[i].norm();
  }

}
void
computeNormal(std::vector<Z3i::RealPoint> &outputNormals,std::vector<Z3i::RealPoint> basicNormals, std::vector<CylindricalPoint> cylPts, std::vector<Z3i::RealPoint> cl, std::vector<Z3i::RealPoint>pcl) {
  trace.info()<<"Start compute normal ..."<<std::endl;
  for (int i = 0; i < basicNormals.size(); i++){
    Z3i::RealPoint vnormal = basicNormals.at(i);
    Z3i::RealPoint vnormalised(vnormal/vnormal.norm());
    //truc
    int segmentID=cylPts.at(i).segmentId;
    //trace.info()<<i<<" "<<segmentID<<std::endl;
    Z3i::RealPoint directionFiber = cl[segmentID + 1] - cl[segmentID];
    directionFiber = directionFiber / directionFiber.norm();
    double dist = directionFiber.dot(pcl[i] - cl[segmentID]);
    Z3i::RealPoint proj = cl[segmentID] + dist*directionFiber;
    Z3i::RealPoint vectorRadial = pcl[i] - proj;
    Z3i::RealPoint vectorRadial_normalized = (vectorRadial/vectorRadial.norm());
    double angle_radial_normal=acos((vnormalised.dot(vectorRadial_normalized))/(vnormalised.norm()*vectorRadial_normalized.norm()));
    double chroma=sin(angle_radial_normal)*vnormal.norm();

    if (vnormalised.dot(vectorRadial_normalized) < 0){
      vnormalised *= -1.0;
    }
    Z3i::RealPoint faceSvect = vectorRadial_normalized.crossProduct(Z3i::RealPoint(0.0,0.0,1.0));
      faceSvect = faceSvect/faceSvect.norm();

    outputNormals[i] = Z3i::RealPoint(Z3i::RealPoint(0.0,0.0,1.0).dot(vnormalised),
                                          faceSvect.dot(vnormalised),
                                          vectorRadial_normalized.dot(vnormalised));
    outputNormals[i] = outputNormals[i]/outputNormals[i].norm();
  }

}
/*Not used : dgtal image can be defined with negatif domain*/
void
csvTo2dVector(std::string csvFile,std::vector<std::vector<std::string>> &array){
  std::ifstream f;
  f.open(csvFile);
  if (! f.is_open()) {
      std::cerr << "error: file open failed '" << csvFile << "'.\n";
      exit(1);
  }
  std::string line, val;
  while (std::getline (f, line)) {                                        /* read each line */
        std::vector<std::string> v;                                       /* row vector v */
        std::stringstream s (line);                                       /* stringstream line */
        while (getline (s, val, ';'))                                     /* get each value (',' delimited) */
            v.push_back (std::string (val));                              /* add to row vector */
        array.push_back (v);                                              /* add row vector to array */
  }
  f.close();
}

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
  std::string logPrefixe="DO5";
  std::string inputCenterline;
  int centerline_voxel_size{30};
  double accRadius {350.0};
  double searchRadius {30.0};
  int nbControlPoint {2};
  double binWidth {0.001};
  double confThreshold {0.0};
  double patchWidth {25.0};
  int patchHeight {100};
  int voxelSize {3};
  int MaxdF {4};
  int minRV {-5};
  int intensityPerCm {5};
  int nbSegment {10};
  double sectorLength {50};
  int pad {200};
  std::string outputPrefix;
  bool hsvRending=false;
  bool vgRending=false;

  DGtal::Mesh<Z3i::RealPoint> scaledMesh(true);

  bool invertNormal {false};

  DGtal::Mesh<Z3i::RealPoint> oriMesh(true);
  DGtal::Mesh<Z3i::RealPoint> testMesh(true);
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Allowed options are: ");
  app.add_option("-p,--input,1", logPrefixe , "log prefixe");
  app.add_option("--accRadius,-r", accRadius, "accumulation radius.", true);
  app.add_option("--searchRadius,-R", searchRadius, "search for neighbor in radius.", true);
  app.add_option("--nbControl", nbControlPoint, "Number of control points for bsplines", true);
  app.add_flag("--invertNormal,-n", invertNormal, "invert normal to apply accumulation.");
  app.add_option("--binWidth,-b", binWidth, "bin width used to compute threshold", true);
  app.add_option("--confThreshold,-t", confThreshold, "threshold in the confidence estimation.", true );
  app.add_option("--patchWidth,-a",patchWidth,"Arc length/ width of patch", true);
  app.add_option("--patchHeight,-e",patchHeight,"Height of patch", true);
  app.add_option("--voxelSize,-v", voxelSize, "Voxel size", true);
  app.add_option("--centerline,-c", inputCenterline, "Centerline file of log or trunk", true);
  app.add_option("--centerlineSubSample,-s", inputCenterline, "Centerline file of log or trunk", true);
  app.add_option("--MaxdF,-d", MaxdF, "Max decrease factor for multi resolution search", true );
  app.add_option("--minRV", minRV, "relief value for 0 level in grayscale intensity", true );
  app.add_option("--intensityPerCm", intensityPerCm, "number of grayscale intensity to represente 1cm of relief", true);
  app.add_option("--nbSegment", nbSegment, "Number of segments to compute the centerline", true);
  app.add_option("--sectorLength", sectorLength, "used to segment branch", true);
  app.add_option("--pad", pad, "used to discretized angle", true);
  app.add_flag("--hsvRending", hsvRending, "use  HSV rending (given from the normal vector)");
  app.add_flag("--vgRending", vgRending, "use  Video Game rending (given from the normal vector)");
  app.add_option("--output,-o,2", outputPrefix, "output prefix: output-defect.off, output-def-faces-ids, ...");
  //->required();

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  
  
  assert(voxelSize > 0);
  double scale=400;
  //trace.info()<<CENTERLINE+logPrefixe+"_centerline.xyz"<<std::endl;
  //trace.info()<<TRUNK_MESH_CLEAN+logPrefixe+".off"<<std::endl;
  
  //std::vector<Z3i::RealPoint> centerline = PointListReader<Z3i::RealPoint>::getPointsFromFile(CENTERLINE+logPrefixe+"_centerline.xyz");
  std::vector<Z3i::RealPoint> centerline = PointListReader<Z3i::RealPoint>::getPointsFromFile("/volWork/these/source/BERTRAND/SampleSujetM2/SampleSujetM2/chene1_centerline.xyz");
  DGtal::Mesh<Z3i::RealPoint> mesh(true);
  //MeshReader<Z3i::RealPoint>::importOFFFile(TRUNK_MESH_CLEAN+logPrefixe+".off", mesh, false);
  MeshReader<Z3i::RealPoint>::importOFFFile("/volWork/these/source/BERTRAND/SampleSujetM2/SampleSujetM2/chene1.off", mesh, false);
  std::vector<Z3i::RealPoint> pointCloud_trunk=std::vector<Z3i::RealPoint> ();
  for (int i = 0; i < mesh.nbFaces(); i++){
    
  }
  for (int i = 0; i < mesh.nbVertex(); i++){
      pointCloud_trunk.push_back(mesh.getVertex(i));
  }
  /******************/
  /*SCALE CENTERLINE*/
  /******************/
  //scaleCloud(centerline,scale);
  //scaleCloud(pointCloud_trunk,scale);
  /***********************/
  /*CENTERLINE SUB SAMPLE*/
  /***********************/
  std::vector<Z3i::RealPoint> centerline_subSampled;
  SegmentationHelper::simpleSubSample(centerline, centerline_voxel_size, centerline_subSampled);
  /********************/
  /*COMPUTE CYL POINTS*/
  /********************/
  std::vector<CylindricalPoint> cylPoints= std::vector<CylindricalPoint>();
  //compute
  /*CylindricalCoordinateSystem ccs(centerline_subSampled, Z3i::RealPoint(0.0,0.0,0.0));
  trace.info()<<"centerline_subSampled : "<<  centerline_subSampled.size()<<std::endl;
  for(unsigned int i = 0; i < pointCloud_trunk.size(); i++){
    trace.progressBar(i, pointCloud_trunk.size());
    Z3i::RealPoint xyzP=pointCloud_trunk.at(i);
    CylindricalPoint cylP = ccs.xyz2Cylindrical(xyzP);
    cylPoints[i]=cylP;
  }*/
  //read
  IOHelper::readCylindricalPointFromText("/volWork/these/source/BERTRAND/SampleSujetM2/SampleSujetM2/chene1-cyl",cylPoints);
  /************************/
  /*COMPUTE DELTA DISTANCE*/
  /************************/
  std::vector<double> distances (pointCloud_trunk.size());
  //trace.info()<<distances.size()<<std::endl;
  /*
  Patch pa(pointCloud_trunk,cylPoints,centerline_subSampled,patchWidth,patchHeight,binWidth,sectorLength,25);
  pa.init2();
  distances = pa.getDeltaDist();
  */
  
  /****************************/
  /*COMPUTE MAP DISCRETISATION*/
  /****************************/
  trace.info()<<pointCloud_trunk.size()<<std::endl;
  trace.info()<<distances.size()<<std::endl;
  trace.info()<<cylPoints.size()<<std::endl;
  trace.info() << pointCloud_trunk[50]<<std::endl;
  trace.info() << cylPoints[50].height<<" "<<cylPoints[50].angle<< " " <<cylPoints[50].radius<<std::endl;
  
  
  trace.info() << centerline_subSampled[4]<<std::endl;
  UnrolledMap unrolled_map(distances,MaxdF,minRV,intensityPerCm, pad);
  unrolled_map.computeDicretisation(cylPoints);
  //unrolled_map.majBinaryMask(6);
  unrolled_map.computeNormalizedImageMultiScale();
  unrolled_map.computeGRAYImage();
  Image2dGrayScale relief = unrolled_map.getReliefImageGrayScale();
  //PGMWriter<Image2dGrayScale>::exportPGM(OUTPUT_MAP_RELIEF+logPrefixe+".pgm", relief);
  PGMWriter<Image2dGrayScale>::exportPGM("output.pgm", relief);
  //relief>>OUTPUT_MAP_RELIEF+logPrefixe+".pgm";
  /**********************************/
  /*COMPUTE NORMALS MAP (HSV and VG)*/
  /**********************************/ 
  std::vector<Z3i::RealPoint> localNormal(pointCloud_trunk.size());
  computeNormal(localNormal,mesh, cylPoints, centerline_subSampled,pointCloud_trunk);
  std::vector<std::vector<std::vector<unsigned int>>> discretisationMap=unrolled_map.getDiscretisation();
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned int> Image2DC;
  typedef ImageContainerBySTLVector < Z2i::Domain, Color> Image2DColor;

  Image2DColor map_normals_hsv (relief.domain());
  Image2DColor map_normals_vg (relief.domain());
  

  for (int i=0; i < discretisationMap.size(); ++i){
    for (int j=0; j < discretisationMap[0].size(); ++j){
      std::vector<unsigned int> v =discretisationMap.at(i).at(j) ;
      if(!v.empty()){
        int id=discretisationMap.at(i).at(j).at(0);
        Z3i::RealPoint n=localNormal.at(id);
       
        for(int h=1;h<discretisationMap.at(i).at(j).size();h++){
          id=discretisationMap.at(i).at(j).at(0);
          n+=localNormal.at(id);
        }
        n=n/discretisationMap.at(i).at(j).size();
         
        //HSV
        double sat = 1.0*( sin(acos(Z3i::RealPoint(0.0,0.0,1.0).dot(n))));
        double value = 1.0;
        double hue = ((int)(((2.0*M_PI+atan2(Z3i::RealPoint(0.0,1.0,0.0).dot(n),Z3i::RealPoint(1.0,0.0,0.0).dot(n)))/(2.0*M_PI))*360.0+100))%360;
        DGtal::Color colCode = colorFromHSB(hue, sat, value);
        map_normals_hsv.setValue(Z2i::Point(j,i), colCode);
        //VG
        DGtal::functors::Rescaling<double, unsigned int> rgRescale (-1.0, 1.0, 0, 255);
        DGtal::functors::Rescaling<double, unsigned int> bRescale (0.0, -1.0, 128, 255);
        DGtal::Color c (rgRescale(n[0]), rgRescale(n[1]), bRescale(n[2]) );
        map_normals_vg.setValue(Z2i::Point(j,i), c);
      }
    }
  }
  //IdColor id;
  //relief>>OUTPUT_MAP_RELIEF+logPrefixe+".pgm";
  //PPMWriter<Image2DC, IdColor  >::exportPPM(OUTPUT_MAP_HSV+logPrefixe+".ppm", map_normals_hsv, id);
  //PPMWriter<Image2DC, IdColor  >::exportPPM(OUTPUT_MAP_VG+logPrefixe+".ppm", map_normals_vg, id);
  
  //TODO : regarder pourquoi sa inverse le sens.
  //STBWriter<Image2DColor>::exportPNG(OUTPUT_MAP_HSV+logPrefixe+".png", map_normals_hsv);
  STBWriter<Image2DColor>::exportPNG("output_hsv.png", map_normals_hsv);
  //STBWriter<Image2DColor>::exportPNG(OUTPUT_MAP_VG+logPrefixe+".png", map_normals_vg);
  STBWriter<Image2DColor>::exportPNG("output_vg.png", map_normals_vg);
  
  /*******************************/
  /*WRITE MAP INFO DISCRETISATION*/
  /*******************************/
  //IOHelper::writeDiscretisationToFile(unrolled_map.getDiscretisation(),0,0,OUTPUT_DISCRETISATION+logPrefixe+"_discretisation");
  //IOHelper::export2Text(unrolled_map.getinfo(),OUTPUT_DISCRETISATION+logPrefixe+"_mapInfo");
  return 0;


}