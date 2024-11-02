#include <iostream>
#include <set>
#include <map>
#include <utility>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "IOHelper.h"

using namespace DGtal;
//typedef DGtal::PointVector<3U, double, std::array<double, 6UL>> LongRealPoint;

void IOHelper::export2Text(const std::vector<DGtal::Z3i::RealPoint> &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        Z3i::RealPoint p = v[i];
        outStream << p[0] << " "<< p[1] << " "<< p[2]<<std::endl;
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<unsigned int> &idPcl, const std::vector<unsigned int> &idSector, const std::string &filename){
  std::ofstream outStream;
  outStream.open(filename.c_str(), std::ofstream::out);
  for(unsigned int i = 0; i < idPcl.size();i++){
    unsigned int idpoint=idPcl[i];
    unsigned int idsector=idSector[i];
    outStream << idpoint << " "<< idsector<<std::endl;
  }
  outStream.close();
}
void IOHelper::export2Text(const std::vector<DGtal::Z3i::RealPoint> &v,std::vector<Z3i::RealPoint>c, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(int i = 0; i < v.size();i++){

        Z3i::RealPoint p = v[i];
        Z3i::RealPoint col= c[i];
        outStream << p[0] << " "<< p[1] << " "<< p[2]<< " "
                  << col[0]<< " "<<col[1]<< " "<<col[2]<<std::endl;
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<std::vector< std::vector<unsigned int> > > &sectorDefectClusters, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    int defect_id=0;
    //for each sector
    for(int i = 0; i < sectorDefectClusters.size();i++){
      //for each defects
      for(int j = 0; j < sectorDefectClusters[i].size();j++){
        //out ID_sector, ID_ defect
        outStream<<i+1<< ","<<defect_id<< ",";
        //for each lobal ID
        for(int k = 0; k < sectorDefectClusters[i][j].size();k++){
          //all ID_points
          outStream<< sectorDefectClusters[i][j][k] << " ";
        }
        outStream<<std::endl;
        defect_id+=1;
      }
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<DGtal::Z3i::RealPoint> &pointCloud,
        const std::vector<unsigned int> &indices, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < indices.size();i++){
        Z3i::RealPoint p = pointCloud.at(indices.at(i));
        outStream << p[0] << " "<< p[1] << " "<< p[2]<<std::endl;
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<CylindricalPoint> &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        CylindricalPoint p = v[i];
        outStream << p.radius << " "<< p.angle << " "<< p.height <<std::endl;
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<Z3i::Point> &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        Z3i::Point p = v[i];
        outStream << p[0] << " "<< p[1] << " "<< p[2] <<std::endl;
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<int> &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        int p = v[i];
        outStream << p <<std::endl;
    }
    outStream.close();
}

void IOHelper::readCylindricalPointFromText(const std::string &filename,std::vector<CylindricalPoint> &v){
  std::ifstream infile;
  infile.open(filename.c_str(), std::ifstream::in);
  std::string delimiter=" ";
  std::string currentLine;
  CylindricalPoint mpCurrent;
  getline(infile, currentLine );
  while ( infile.good() ){
    std::vector<double> CylValue;
    std::stringstream linestream(currentLine);
    std::string value;

    while(getline(linestream,value,' ')){
      CylValue.push_back(std::stod(value));
    }
    if(CylValue.size()>0){
      mpCurrent.radius=CylValue.at(0);
      mpCurrent.angle=CylValue.at(1);
      mpCurrent.height=CylValue.at(2);
      v.push_back(mpCurrent);
    }

    getline(infile, currentLine);
  }

}
void IOHelper::export2Text(const std::vector<std::pair<double, double> > &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        std::pair<double, double> pa = v.at(i);
        outStream << pa.first << " "<< pa.second <<std::endl;
    }
    outStream.close();
}

/*
void IOHelper::export2Text(const std::vector<double> &xs, const std::vector<double> &ys, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    if( xs.size() == ys.size() ){
        for(unsigned int i = 0; i < xs.size();i++){
            outStream << xs.at(i) << " "<< ys.at(i) <<std::endl;
        }
    }
    outStream.close();

}
*/
void IOHelper::writeDiscretisationNormalesToFile(const std::vector<std::vector<std::vector<unsigned int>>> &discretisation,std::vector<Z3i::RealPoint> &normales,const std::string &fileName){
  trace.info()<<"Writting..."<< fileName<<std::endl;
  std::ofstream outStream;
  outStream.open(fileName.c_str(), std::ofstream::out);
  //first line is for dimension X/Y of discretisation
    outStream<<"#height and widht of discretisation"<<std::endl;
    outStream<<"#then each line for 1 pixel (-1 for empty celss)"<<std::endl;
    outStream<<discretisation[0].size()<<" "<<discretisation.size()<<std::endl;
  //then the discretisation
  int rows=discretisation[0].size();
  int cols=discretisation.size();

  for (int i=0; i < cols; ++i){
    for (int j=0; j < rows; ++j){

      std::vector<unsigned int> v =discretisation.at(i).at(j) ;
      if(!v.empty()){
        for(auto it = v.begin(); it != v.end(); ++it) {
          outStream << normales[*it][0]<<" "<<normales[*it][1]<<" "<<normales[*it][2]<<" ";
        }
        outStream<<std::endl;
      }else{
        outStream<<"-1 -1 -1 "<<std::endl;
      }
    }
  }
  outStream.close();
}
void IOHelper::readDistanceFromFile(const std::string &fileName, std::vector<double> &vectDistances){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
      if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
          vectDistances.push_back(std::stod(str));
      }
      getline(infile, str);
    }
  }


void IOHelper::export2OFF(const Mesh<Z3i::RealPoint> &mesh, std::string fileName){
    std::ofstream offMesh (fileName.c_str());
    DGtal::MeshWriter<Z3i::RealPoint>::export2OFF(offMesh, mesh);
    offMesh.close();
}


void IOHelper::export2OFF(const Mesh<Z3i::RealPoint> &mesh, const std::vector<unsigned int> &faceIds, const std::string &fileName){
    //slow
    std::set<unsigned int> vertexes;
    std::vector<unsigned int> vertexList;

    //map between old vertex id and new one
    std::map<unsigned int, unsigned int> vertexIdMap;
    unsigned int currentId = 0;
    for(unsigned int fId: faceIds){
        std::vector<long unsigned int>  aFace = mesh.getFace(fId);
        for(unsigned int indexVertex : aFace){
            std::pair<std::set<unsigned int>::iterator, bool> ret = vertexes.insert(indexVertex);
            //new item
            if( ret.second == true ){
                vertexIdMap[indexVertex] = currentId;
                vertexList.push_back(indexVertex);
                currentId++;
            }
        }
    }
    assert(vertexIdMap.size() == vertexList.size());
    std::ofstream out (fileName.c_str());
    out << "OFF"<< std::endl;
    out << vertexes.size()  << " " << faceIds.size() << " " << 0 << " " << std::endl;

    for(unsigned int vId: vertexList){
        out << mesh.getVertex(vId)[0] << " " << mesh.getVertex(vId)[1] << " "<< mesh.getVertex(vId)[2] << std::endl;
    }

    //print faces!!
    for(unsigned int fId: faceIds){
        std::vector<long unsigned int>  aFace = mesh.getFace(fId);
        out << aFace.size() << " " ;
        for(unsigned int indexVertex : aFace){
            out << vertexIdMap[indexVertex] << " " ;
        }
        DGtal::Color col = mesh.getFaceColor(fId);
        if( mesh.isStoringFaceColors() ){
            out << " ";
            out << ((double) col.red())/255.0 << " "
                << ((double) col.green())/255.0 << " "<< ((double) col.blue())/255.0
                << " " << ((double) col.alpha())/255.0 ;
        }
        out << std::endl;
    }
    out.close();

}

void IOHelper::readIntsFromFile(const std::string &fileName, std::vector<int> &rs){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
        if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
            rs.push_back(std::stoi(str));
        }
        getline(infile, str);
    }
}

void IOHelper::readIntsFromFile(const std::string &fileName, std::vector<int> &id, std::vector<int> &sector){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string line, val;
    std::getline (infile, line);
    while (std::getline (infile, line)) {                                        /* read each line */
      std::vector<std::string> v;                                       /* row vector v */
      std::stringstream s (line);                                       /* stringstream line */
      while (getline (s, val, ' '))                                     /* get each value (',' delimited) */
        v.push_back (std::string (val));                              /* add to row vector */
      id.push_back(std::stoi(v[0]));
      sector.push_back(std::stoi(v[1]));
   }
}

void IOHelper::readStringsFromFile(const std::string &fileName, std::vector<std::string> &rs){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
        if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
            rs.push_back(str);
        }
        getline(infile, str);
    }
}




void IOHelper::readFeatures(const std::string &fileName, std::vector<std::vector<float> > &rs){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
        if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
            std::vector<float> feature;
            std::istringstream s(str);
            float f;
            while (s >> f) {
                feature.push_back(f);
            }
            rs.push_back(feature);
        }
        getline(infile, str);
    }
}

void IOHelper::readXYZFile(const std::string &fileName, std::vector<DGtal::Z3i::RealPoint> &centerline){
    std::ifstream input{fileName};
    std::string line;
    std::getline(input, line);
    while(std::getline(input, line)){
      std::istringstream s(line);

      std::string field;
      int i =0;
      std::vector<float> tempVecPoints;
      while (getline(s, field,' ') && i < 3) {

        tempVecPoints.push_back(std::stof(field));
        i++;
      }
      centerline.push_back(Z3i::RealPoint (tempVecPoints[0],tempVecPoints[1],tempVecPoints[2]));
    }
}

void IOHelper::readDEFECTFile(const std::string &fileName, std::vector<std::vector<int> > &defs){
  std::ifstream input{fileName};
  std::string str;
  getline(input, str );

  while(input.good()){

    if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
      int i =0;
      std::istringstream s(str);
      std::string field;
      std::vector<int> defect_id;
      while (getline(s, field,',')) {
        if(i>1){
          std::istringstream s2(field);
          std::string field2;
          while (getline(s2, field2,' ')) {
            defect_id.push_back(std::stoi(field2));
          }
        }
        i+=1;
      }
      defs.push_back(defect_id);
    }
    getline(input, str);
  }
}
void IOHelper::writeDiscretisationToFile(const std::vector<std::vector<std::vector<unsigned int>>> &discretisation,const int &rowcroppedBot,const int &rowcroppedTop, const std::string &fileName){
  trace.info()<<"Writting discretisation ..."<<std::endl;
  std::ofstream outStream;
  outStream.open(fileName.c_str(), std::ofstream::out);
  //first line is for dimension X/Y of discretisation
  outStream<<discretisation.size()<<" "<<discretisation[0].size()<<std::endl;
  //Second line is for nb line cropped
  outStream<<rowcroppedBot<<" "<<rowcroppedTop<<std::endl;
  //then the discretisation
  int rows=discretisation.size();
  int cols=discretisation[0].size();

  for (int i=0; i < cols; ++i){
    for (int j=0; j < rows; ++j){
      std::vector<unsigned int> v =discretisation.at(j).at(i) ;
      if(!v.empty()){
        for(auto it = v.begin(); it != v.end(); ++it) {
          outStream << *it;
          outStream << " ";
        }
        outStream<<std::endl;
      }else{
        outStream<<"-1 "<<std::endl;
      }
    }
  }
  outStream.close();
}
void IOHelper::readDiscretisationFromFile(const std::string &fileName, std::vector<std::vector<std::vector<unsigned int>>> &discretisation, int &rowcroppedBot,int &rowcroppedTop){
  std::ifstream infile;
  infile.open(fileName.c_str(), std::ifstream::in);
  std::string currentLine;
  std::string delimiter=" ";
  std::string subLine;

  getline(infile, currentLine);
  //Here we want to read the dimension of the image. In discretisation.txt, dimension are delimited by a single space : dimY dimX. Dimension is located at the first line of the file
  int cols,rows;

  subLine=currentLine.substr(0, currentLine.find(delimiter));
  //trace.info()<<fileName<<std::endl;
  if(subLine.empty()){
    trace.info()<<"Problem in discretisation.txt. First line shoulde be : dimX dimY"<<std::endl;
  }else{
    rows=std::stoi(subLine);
  }
  subLine=currentLine.substr(currentLine.find(delimiter)+1,currentLine.size() );
  if(subLine.empty()){
    trace.info()<<"Problem in discretisation.txt. First line shoulde be : dimX dimY"<<std::endl;
  }else{
    cols=std::stoi(subLine);
  }
  //trace.info()<<"cols "<<cols<<std::endl;
  //trace.info()<<"rows "<<rows<<std::endl;
  //Here we want to read the number of rowcropped (second line of the file discretisation.txt)
  getline(infile, currentLine);
  if(currentLine.empty()){
    trace.info()<<"Problem in discretisation.txt. Second line shoulde be : numberRowCroppedFromBot numberRowCroppedFromTop"<<std::endl;
  }else{
    rowcroppedBot=std::stoi(currentLine);
  }
  subLine=currentLine.substr(currentLine.find(delimiter)+1,currentLine.size() );
  if(subLine.empty()){
    trace.info()<<"Problem in discretisation.txt. Second line shoulde be : numberRowCroppedFromBot numberRowCroppedFromTop"<<std::endl;
  }else{
    rowcroppedTop=std::stoi(subLine);
  }
  //Here we want to fill input discretisation with value in discretisation.txt starting at the trhird line
  //init dicretisation size with row and cols read before
  discretisation.resize(rows);
  for (int i = 0; i < rows; ++i){
      discretisation[i].resize(cols);
  }
  //loop on other line in the discretisation.txt
  size_t pos = 0;
  //actual position to add in discretisation
  int i=0,j=0;//i for cols and j for rows
  //there are exactly the same number of line in discretisation.txt than the number of cels un discretisation
  for (int i=0; i < cols; ++i){
    for (int j=0; j < rows; ++j){
      getline(infile, currentLine);
      std::stringstream linestream(currentLine);
      std::string value;
      while(getline(linestream,value,' ')){
        if(std::stod(value)!=-1){
          discretisation[j][i].push_back( std::stod(value));
        }
      }
    }
  }


}

void IOHelper::export2Text(info_groundtruth info, const std::string &filename){
  std::ofstream outStream;
  /*write reliefMap info*/
  outStream.open(filename, std::ofstream::out);
  outStream << "height" << " "<< "width" << " "<< "upperBA"<< " " << "shift" << "minH"<<std::endl;
  outStream << info.height << " "<< info.width << " "<< info.upperBA<< " " << info.shift << " " << info.minH<<std::endl;
  outStream.close();
}
 void IOHelper::readGroundTruthInfo( const std::string &filename, info_groundtruth &info){
    std::ifstream infile;
    infile.open(filename.c_str(), std::ifstream::in);
    std::string currentLine;
    std::string value;
    getline(infile, currentLine );//dont use line 0
    getline(infile, currentLine );
    std::stringstream linestream(currentLine);
    getline(linestream,value,' ');
    info.height=std::stoi(value);
    getline(linestream,value,' ');
    info.width=std::stoi(value);
    getline(linestream,value,' ');
    info.upperBA=std::stod(value);
    getline(linestream,value,' ');
    info.shift=std::stod(value);
    getline(linestream,value,' ');
    info.minH=std::stod(value);
    infile.close();
 }
