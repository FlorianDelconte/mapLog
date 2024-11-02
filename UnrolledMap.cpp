//  UnrolledMap.cpp
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

#include "UnrolledMap.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ArrayImageAdapter.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"



using namespace DGtal;

/*
See Alg 1.
*/
void
UnrolledMap::computeDicretisation(std::vector<Z3i::RealPoint> CPoints){
    trace.info()<<"begin discretisation"<<std::endl;


    /*Declaration*/
    auto minMaxHeight = std::minmax_element(CPoints.begin(), CPoints.end(),
                                       [](const Z3i::RealPoint& a, const Z3i::RealPoint& b) {
                                           return a[2] < b[2]; 
                                       });
   
    minH = (*minMaxHeight.first)[2];
    double maxH = (*minMaxHeight.second)[2];

    Z3i::RealPoint mpCurrent;
    double meanR=0.;
    for(unsigned int i = 0; i < CPoints.size(); i++){
        mpCurrent=CPoints.at(i);
        meanR+=mpCurrent[0];
    }
    meanR/=CPoints.size();
    std::vector<unsigned int>T(pad);
    /***************/
    height=round(maxH-minH+1);                                                            //Lenght of relief map, bound to nearest

    double IIPi = 2*M_PI;
    int i;
    for (Z3i::RealPoint p : CPoints){
        i = round(((p[1] * (pad-1)) / IIPi)) ;
        T.at(i)=1;
    }
   
    bool find;                                       //Boolean value to search '01' pattern in T
    int firstIndexNoNull=0;                         //Store the 1 index in '01' pattern found
    unsigned int id=0;                            //Current index to loop over T
    int nbZero=0;                               //Count zero values
    while ( id< pad && !find){
        int currentValue=T.at(id);
        int nextValue=T.at((id+1)%(pad-1));
        if(currentValue==0 && nextValue==1){
            firstIndexNoNull=id+1;
            find=true;
        }
        //count nb zero
        if(currentValue==0){
            nbZero+=1;
        }
        id+=1;
    }

    width=round(2*M_PI*(double(pad-nbZero)/pad)*meanR);                     //Width of the relief map
    
    upperBA=double(pad-nbZero)/pad*(2*M_PI);                                     //Upper bound of angle in CPoints
    shift=(2*M_PI)-((double(firstIndexNoNull)/pad)*(2*M_PI));                    //Shift to apply on CPoints

    //preallocate size for unrolled map
    discretisationMap.resize(height);
    discretisationMap_full.resize(height);
    for (int i = 0; i < height; ++i){
        discretisationMap[i].resize(width);
        discretisationMap_full[i].resize(width);
    }
    int posX, posY;
    double shiftedCurrentAngle;
    trace.info()<<"size of discretisation : "<<"["<<height-1<<" ; "<<width-1<<" ]"<<std::endl;
    for(unsigned int i = 0; i < CPoints.size(); i++){
        mpCurrent=CPoints.at(i);
        shiftedCurrentAngle= std::fmod(mpCurrent[1]+shift,2*M_PI);                      //Apply shift and use modulus
        if(shiftedCurrentAngle>upperBA ||shiftedCurrentAngle<0){
            //std::cout<<"currentAngle :"<<shiftedCurrentAngle<<std::endl;
            //TODO : check why sometime shift is'nt good (no need to make a IF ELSE NORMALY)
        }else{
          posX=round(((width-1)/(upperBA))*(shiftedCurrentAngle-upperBA)+(width-1));          //bound to nearest
          posY=round(mpCurrent[2]-minH);                                                        //bound to nearest
          discretisationMap[posY][posX].push_back(i);
          discretisationMap_full[posY][posX].push_back(i);                                          //push current index in the cell
        }
    }
    trace.info()<<"end discretisation"<<std::endl;
}
/*
See Alg 2.
*/
void
UnrolledMap::computeNormalizedImageMultiScale(){
    trace.info()<<"Compute normalized image in multi scale ..."<<std::endl;
    time_t begin,end; 
    time (&begin);
    
    Z2i::Domain domain(Z2i::Point(0,0), Z2i::Point(width-1,height-1));
    reliefImage = Image2dNormalized(domain);
    int dF=1;
    bool needDisplay=false;
    double relief;
    trace.info()<<"MaxdF : "<<MaxdF<<std::endl;
    trace.info()<<"domain : "<<domain<<std::endl;
    relief=-1;
    for(unsigned int i = 0; i < height; i++){
        trace.progressBar(i, height);
        for(unsigned int j = 0; j < width; j++){
                relief=maxReliefRepresentation(i,j,0);         //relief of the current cell at 1/2^0
                dF=1;
                while(relief==-1 && dF<MaxdF){
                    relief=maxReliefRepresentation(i,j,dF);     //relief of the current cell at 1/2^dF
                    dF+=1;
                }

                reliefImage.setValue(Z2i::Point(j,i),relief);
        }
    }
    trace.info()<<" "<<std::endl;
    time (&end);
    double difference = difftime (end,begin);
    printf ("--------time : MULTI RESOLUTION PROCESS %.2lf seconds.\n", difference );
}

/*
See Alg 4.
*/
double
UnrolledMap::maxReliefRepresentation(unsigned int i, unsigned int j,int dF){
    //std::cout<<"------------maxReliefRepresentation()-----------------"<<std::endl;
    double maxRelief=INT_MIN;
    double currentRelief;
    std::vector<unsigned int> IndexesInCell = getIndPointsInLowerResolution(i,j,dF);
    //std::cout<<"IndexesInCell : "<<"OK"<<std::endl;
    if(!IndexesInCell.empty()){
        unsigned int IndP;
        //std::cout<<"-"<<std::endl;
        //std::cout<<IndexesInCell.size()<<std::endl;
        //std::cout<<discretisationMap.at(i).at(j).size()<<" -> ";
        for(std::vector<unsigned int>::iterator it = std::begin(IndexesInCell); it != std::end(IndexesInCell); ++it) {
            IndP=*it;
            currentRelief=reliefRepresentation.at(IndP);

            if(currentRelief>maxRelief){
                maxRelief=currentRelief;
            }
            discretisationMap_full[i][j].push_back(IndP);
        }
        //std::cout<<discretisationMap.at(i).at(j).size()<<std::endl;
        /*std::cout<<discretisationMap.at(i).at(j).size()<<std::endl;
        std::cout<<IndexesInCell.size()<<std::endl;
        for(int z =0;z<IndexesInCell.size();z++){
            discretisationMap.at(i).at(j).push_back(IndexesInCell.at(z));
        }
        std::cout<<discretisationMap.at(i).at(j).size()<<std::endl;exit(1);*/
        //discretisationMap.at(i).at(j)=IndexesInCell;
    }else{
        maxRelief=-1;                                                                               //conditional value
    }
    return maxRelief;
}


/*
See Alg 5.
*/
std::vector<unsigned int >
UnrolledMap::getIndPointsInLowerResolution(unsigned int i,unsigned int j,int dF){
    int down_pad=pow(2,dF);
    std::vector<unsigned int > outPutInd;
    int i_topLeftCorner,j_topLeftCorner;
    i_topLeftCorner=(i/down_pad)*down_pad;
    j_topLeftCorner=(j/down_pad)*down_pad;                                                                                            //round to lowest integer
    //loop on cells (of unrolledSurace) in region containing (i,j)
    //std::cout<<"down_pad : "<<down_pad<<std::endl;
    //std::cout<<"i_topLeftCorner : "<<i_topLeftCorner<<std::endl;
    //std::cout<<"j_topLeftCorner : "<<j_topLeftCorner<<std::endl;
    for( int k = i_topLeftCorner; k < (i_topLeftCorner+down_pad); k++){//+down_pad
        for( int l = j_topLeftCorner; l < (j_topLeftCorner+down_pad); l++){//+down_pad
            //check if cells of region is in unrolled surface -> no segmentation fault
            if((k < height)&&(l < width)){
                std::vector<unsigned int > temp =discretisationMap[k][l];
                for(int z =0;z<temp.size();z++){
                    outPutInd.push_back(temp[z]);
                }
                //unlarge  vector size
                //outPutInd.reserve(outPutInd.size() + discretisationMap[k][l].size());
                //concat vector
                //outPutInd.insert(outPutInd.end(), discretisationMap[k][l].begin(),  discretisationMap[k][l].end());         //push indexes at pos (k,l)
            }
        }
    }
    return outPutInd;
}

/*
See Alg 6.
*/
Image2dGrayScale
UnrolledMap::toGrayscaleImageFixed(int intensityPerCm, double minRV){

    double maxRV=(255*(1./intensityPerCm))+minRV;
    int grayscaleValue;
    double reliefValue;
    //CREATE GRAYSCALE IMAGE
    Image2dGrayScale grayScaleReliefImage=Image2dGrayScale(reliefImage.domain());
    std::cout<<"GRAY : "<<reliefImage.domain()<<std::endl;
    //FILL THE GRAYSCALLE IMAGE
    for ( auto pixel : reliefImage.domain() ){
      //if(mask(pixel)){
        reliefValue=reliefImage(pixel);
        if(reliefValue<minRV){
            grayscaleValue=0;
        }else if(reliefValue>maxRV){
            grayscaleValue=255;
        }else{
            grayscaleValue=((reliefValue-minRV)/(maxRV-minRV))*255;
        }
        grayScaleReliefImage.setValue(pixel,grayscaleValue);
      //}
    }
    return grayScaleReliefImage;
}

Image2dGrayScale
UnrolledMap::toGrayscaleImageMinMax(){
    int grayscaleValue;
    double reliefValue;
    /*DGtal*/
    Image2dGrayScale grayScaleReliefImage=Image2dGrayScale(reliefImage.domain());
    //SEARCH MIN MAX IN RELIEFIMAGE
    double minDG, maxDG;
    auto minMaxHeight = std::minmax_element(reliefRepresentation.begin(), reliefRepresentation.end());
    minDG = (*minMaxHeight.first);
    maxDG = (*minMaxHeight.second);
    trace.info()<<minDG<<std::endl;
    trace.info()<<maxDG<<std::endl;
    for ( auto point : reliefImage.domain() ){
        
        reliefValue=reliefImage(point);
       
        grayscaleValue=((reliefValue-minDG)/(maxDG-minDG))*255;
        
        grayScaleReliefImage.setValue(point,grayscaleValue);
        
    }

    return grayScaleReliefImage;
}

void
UnrolledMap::computeRGBImage(){
    double reliefValue;
    Color reliefColor;
    //first : compute in grayScale
    Image2dGrayScale reliefGray = toGrayscaleImageFixed(intensity_per_cm,zero_level_intensity);
    //second : convert grayscale to rgb
    imageRGB imagergb=imageRGB(reliefGray.domain());

    float min=*min_element(reliefGray.range().begin(), reliefGray.range().end());
    float max=*max_element(reliefGray.range().begin(), reliefGray.range().end());
    GradientColorMap<float,CMAP_COPPER> gradient( min, max);
    for ( auto point : reliefGray.domain() ){
        reliefValue=reliefGray(point);
        
        reliefColor=gradient(reliefValue);

        imagergb.setValue(point,reliefColor);
    }

    reliefImageRGB=imagergb;
}
void
UnrolledMap::computeGRAYImage(){
    if(intensity_per_cm==-1){
      trace.info()<<"make grayscale image with [min ;max]"<< std::endl;
      reliefImageGrayScale=toGrayscaleImageMinMax();
    }else{
      trace.info()<<"make grayscale with zero = "<<zero_level_intensity<< " and pad = "<<intensity_per_cm<< std::endl;
      reliefImageGrayScale=toGrayscaleImageFixed(intensity_per_cm,zero_level_intensity);
    }

}

Image2dGrayScale
UnrolledMap::getReliefImageGrayScale(){
    return reliefImageGrayScale;
}


imageRGB
UnrolledMap::getReliefImageRGB(){
    return reliefImageRGB;
}

std::vector<std::vector<std::vector<unsigned int>>>
UnrolledMap::getDiscretisation(){
  return discretisationMap_full;
}
