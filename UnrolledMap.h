//  UnrolledMap.h
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

#ifndef UNROLL_MAP
#define UNROLL_MAP

#include <utility>
#include <iostream>
#include <vector>

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/boards/Board2D.h"

using namespace DGtal;
  //Image of double to store normalized values
  typedef ImageContainerBySTLVector<Z2i::Domain, float> Image2dNormalized;
  //Image of Char to make a grayscale image
  typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2dGrayScale;
  //Image of Color to make a rgb image
  typedef ImageContainerBySTLVector<Z2i::Domain, Color > imageRGB;
  //Mask image binary value, use by multi resolution process
  typedef ImageContainerBySTLVector<Z2i::Domain, bool> Image2dBin;


/*use to creagte groundtruth*/
struct info_groundtruth{
	int height;     //height of the reliefmap
  int width;      //width of the reliefmap
  double upperBA; //upperbound of angle
  double shift;   //shift to aply on the angle (=2pi if complete trunk)
  double minH;    //minimal height
};

class UnrolledMap{
  public:
    /**
    Empty constructor
    **/
    UnrolledMap();
    /**
    Constructor.
   **/
    UnrolledMap(std::vector<double> DeltaDistance,int dF,int gs_ori,int gs_intensity,int pad)://std::vector<CylindricalPoint> CylindricalPoints,
      MaxdF(dF),
      reliefRepresentation(DeltaDistance),
      reliefImage(Z2i::Domain()),
      reliefImageGrayScale(Z2i::Domain()),
      reliefImageRGB(Z2i::Domain()),
      mask(Z2i::Domain()),
      maxIndTop(0),
      zero_level_intensity(gs_ori),
      intensity_per_cm(gs_intensity),
      minIndBot(0),
      pad(pad){};//CPoints(CylindricalPoints),

    /**
    Copy constructor
    **/
    UnrolledMap(const UnrolledMap &um):
      MaxdF(um.MaxdF),

      reliefRepresentation(um.reliefRepresentation),
      reliefImage(um.reliefImage),
      reliefImageGrayScale(um.reliefImageGrayScale),
      reliefImageRGB(um.reliefImageRGB),
      mask(um.mask),
      maxIndTop(um.maxIndTop),
      minIndBot(um.minIndBot),
      height(um.height),
      width(um.width){};//CPoints(um.CPoints),

    /**
    Compute discretisation
     **/
    void computeDicretisation(std::vector<Z3i::RealPoint> CPoints);
    /**
    compute normalized image of unrolled map with the black pixel complete by a multi scale analyse of unrolled map
     **/
    void computeNormalizedImageMultiScale();
    /**
    write rgb image from normalized image
    **/
    void computeRGBImage();
    /**
    compute rgb image from normalized image
    **/
    void computeGRAYImage();

    /**
    return the vector of ind at pos (i,j) in unrolled_surface with the decrease factor specify by df : 1/dF.
    **/
    std::vector<unsigned int > getIndPointsInLowerResolution(unsigned int i,unsigned int j,int dF);
    /**
    return the vector of ind at pos (i,j) in unrolled_surface. Need unrolled_sruface to be build
    **/
    std::vector<unsigned int > getPointsUnrolled_surface(unsigned int i,unsigned int j);

    /**
    * Crop bot and top of normalizd image
    **/
    void cropTopBotImage();
    /*********
    **GETTERS**
    **********/
    /**
    return rgb image
    **/
    Image2dGrayScale getReliefImageGrayScale();
    /**
    return rgb image
    **/
    imageRGB getReliefImageRGB();
    /**
    return discretisation vectors
    **/
    std::vector<std::vector<std::vector<unsigned int>>> getDiscretisation();
  protected:

    /**
     transforme reliefImage betwen [0,255] with min and max reliefREP
     **/
    Image2dGrayScale toGrayscaleImageMinMax();
    /**
     Normalize reliefImage betwen [0,255] with a fixed start and padding
     **/
    Image2dGrayScale toGrayscaleImageFixed(int intensityPerCm,double reliefValueforZero);

    
    /**
    return the max of relief representation. You can specify the resolution by dF : 1/dF.
    if df=1 return the mean at position (i,j) in unrolled surface.
    **/
    double maxReliefRepresentation(unsigned int i, unsigned int j,int dF);


    //unrolled surface representation : each cells contain some index points.
    std::vector<std::vector<std::vector<unsigned int>>> discretisationMap;
    std::vector<std::vector<std::vector<unsigned int>>> discretisationMap_full;
    //Represenation of the relief, radius of deltadiff.
    std::vector<double> reliefRepresentation;
    //Cylindricales Points
    //std::vector<CylindricalPoint> CPoints;
    //discretisation
    int height, width;
    //the maximum decrease factor for multi resolution research (2^n) with n = decreaseFactor
    int MaxdF;
    // angle discretisation to compute partial circumference
    int pad;
    double upperBA;
    double shift;
    double minH;
    //index of lines to be cropped
    unsigned int maxIndTop,minIndBot;
    //normalizedimage with DGtal
    Image2dNormalized reliefImage;
    //grayscale image with DGTAL
    Image2dGrayScale reliefImageGrayScale;
    //rgb image with DGTAL
    imageRGB reliefImageRGB;
    //relief value for the 0 level intensity in grayscale map
    int zero_level_intensity;
    //number of intensity value to represent 1 cm of relief
    int intensity_per_cm;
    //mask to counter cellule who don't need to be fill by multi resolution analysis
    Image2dBin mask;
};
#endif
