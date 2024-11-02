# TLDDC_IPOL
Repository of work submitted to IPOL 2021:
**CNN-based Method for Segmenting Tree SurfaceSingularites**
## Dependencies
There is an installation procedure (tested on Ubuntu 20.04) below. July 2, 2021.
First, you have to install the dependencies linked to the c++ code:

###### DGtal version 1.1.0 or later, see [DGtal installation] (https://github.com/DGtal-team/DGtal/blob/master/README.md) an source source https://github.com/DGtal-team/DGtal
###### Update system
```
sudo apt update
sudo apt upgrade
```
###### Install DGtal dependencies (for example) :

```
sudo apt install build-essential
sudo apt install zlib1g
sudo apt-get install libboost-all-dev
sudo apt install cmake
sudo apt install git
```

###### Build DGtal :
```
git clone https://github.com/DGtal-team/DGtal.git
cd DGTALSOURCES
mkdir build
cd build
cmake ..
make
```
###### GNU GSL
```
sudo apt install libgsl-dev
```
###### PCL
```
sudo apt install libpcl-dev
```
###### Build our c++ code :
```
mkdir build
cd build
cmake ..  -DDGtal_DIR=/path/to/DGtal
make
```
Then, you have to install the dependencies linked to the python code (we recommend the reviewer to use un virtual environnement like described in the tensorflow guide):
###### python3
```
sudo apt install python3-dev python3-pip python3-venv
```
###### create a virtual environnement
```
cd ~
python3 -m venv --system-site-packages ./venv
source ./venv/bin/activate
pip install --upgrade pip
```
###### Tensorflow (2.2.0 version)
```
pip install --upgrade tensorflow==2.2.0
```
###### Tensorflow-addons
```
pip install tensorflow-addons==0.10.0
```
###### OpenCv
```
pip install opencv-python
```

## Use
command to generate the reliefMap:
```
./segunroll -i InputMesh [-h] [CenterlineParameters] [ReliefMapParameters] -o  outputName
```

command to generate the label image :
```
./groundtruth -i InputGrouthTruth  -o outputName
```

command to extract 3D point from segmentation (already computed):
```
./segToMesh -i InputMesh -s segmentationMapName -o outputName
```

command to call segunroll and groundtruth on all mesh file in a directory:
```
./makeAllPair.sh PathToMeshDirectory OutputPath
```

command use our trained model on one reliefmap :
```
python3 predict.py InputReliefMap PathToModel Threshold
```

command to test to complete process of our methode (call segunroll then predict.py then segToMesh):
```
./deep-segmentation.sh PathToModel PathToMesh
```

command to reproduce the result:
```
./testK_folds.sh PathToMeshDirectory PathToModeleDirectory
```
