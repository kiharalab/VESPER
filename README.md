# VESPER
VESPER is a computational tool using local vector based algorithm that can accurately identify the global and local alignment of cryo-electron microscopy (EM) maps.

Copyright (C) 2019 Xusi Han, Genki Terashi, Daisuke Kihara, and Purdue University.

License: GPL v3 for academic use. (For commercial use, please contact us for different licensing)

Contact: Daisuke Kihara (dkihara@purdue.edu)

## Pre-required software
Python 3: https://www.python.org/downloads/
Numpy: pip/conda install numpy
SciPy: pip/conda install scipy
FFTW: http://www.fftw.org/download.html
Pymol (for visualization): https://pymol.org/2/

## VESPER protocol
VESPER protocol consists of three steps:
(1) To identify the best superimposition of two EM maps, each map is firstly converted to a set of unit vectors using the mean shift algorithm.
(2) After conversion of maps into unit vectors, the best superimposition of two maps is identified using the fast Fourier transform (FFT). For each rotation sampled, a translation scan is performed using FFTs to optimize the summation of dot products of matched vectors (DOT score). 
(3) For each of the top 10 models from FFT search, VESPER would perform 5° local refinement along each axis and then write top 10 models after the refinement into the output file. 

Given a query map, all other maps in the database are ranked by their similarity to the query, which is measured by the normalized z-score of the best superimposition. To calculate the normalized z-score, we firstly cluster the best DOT score for each of the sampled rotations using single-linkage clustering at cutoff = 20% * (Maximum DOT score – Minimum DOT score). Normalized z-score for the best superimposition is then calculated using mean and standard deviation of the largest cluster. 

## Usage
(1) Identify the best superimpostion between two EM maps by VESPER
```
Usage: EMVEC_FIT -a [MAP1.mrc (large)] -b [MAP2.mrc (small)] [(option)]

---Options---
-c [int  ] :Number of cores for threads def=2
-t [float] :Threshold of density map1 def=0.000
-T [float] :Threshold of density map2 def=0.000
-g [float] : bandwidth of the gaussian filter
             def=16.0, sigma = 0.5*[float]
-s [float] : sampling grid space def=2.0
-A [float] : sampling Angle interval def=15.0
-N [int  ] : Refine Top [int] models def=10
-S         : Show topN models in PDB format def=false
-V         : Vector Products Mode def=true
-L         : Overlap Mode def=false
-C         : Cross Correlation Coefficient Mode def=false
             Using normalized density Value by Gaussian Filter
-P         : Pearson Correlation Coefficient Mode def=false
             Using normalized density Value by Gaussian Filter and average density
 ```
 
(2) Calculate the normalized z-score for each of top 10 models in VESPER output
```
python cluster_score.py:
usage: cluster_score.py [-h] -i INPUT_FILE [-c CUTOFF] [-o OUT_NAME]

Calculate the normalized z-score for top 10 models from VESPER. Normalized
z-scores for top 10 models are written into the output file.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Required. Name of input file.
  -c CUTOFF             Optional. Clustering cutoff ranging from 0 to 1.
                        Default = 0.2.
  -o OUT_NAME, --output OUT_NAME
                        Optional. Name of output file. If not specified, the
                        output file would be named as input filename followed
                        by .normzscore.
```                       
                        
## 1. Compile VESPER source code.
Firstly download VESPER code from github.
```
git clone https://github.com/kiharalab/VESPER
```
Next, change into VESPER_code directory and compile from the source codes.
```
cd /your_path_to_VESPER/VESPER_code/
Make
cp VESPER ../
```
## 2. Identify the best fitting of two EM maps.
```
VESPER 
