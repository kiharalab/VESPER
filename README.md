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

(2) The best superimposition of two maps is identified using the fast Fourier transform (FFT). For each rotation sampled, a translation scan is performed using FFTs to optimize the summation of dot products of matched vectors (DOT score). 

(3) For each of the top 10 models from FFT search, VESPER performs 5° local refinement along each axis and then write top 10 models after the refinement into the output file. 

Given a query map, all other maps in the database are ranked by their similarity to the query, which is measured by the normalized z-score of the best superimposition. To calculate the normalized z-score, we firstly cluster the best DOT score for each of the sampled rotations using single-linkage clustering at cutoff = 20% * (Maximum DOT score – Minimum DOT score). Normalized z-score for the best superimposition is then calculated using mean and standard deviation of the largest cluster. 

## Usage
(1) Identify the best superimpostion between two EM maps by VESPER.
```
Usage: EMVEC_FIT -a [MAP1.mrc (large)] -b [MAP2.mrc (small)] [(option)]

---Options---
-t [float] : Threshold of density map1 def=0.000
-T [float] : Threshold of density map2 def=0.000
-g [float] : Bandwidth of the gaussian filter
             def=16.0, sigma = 0.5*[float]
-s [float] : Sampling grid space def=7.0
-A [float] : Sampling Angle interval def=30.0
-c [int  ] : Number of cores for threads def=2
-N [int  ] : Refine Top [int] models def=10
-S         : Show topN models in PDB format def=false
-V         : Vector Products Mode def=true
-L         : Overlap Mode def=false
-C         : Cross Correlation Coefficient Mode def=false
             Using normalized density Value by Gaussian Filter
-P         : Pearson Correlation Coefficient Mode def=false
             Using normalized density Value by Gaussian Filter and average density
 ```
 
(2) Calculate the normalized z-score for each of top 10 models in VESPER output.
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
make
cp VESPER ../
```
## 2. Identify the best fitting of two EM maps.
```
VESPER -a [MAP1.mrc] -b [MAP2.mrc] [options] > [VESPER_output_filename]
```

**Inputs:**

VESPER expects MAP1.mrc and MAP2.mrc to be valid filenames. Supported file formats are MRC and CCP4.

**Sample Inputs:**

Two sample map files can be found in example_data/ folder.


**Options:**

-a: Name of the first map MAP1. 

-b: Name of the second map MAP2.

-t: Contour level for MAP1. Default = 0.

-T: Contour level for MAP2. Default = 0.

-s: Grid spacing for both maps. We recommend to use 7. Default = 7.0.

-A: Angle spacing for fast Fourier transform (FFT) search. We recommend to use 30. Default = 30.0.

-g: Bandwidth of the Gaussian filter. Default = 16.

-c: Number of cores to use. Default = 2.

-N: The number of top models to perform local refinement. Default = 10.

-S: Show top N models in PDB format. 

-V: Vector products mode. Default = true.

-C: Cross correlation coefficient mode. Density values are not normalized around the mean. Default = false.

-P: Pearson correlation coefficient mode. Default = false.

-L: Overlap mode. Default = false.


**Output format:**
By default, VESPER writes the vector information for each of top 10 models after local refinement into VESPER output. Vector information for the first model starts with two lines like the ones shown below.

```
#0 R={0.0 0.0 5.0} MTX={0.996194700 -0.087155725 0.000000000 0.087155725 0.996194700 0.000000000 -0.000000000 0.000000000 1.000000000} T={37.335 24.247 93.200} sco= 144.563 zsco= 10.477247
Overlap= 0.0871 249/2859 CC= 0.206020 PCC= 0.068408 Score= 144.6
```

In the first line, 0 is the index of the first model. Here the model index starts from 0. MTX shows the rotation matrix of MAP2 relative to MAP1. T shows the translation vector of MAP2 relative to MAP1. sco shows the DOT score, which is the summation of dot products of matched vectors. zsco shows the non-normalized z-score. In the second line, Overlap shows the percentage of overlapped pixels between two maps. CC shows the correlation coefficient where density values are not normalized around the mean. PCC shows the Pearson correlation coefficient. Score is the same as sco, which shows the DOT score.

After these two lines, the output file shows the vector information for the first model. Each vector is represented by two atoms, one for start position (CA) and one for end  (CB). The number in the last column shows the fitness score, which is dot product between this vector and the matched vector in MAP1. Fitness score ranges from -1 to 1: 1 for a perfect match, 0 if two vectors are perpendicular or there is no matched vector, and -1 if two vectors are in the opposite direction. One example is shown below. In this case, there is no matched vector in MAP1 for this specific vector. Thus, fitness score in the last column is 0.00.

```
ATOM      1  CA  ALA     1     101.600 164.600 248.600  1.00  0.00
ATOM      2  CB  ALA     1     106.065 166.626 243.604  1.00  0.00
```


**Usage:**
```
./VESPER ./example_data/emd_8724.map ./example_data/emd_8409.map -t 0.04 -T 0.048 -s 7 -A 30 -c 5 -S > 8724_8409_s7a30.pdb
```

## 3. Calculate normalized z-score for top 10 models in VESPER output.
```
python cluster_score.py -i [VESPER_output_file] -c [Clustering cutoff] -o [Output_file]
```

**Input:**

This program takes the VESPER output file as the input.


**Options:**

-i: Name of input file.

-c: Clustering cutoff (c * (Maximum DOT score – Minimum DOT score)) ranging from 0 to 1. Default = 0.2.

-o: Output filename. By default, output file would be named as input filename followed by .normzscore.


**Output format:**
It shows the normalized z-score for each of top 10 models. One line for each model.


**Usage:**
```
python cluster_score.py -i 8724_8409_s7a30.pdb
```

