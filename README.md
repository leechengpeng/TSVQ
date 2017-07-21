# TSVQ
[![License](https://img.shields.io/badge/license-BSD-blue.svg)](LICENSE) 

**Tree Structure Vector Quantization (TSVQ)** is a classical quantization technique from signal processing that allows the modeling of probability density functions by the distribution of prototype vectors. It was originally used for data compression. It works by dividing a large set of points (vectors) into groups having approximately the same number of points closet to them. Each group is represented by its centroid point, as in k-means and some other clustering algorithms.

## Usage
1. **Import TSVQ**
```C++
#include "TSVQ.hpp"
```
2. **Build Tree Structures**
```C++
void build(const std::vector<T*>& vVectorSet);
```
3. **Find Code Vector**
```C++
const T* quantizeVector(const T *vVector) const;
```

## Applications
**TSVQ** is used for lossy data compression, lossy data correction, pattern recognition, density estimation and clustering.

### Image Compression
Project: [Applications/Image](https://github.com/leechengpeng/TSVQ/tree/master/Applications/Image)

Requirement: `OpenCV`, `Visual Studio`
#### VQ: 2
![VQ2](Resources/Image/VQ2.jpg)
#### VQ: 8
![VQ8](Resources/Image/VQ8.jpg)
#### VQ: 16
![VQ16](Resources/Image/VQ16.jpg)
#### VQ: 64
![VQ64](Resources/Image/VQ64.jpg)
#### Original Image
![original](Resources/Image/waterfall.jpg)
