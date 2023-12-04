# FLAMINGOrLite
Lite version of FLAMINGO, High-resolution 3D chromosome structures reconstruction based on Hi-C. Memory-optimized, scability 

---
## Dependencies
The implementation of the algorithm is based on R. It depends on 3 R packages: strawr, parallel, and Matrix.

---
## Installation of the FLAMINGOrLite package
```
install.packages("devtools")
library(devtools)
install_github('JiaxinYangJX/FLAMINGOrLite',ref='HEAD')
```
---
## Reconstruct the 3D genome structure with FLAMINGOrLite
```
library(FLAMINGOrLite)
res = flamingo_main('../FLAMINGO/4DNFI1UEG1HD.hic','hic',1e6,5e3,'chr21','KR',20)

```

Arguments
hic_data	
Input Hi-C data. Now only support .hic

file_format	
Foramt of the input hic data. Now only support .hic.

domain_res	
Size of the domains in bps, e.g. 1e6. Try strawr::readHicBpResolutions() to see available resolutions.

frag_res	
Size of the fragment in bps, e.g. 5e3. Try strawr::readHicBpResolutions() to see available resolutions.

chr_name	
Name of the chromosome, e.g. chr1. Try strawr::readHicChroms() to see available chromosomes.

normalization	
Normalization method in .hic file. Try strawr::readHicNormTypes() to see available methods. Could be 'NONE'.

nThread	
Number of thread avalable for the reconstruction. Default = 1.

sample_rate	
Fraction of available entries in Hi-C to be used during the reconstruction. Default = 0.75.

lambda	
Weights for all sampled entries. Default = 10.

r	
Weights for distance between consecutive points. Default = 1.

max_dist	
Maximum allowed distance betwee two consecutive points. Default = 0.01

alpha	
Convertion factor between interaction frequency and pairwise distance. Default = -0.25.

inf_dist	
Maximun allowed distance betwee any two points. Default = 2.

error_threshold	
Error thresholds for reconstruction. Default = 1e-3.

max_iter	
Maximum iterations. Default = 500.