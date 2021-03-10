# RecombineX

<p align="center">
  <img src="https://github.com/yjx1217/RecombineX/blob/master/RecombineX.logo.png" alt="RecombineX logo" width="547" height="162"/>
</p>

**RecombineX: a computational framework for tetrad-based meiotic recombination analysis**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description
<div style="text-align: justify"> 
RecombineX is a computational framework for tetrad-based genotyping and meiotic recombination analysis. It handles the full workflow of marker identification, tetrad genotyping, as well as recombination events profiling and classification and produces publication-quality plots. In addition to the conventional reference-genome-based approach, RecombineX also supports the analysis based on the native parent genomes, therefore permitting the close examination on how native parent genomic backgrounds may affect meiotic recombination landscapes of the resulting tetrads. Moreover, RecombineX can also handle partially viable tetrads (e.g. the tetrad with only 3 viable gametes) with its genotype inference feature, which is very useful for studying genome incompatibility. Also, RecombineX shines in its high scalability, capable of processing thousands of sequenced tetrads. Finally, we also developed a tetrad simulation module for RecombineX, which provides rich parameters for users to simulate recombinant tetrads with all introduced recombination events recorded in detail, which can be very useful for downstream hypothesis testing and software development.
</div>

<p align="center">
  <img src="https://github.com/yjx1217/RecombineX/blob/master/RecombineX.overview.png" alt="RecombineX logo" width="915" height="888"/>
</p>

Under the hood, a series of task-specific modules are provided to carry out the full workflow of RecombineX:

* **00.Reference_Genome**
  * preparing the reference genome
* **00.Parent_Genomes**
  * preparing the genomes of the two native crossing parents (for the "parents-based" mode only)
* **00.Parent_Reads**
  * downloading (by SRA tools) the Illumina reads of the two native crossing parents
* **00.Gamete_Reads**
  * downloading (by SRA tools) the Illumina reads of labeled gametes
* **01.Reference_Genome_Preprocessing**
  * preprocessing the reference genome (for the "reference-based" mode only)
* **02.Polymorphic_Markers_by_Reference_based_Read_Mapping**
  * identifying polymorphic markers between the two crossing parents based on the reference genome (for the "reference-based" mode only)
* **03.Gamete_Read_Mapping_to_Reference_Genome**
  * mapping the reads of labeled gametes to the reference genome (for the "reference-based" mode only)
* **04.Tetrad_Genotyping_by_Reference_Genome**
  * assigning genotypes to labeled gametes from the same tetrad based on the reference genome (for the "reference-based" mode only)
* **05.Recombination_Profiling_by_Reference_Genome**
  * profiling and classifying recombination events for each tetrad based on the reference genome (for the "reference-based" mode only)
* **11.Parent_Genome_Preprocessing**
  * preprocessing the parent genomes (for the "parent-based" mode only)
* **12.Polymorphic_Markers_by_Cross_Parent_Genome_Alignment**
  * identifying polymorphic markers between the two crossing parents based on whole genome alignment of the two parents (for the "parent-based" mode only)
* **13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping**
  * identifying polymorphic markers between the two crossing parents based on cross-parent read mapping (for the "parent-based" mode only)
* **14.Polymorphic_Markers_by_Consensus**
  * identifying consensus polymorphic markers between the two crossing parents based on both whole genome alignment and cross-parent read mapping (for the "parent-based" mode only)
* **15.Gamete_Read_Mapping_to_Parent_Genomes**
  * mapping the reads of labeled gametes to the genomes of two native parents (for the "parent-based" mode only)
* **16.Tetrad_Genotyping_by_Parent_Genomes**
  * assigning genotypes to labeled gametes from the same tetrad based on parent genomes (for the "parent-based" mode only)
* **17.Recombination_Profiling_by_Parent_Genomes**
  * profiling and classifying recombination events for each tetrad based on parent genomes (for the "parent-based" mode only)
* **20.Recombinant_Tetrad_Simulation**
  * simulating recombinant tetrads with defined numbers of COs and NCOs.

## License
RecombineX is distributed under the MIT license.

## Installation
```sh
git clone https://github.com/yjx1217/RecombineX.git
cd RecombineX
bash ./install_dependencies.sh
```

## Requirements
### Hardware, operating system and network requirements
RecombineX is designed for a desktop or computing server running an x86-64-bit Linux operating system. Multithreaded processors are preferred to speed up the process since many modules support multithreaded processing. 

### Software requirements
* bash (https://www.gnu.org/software/bash/)
* bzip2 and libbz2-dev (http://www.bzip.org/)
* curl (https://curl.haxx.se/)
* gcc and g++ (https://gcc.gnu.org/)
* git (https://git-scm.com/)
* GNU make (https://www.gnu.org/software/make/)
* gzip (https://www.gnu.org/software/gzip/)
* libopenssl-devel
* libcurl-devel
* java runtime environment (JRE) v1.8.0 (https://www.java.com)
* perl v5.12 or newer (https://www.perl.org/)
* tar (https://www.gnu.org/software/tar/)
* unzip (http://infozip.sourceforge.net/UnZip.html)
* wget (https://www.gnu.org/software/wget/)
* zlib and zlib-devel (https://zlib.net/)
* xz and xz-devel (https://tukaani.org/xz/)
