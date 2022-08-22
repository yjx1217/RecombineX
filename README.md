# RecombineX

<p align="center">
  <img src="https://github.com/yjx1217/RecombineX/blob/master/RecombineX.logo.png" alt="RecombineX_logo" width="547" height="162"/>
</p>

**RecombineX: a computational framework for high-throughput gamete genotyping and meiotic recombination analysis**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description
<div style="text-align: justify"> 
Meiotic recombination is an essential biological process that ensures faithful chromosome segregation and promotes parental allele reshuffling. Tetrad analysis is a powerful approach to quantify the genetic makeups and recombination landscapes of meiotic products. Here we present RecombineX, an integrated computational framework that automates the full workflow of marker identification, gamete genotyping, as well as tetrad- and gamete-based recombination profiling in a high-throughput fashion, capable of processing hundreds of tetrads in a single batch. Aside from conventional reference-based analysis, RecombineX can also perform analysis based on parental genome assemblies, which enables analyzing meiotic recombination landscapes in their native genomic contexts. Additional features such as copy number variation profiling and missing genotype inference further enhance downstream analysis. RecombineX also includes a dedicate module for simulating the genomes and reads of recombinant tetrads for any given organisms, which enables fine-tuned simulation-based hypothesis testing. 
</div>

<p align="center">
  <img src="https://github.com/yjx1217/RecombineX/blob/master/RecombineX.overview.png" alt="RecombineX_overview" width="915" height="888"/>
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
* **04.Gamete_Genotyping_by_Reference_Genome**
  * assigning genotypes to a list of pre-defined gametes based on the reference genome (for the "reference-based" mode only)
* **05.Tetrad_based_Recombination_Profiling_by_Reference_Genome**
  * profiling and classifying recombination events (both COs and GCs) for each tetrad based on the reference genome (for the "reference-based" mode only)
* **06.Gamete_based_Recombination_Profiling_by_Reference_Genome**
  * profiling and classifying recombination events (only COs) for each gamete based on the reference genome (for the "reference-based" mode only)
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
* **16.Gamete_Genotyping_by_Parent_Genomes**
  * assigning genotypes to a list of pre-defined gametes from the same tetrad based on parent genomes (for the "parent-based" mode only)
* **17.Tetrad_based_Recombination_Profiling_by_Parent_Genomes**
  * profiling and classifying recombination events (both COs and GCs) for each tetrad based on parent genomes (for the "parent-based" mode only)
* **18.Gamete_based_Recombination_Profiling_by_Parent_Genomes**
  * profiling and classifying recombination events (only COs) for each gamete based on parent genomes (for the "parent-based" mode only)
* **20.Recombinant_Tetrad_Simulation**
  * simulating recombinant tetrads with defined numbers of COs and NCOs.

## Citation
Jing Li, Bertrand Llorente, Gianni Liti, Jia-Xing Yue. (2022) RecombineX: a computational framework for high-throughput gamete genotyping and tetrad-based meiotic recombination profiling. PLoS Genetics 18(5): e1010047. (https://doi.org/10.1371/journal.pgen.1010047) 

## License
RecombineX itself is distributed under the MIT license but some of its dependencies might have more strict license for commercial use. Please check the licensing details of those dependencies.

## Release history
* v1.1.0 Released on 2022/08/22
  ** new feature highlight: support for random-gamete-based recombination profiling analysis
* v1.0.0 Released on 2022/05/07

## Installation
```sh
git clone https://github.com/yjx1217/RecombineX.git
cd RecombineX
bash ./install_dependencies.sh
```
If the installation succeeds, you should see the following massage:
“RecombineX message: This bash script has been successfully processed! :)”
This signifies the success of the installation process. 

Upon the success of the installation, a subdirectory named build and a file named env.sh will be generated. The build subdirectory holds all the installed dependencies, while the env.sh file contains the execution paths of these dependencies. This file will be automatically loaded to set up the working environment for RecombineX’s various modules. The base directory of RecombineX is defined as $RECOMBINEX_HOME in this file.
If unexpected error occurs during installation, normally you can just re-do “bash ./install_dependencies.sh” step and the installation should be able to automatically resume from the previous interruption point. 


## Requirements
### Hardware, operating system and network requirements
RecombineX is designed for a desktop or computing server running an x86-64-bit Linux operating system. Multithreaded processors are preferred to speed up the process since many modules support multithreaded processing. A stable internet connection is required for its installation. 

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
