[<img src="img/refpop.png"/>]()

# refpop

##### Licence, status and metrics
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)]()
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub repo size](https://img.shields.io/github/repo-size/ljacquin/refpop)
![GitHub language count](https://img.shields.io/github/languages/count/ljacquin/refpop)
![GitHub top language](https://img.shields.io/github/languages/top/ljacquin/refpop)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/ljacquin/refpop)  
![GitHub all releases](https://img.shields.io/github/downloads/ljacquin/refpop/total)
![GitHub stars](https://img.shields.io/github/stars/ljacquin/refpop)  

##### Languages and technologies
[![R Badge](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)
[![Shell Script Badge](https://img.shields.io/badge/Shell_Script-121011?style=for-the-badge&logo=gnu-bash&logoColor=white)](https://en.wikipedia.org/wiki/Bash_(Unix_shell))
[![Python Badge](https://img.shields.io/badge/Python-FFD43B?style=for-the-badge&logo=python&logoColor=blue)](https://www.python.org/)
[![Plotly Badge](https://img.shields.io/badge/Plotly-239120?style=for-the-badge&logo=plotly&logoColor=white)](https://plotly.com/)
[![RStudio Badge](https://img.shields.io/badge/RSTUDIOAPI-blue?style=for-the-badge&logoColor=white)](https://www.rstudio.com/)

## Objective, automated tasks & analyses, and instructions :

### 🎯 Objective

This repository hosts ```R``` scripts designed to ensure reproducible data analysis and results, aimed at characterizing the refpop population, according to the FAIR principles defined as follows :

* **F**indable.

* **A**ccessible.

* **I**nteroperable.

* **R**eusable.

### 🔄 Automated tasks & analyses 

#### 📂 Repository overview

The ```refpop/``` repository hosts ```R``` scripts that aim to automate the following tasks and analyses to the best of their ability:

* **Phenotype outlier detection**: Intended to automatically produce a clean phenotype dataset, i.e., one that is free of outliers, using robust multivariate and knowledge-based rule methods.

* **Spatial heterogeneity correction**: Aims to address spatial heterogeneity for each trait in the cleaned phenotype dataset, and then identify environments with low heritability using a robust univariate outlier detection method and exclude these for further analysis. For each trait, heritability distributions across all environments are computed before and after spatial heterogeneity correction; these distributions are also computed by management type.

* **Adjusted LS-means computation**: Focuses on computing adjusted phenotypic LS-means for each genotype across environments, using the adjusted phenotypes obtained from spatial heterogeneity correction in each environment.

* **Adjusted BLUP computation**: Computes the adjusted phenotypic BLUP for each genotype across environments, using the adjusted phenotypes corrected for spatial heterogeneity, and principal component coordinates of genotypes derived from genomic data which are treated as fixed-effect covariates.

* **WISER estimation**: Estimates phenotypes using the WISER approach, which corrects fixed effects for the genetic covariance structure (i.e., population structure) embedded within the experimental design: https://github.com/ljacquin/wiser. 

* **Data structure inference**: Focuses on inferring the genomic, phenotypic, and pedigree data structures of the refpop population.

* **LS-means genomic predictive ability evaluation**: Seeks to evaluate, for each trait, the distributions of genomic predictive abilities associated to several prediction methods, for the adjusted phenotypic LS-means.

* **BLUP genomic predictive ability evaluation**: Seeks to evaluate, for each trait, the distributions of genomic predictive abilities associated to several prediction methods, for the adjusted BLUP phenotypes.

* **Wiser genomic predictive ability evaluation**: Seeks to evaluate, for each trait, the distributions of genomic predictive abilities associated to several prediction methods, for the wiser estimated phenotypes.


* **G x E x M evaluation**: Directed towards evaluating, for each trait, the interactions between genotypes, their environments and management types.

#### └ 📁 Repository details 

The ```refpop/``` repository contains three main folders that are central to the tasks and analyses performed by the ```refpop/``` ```R``` scripts:

* ```src/```: Houses the ```R``` scripts.

* ```data/```: Contains the data sourced, processed and analyzed by the ```R``` scripts. Note that some processed data, such as the clean phenotype dataset for example, are written and stored in this folder. Note that ```data/genotype_data/refpop_genotype.zip``` and ```data/genotype_data/genotype_data.zip``` should be decompressed before using the ```R``` scripts in ```src/```.

* ```results/```: Stores the results generated by the ```R``` script tasks and analyses.

  Furthermore, the ```R``` scripts in ```src/``` are distributed accross several subfolders: ```refpop_data_structure_analysis/```, ```refpop_data_treatment_and_analysis/```, ```refpop_gem_analysis/```, ```refpop_genomic_prediction_and_analysis/``` and ```refpop_wiser_genomic_prediction_and_analysis/```. These folders are organized according to the specific tasks and analyses performed by the user, such as phenotypic outlier identification and removal, spatial heterogeneity correction, phenotype estimation, and G x E x M analyses, as well as more specialized tasks like genomic prediction modeling and analysis, etc.

* #### 🧩 Interdependencies of scripts for analyses
Most scripts rely on the preprocessed data generated by the scripts in ```refpop_data_treatment_and_analysis/```. Therefore, it is strongly recommended to execute the data treatment steps first before proceeding with any further analyses. Otherwise, there are no interdependencies between the ```R``` scripts.

* #### 📜 Details of scripts analyses

  The analyses performed by the different ```R``` scripts in the subfolders of ```src/``` are described as follows:
  
  * ```refpop_data_treatment_and_analysis/```
  
    * ```refpop_0_raw_phenotype_data_manage_type_correction.R```: This script corrects manual annotation errors related to management types in the raw phenotypic dataset.
  
    * ```refpop_1_raw_phenotype_data_outlier_detection_per_env.R```: This script conducts outlier detection for phenotypes using two methods. The first method employs robust multivariate outlier detection, which  considers the covariance structure among phenotypic measurements of traits. This method utilizes the Mahalanobis distance (MD) with the minimum covariance determinant (MCD) estimator to mitigate contamination points (i.e., outliers) during the detection process. Given that MD can be sensitive to the curse of dimensionality - resulting in diminished outlier detection accuracy with an increasing number of variables - a principal component analysis (PCA) is performed beforehand to reduce dimensionality. The second approach is a knowledge-based rule univariate method that applies specific criteria according to the refpop phenotyping protocol:
    
      * An outlier is identified if the sample size exceeds 20. 
      * An outlier is identified if the sample size exceeds the total number of fruits per tree (at harvest date).  
        <p> </p>
  
    * ```refpop_2_raw_phenotype_data_spat_hetero_correct_and_h2_estim.R```: For each trait, this script applies spatial heterogeneity correction to the cleaned phenotype dataset and computes the distributions of heritability values, which are derived from estimates in each environment, both before and after the correction. It also calculates the distributions of heritabilities by management type before and after the correction. The spatial heterogeneity correction is carried out using the spatial analysis of field trials with splines (SpATS). Additionally, this script identifies environments with low heritability using the median absolute deviation (MAD) as a unilateral test and excludes them for further analysis and computations, such as adjusted least squares means (LS-means) of phenotypes or the computed heritability using pooled data from all environments.
    
    * ```refpop_3_adjusted_blups_lsmeans_phenotypes_computation.R```: This script computes the adjusted phenotypic LS-means for each genotype across environments, using the adjusted phenotypes obtained from spatial heterogeneity correction in each environment. The least-squares means (LS-means) are computed after fitting a multiple linear regression model that includes genotype and environment effects as covariates, along with the overall mean. From the adjusted phenotypes obtained from spatial heterogeneity correction in each environment, this script also computes the BLUP for each genotype based on a linear mixed model. In this model, genotype effects are treated as random and assumed to be independent and identically distributed (IID), while environmental effects are treated as fixed. The model also includes an overall mean and principal component coordinates for each genotype, derived from genomic data, which are incorporated as fixed effects.
    
       <p> </p>
       
    Note that the ```R``` script in this subfolder are prefixed with ```refpop_0```, ```refpop_1```, ```refpop_2```, and ```refpop_3```, indicating their sequential execution order.
  
  * ```refpop_data_structure_analysis/```
    * ```refpop_genomic_data_structure_analysis.R```: This script performs structure analyses for the refpop genomic data using both uniform manifold approximation and projection (UMAP) and principal component analysis (PCA).
    * ```refpop_pedigree_and_phenotype_data_structure_analysis.R```: This script performs structure analyses for the refpop pedigree data, phenotype data, and their combination, using both uniform manifold approximation and projection (UMAP) and principal component analysis (PCA).
    * ```refpop_gem_data_structure_analysis.R```: This script is a first and small attempt to perform structure type analyses for G x E x M, which is not conclusive. A regression approach is better suited and recommended for this.
    <p> </p>
        
  * ```refpop_gem_analysis/```
    * ```refpop_gem_analysis.R```: This script estimates the variance components for factors such as genotype, country, year, and management in the refpop, along with their second- and third-order interaction terms. Additionally, it calculates the distributions of several statistics, including the phenotypic mean for each trait across different countries, years, and management types, as well as the correlations between phenotypic values of traits across countries.
    <p> </p>
    
  * ```refpop_genomic_prediction_and_analysis/``` (**superseded** by ```refpop_wiser_genomic_prediction_and_analysis/``` in the workflow, see ```execute_refpop_tasks_and_analyses.sh``` in the main ```refpop/``` repository)
  
    * ```refpop_genomic_prediction_trait.R```: This script evaluates, for a defined trait, the distributions of genomic predictive abilities associated to several prediction methods, for the adjusted phenotypic LS-means. These distributions are evaluated using K-folds cross-validation, for n shuffling scenarios of the refpop population, and using the pearson correlation as a measure of predictive ability between the predicted and adjusted phenotypic LS-means. The implemented prediction methods are random forest (RF), support vector regression (SVR), reproducing kernel hilbert space regression (RKHS), genomic blup (GBLUP) and least absolute shrinkage and selection operator (LASSO).
        <p> </p>

  * ```refpop_wiser_genomic_prediction_and_analysis/```
    * ```refpop_wiser_genomic_prediction_trait.R```: This script evaluates, for a defined trait, the distributions of genomic predictive abilities associated to several prediction methods for wiser estimated phenotypes, the adjusted phenotypic LS-means and BLUPs. These distributions are evaluated using K-folds cross-validation, for n shuffling scenarios of the refpop population, and using the pearson correlation as a measure of predictive ability between : 
      * 1) predicted and wiser estimated phenotypes.
      * 2) predicted and adjusted phenotypic LS-means. 
      * 3) predicted and adjusted phenotypic BLUPs. 
       <p> </p>
       
      The implemented prediction methods are random forest (RF), support vector regression (SVR), reproducing kernel hilbert space regression (RKHS), genomic blup (GBLUP), and least absolute shrinkage and selection operator (LASSO). This script also evaluates the distributions of genomic heritabilities associated to the phenotypes estimated using WISER, LS-means, and BLUP, using the GBLUP approach.


### 💻 Instructions

Download the ```refpop``` repository in the current user's directory on a computing cluster or personal computer using one of the following commands :

  *  ```git clone git@github.com:ljacquin/refpop.git``` <p> </p>
    or
  * ```git clone https://github.com/ljacquin/refpop.git``` 
  <p> </p>
  
  ⚠️ Make sure``` git``` is installed beforehand; if not, install it with ```sudo apt install git```.
  <p> </p>

* Given that ```R ≥ 4.1.2``` is already installed, use the following command to install and test ```refpop``` required ```R``` libraries : 

  * ```R -q --vanilla < src/requirements.R```
  * ```R -q --vanilla < src/test_requirements.R```
  <p> </p>

* Within the ```refpop``` folder, execute the following commands to make scripts and programs executable :

  *  ```chmod u+rwx execute_refpop_tasks_and_analyses.sh```
  *  ```chmod u+rwx src/refpop_genomic_prediction_and_analysis/*.sh```
  <p> </p>

* Finally, execute one of the following commands for executing the refpop tasks and analyses :

  *  ```sbatch execute_refpop_tasks_and_analyses.sh```<p> </p>
    or
  * ```./execute_refpop_tasks_and_analyses.sh``` (i.e., interactive execution)
  <p> </p>

⚠️ The tasks and analyses performed by the ```R``` scripts in the ```refpop``` repository can be run in either ```Unix/Linux``` or ```Windows``` environments, as long as ```R``` and the necessary libraries are installed. For local computations in ```RStudio```, ensure that the ```computation_mode``` variable is set to "local" in the ```R``` scripts located in ```src/```. Indeed, while maintaining the required sequential execution order, each ```R``` script can still be run independently for analyses in ```RStudio``` by setting the ```computation_mode``` variable to "local".

## References

* Jung, Michaela, et al. "The apple REFPOP—a reference population for genomics-assisted breeding in apple." Horticulture research 7 (2020).

* Leys, Christophe, et al. "Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median." Journal of experimental social psychology 49.4 (2013): 764-766.

* Hubert, Mia, and Michiel Debruyne. "Minimum covariance determinant." Wiley interdisciplinary reviews: Computational statistics 2.1 (2010): 36-43.

* Higham, Nicholas J. "Computing the nearest correlation matrix—a problem from finance." IMA journal of Numerical Analysis 22.3 (2002): 329-343.

* Breiman, Leo. "Random forests." Machine learning 45 (2001): 5-32.

* Smola, Alex J., and Bernhard Schölkopf. "A tutorial on support vector regression." Statistics and computing 14 (2004): 199-222.

* Jacquin L, Cao T-V and Ahmadi N (2016) A Unified and Comprehensible View of Parametric and Kernel Methods for Genomic Prediction with Application to Rice. Front. Genet. 7:145. doi: 10.3389/fgene.2016.00145

* Tibshirani, Robert. "Regression shrinkage and selection via the lasso." Journal of the Royal Statistical Society Series B: Statistical Methodology 58.1 (1996): 267-288.

