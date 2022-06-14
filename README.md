## Table of contents
* [General info](#general-info)
* [Methods-pipeline](#methods-pipeline)
* [Examples](#examples)
* [Contents](#contents)
* [Getting started](#getting-started)

## General info

This Git repository contains python codes for constructing Grassmannian diffusion maps-based polynomial chaos expansion surrogates (GDM PCE), ideal for complex applications and models generating high-dimensional outputs. 

GDM PCE, an encoder-decoder framework, employs Grassmannian diffusion maps, a dimensionallity reduction technique for compressing quantities of interest (QoI), represented by high-dimensional vectors/matrices. Once a latent representation of the data is identified, polynomial chaos expansion is used to construct a map between input stochastic parameters and coordinates on the Grassmannian diffusion manifold. A decoder framework based on geometric harmonics is employed for decompressing generated QoIs back to the ambient space and enable inexpensive Monte Carlo simulations.

The proposed method is able to handle very high-dimensional datsets, perform succesfully in the small data regime and accelarate uncertainty quantification (UQ) tasks in general.

## Methods-pipeline

Details of the methdology can be found in the published paper [here](10.1615/Int.J.UncertaintyQuantification.2022039936).

*Authors: Katiana Kontolati, Dimitrios Loukrezis, Ketson R. M. dos Santos, Dimitrios G. Giovanis, Michael D. Shields*

Below, a **graphical summary** of the method is provided:

<img src="pipeline.png" width="700">

## Application

Three illustrative examples are provided. The first considers a dielectric cylinder suspended in a homogeneous electric field. The second is the classic Lotka-Volterra dynamical system modeling the evolution of two species interacting with each other, one a predator and one a prey. Finally, the third example considers a system of advection-diffusion-reaction equations which models a first-order chemical reaction between two species. 
 
<img src="applications.png" width="900">
 
## Contents

* ```data``` - Contains files with datasets and scripts required to run the examples.

* ```GDM_PCE.py``` - python code to perform the proposed method

* ```Example-1.ipynb``` - Jupyter notebook for example 1

* ```Example-2.ipynb``` - Jupyter notebook for example 2
 
* ```Example-3.ipynb``` - Jupyter notebook for example 3

**Warning:** Depending on the choice of hyperparameters (latent dimension, number of clusters et.c.) the process of constructing the decoder part of the method may be time-consuming. For fast results and testing of the method, it is advised to choose small values for these parameters.

## Getting started

**1.** Create an Anaconda Python 3.8 virtual environment:
```
conda create -n gdm python==3.8  
conda activate gdm
```

**2.** Clone our repo:

```
git clone https://github.com/katiana22/GDM-PCE.git
```

**3.** Install dependencies via the following commands: 

```
cd GDM-PCE  
pip install -r requirements.txt
``` 

## Citation

If you find this GitHub repository useful for your work, please consider citing this work:

```
@article{kontolati2022manifold,
  title={Manifold learning-based polynomial chaos expansions for high-dimensional surrogate models},
  author={Kontolati, Katiana and Loukrezis, Dimitrios and dos Santos, Ketson RM and Giovanis, Dimitrios G and Shields, Michael D},
  journal={International Journal for Uncertainty Quantification},
  volume={12},
  number={4},
  year={2022},
  publisher={Begel House Inc.}
}
```

### Mainteners
[Katiana Kontolati](https://katiana22.github.io/)

:email: : kontolati@jhu.edu



