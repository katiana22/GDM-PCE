## Table of contents
* [General info](#general-info)
* [Methods](#methods)
* [Examples](#examples)
* [Contents](#contents)

## General info

This Git repository contains python codes for constructing Grassmannian diffusion maps-based adaptive polynomial chaos expansion surrogates, ideal for cases when complex computational models generating high-dimensional outputs are considered. The proposed method is able to handle very high-dimensional datsets, perform succesfully in the small data regime and accelarate uncertainty quantification (UQ) tasks in general.

## Methods
GDM-PCE uses concepts from manifold learning 
(nonlinear dimensionallity reduction) and surrogate modeling (Polynomial chaos expansion or PCE) 

Details of the methdology can be found here: 

## Application

Three illustrative examples are provided. The first considers a dielectric cylinder suspended in a homogeneous electric field. The second is the classic Lotka-Volterra dynamical system modeling the evolution of two species interacting with each other, one a predator and one a prey. Finally, the third example considers a system of advection-diffusion-reaction equations which models a first-order chemical reaction between two species. 
 
## Contents

* _data_ - Contains files with data requires to run the examples.

* _GDM-PCE.py_ - python code to perform the proposed method

* _Example-1.py_ - Jupyter notebook for example 1

* _Example-2.py_ - Jupyter notebook for example 2
 
* _Example-3.py_ - Jupyter notebook for example 3


### Contact
For more information or questions please email me at: kontolati@jhu.edu



