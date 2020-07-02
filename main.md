---
title: "AMP for SLOPE"
author: "Bartłomiej Polaczyk"
date: "02 July 2020"
output: 
  html_document:
    css: AMP-SLOPE.css
    keep_md: true
abstract: |
  The aim of this project is to investigate the the approximate message passing algorithm for SLOPE regularization problem based on [@bu2019algorithmic] and compare it with classical convex optimization methods.
  Some numerical experiments regarding the cases that do not fit into the theoretical framework of [@bu2019algorithmic] are also performed and analyzed.
bibliography: AMP-SLOPE.bib
#header-includes: |
---
%LaTeX commands
\newcommand{\R}{\mathbb{R}}




## Theoretical bacground
### Introduction
We are interested in solving the standard linear inverse problem
$$ y = Ax + w, $$
where $y\in \R^n$ is  and $A\in\R^{n\times p}$ are known parameters of the model, $w\in\R^n$ is a random noise vector and $x\in\R^p$ is an unknown vector of paramteres we wish to estimate. 

The problem

### SLOPE

Regulizing function is of the form.
Introduced -- cite Bogdan.

### AMP
@donoho2009message, @bu2019algorithmic, @zdeborova2016statistical, @MR2810285.

### Convex optimization methods

## Numerical experiments

### Comparison with convex optimization methods

### References
