
<!-- README.md is generated from README.Rmd. Please edit that file -->

QDRS
====

<!-- badges: start -->
<!-- badges: end -->

The R package QDRS provides functions for calculating quantitative
disease risk scores, with application to Eletronic Health Records.

Installation
------------

You can install the *development* version of QDRS from Github:

    # install.packages("devtools")
    devtools::install_github("danqingxu/QDRS")

Usage
-----

    library(QDRS)

    # EHR demo data
    # sample.set: matrix of binary features.
    # sample.group: a grouping factor with two levels "Case" and "Control". 
    # training: an index vector that indicates training rows.
    sample.set = EHR$sample.set
    sample.group = EHR$sample.group
    training = EHR$training

    PheRS.res = PheRS(X = sample.set, group = sample.group)
    Eigen.res = eigen.score(X = sample.set, training = training)
    PC.res = PC(X = sample.set, group = sample.group, training = training, pc.num = 1:2)
    LPC.res = LPC(X = sample.set, group = sample.group, training = training)
    LVS.res = LVS.score(X = sample.set)
    NMF1.res = NMF1(X = sample.set)

License
-------

This package is free and open source software, licensed under GPL
(&gt;=2).