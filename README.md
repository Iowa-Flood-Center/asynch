[![Build Status](https://travis-ci.org/Iowa-Flood-Center/asynch.svg?branch=develop)](https://travis-ci.org/Iowa-Flood-Center/asynch)

# ASYNCH

A numerical library for solving differential equations with a tree structure. Emphasis is given to hillslope-link river basin models.

## Introduction

ASYNCH is a software package for solving large systems of ordinary differential equations with a tree structure. At its heart, ASYNCH is a library that
can be mixed with other software to allow users greater fexibility. It is written entirely in the C programming language. ASYNCH currently supports
an interface for programs written in C/C++ or Python.

Although ASYNCH is a library, it does come prepackaged with simple programs for performing basic simulations. These programs can be easily
modified to increase their fexibility. ASYNCH is designed for a distributed memory computer. The Message Passing Interface (MPI) is used for inter-process communication. ASYNCH
also contains routines for downloading and uploading data from and to PostgreSQL databases.

ASYNCH was designed with hydrological applications in mind. Much of the naming conventions and sample models reflect this. However, there is
no reason why ASYNCH cannot be applied to other problems with ordinary differential equations having tree structures.

The novelty of this solver is that it uses an asynchronous integrator to solve the differential equations. Further, the communication between processes occurs asynchronously. Details about how this works can be found
in *Small, et. al. An Asynchronous Solver for Systems of ODEs Linked by a Directed Tree Structure, Advances in Water Resources, 53, March 2013,
23-32*.


## The ASYNCH folder structure

```
+---docs
+---examples
+---forecasters
+---forecasters_js
+---ide
+---m4
+---src
\---tests
```

## Table of contents

 - [Installation](INSTALL.md)
 - [Getting Started](docs/getting_started.md)
 - [Async]()
   - [libasync]()
   - [asynch CLI]()
 - [Forecasters]()
 - [Forecasters JS]()
 - [Examples]()