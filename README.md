# PART I:  
# Membrane Subtractor

A Python-based computational tool for membrane signal analysis and subtraction in cryo-EM.

## Overview

This software utilizes 2D averages and their corresponding alignment information, employing methods such as Radon transformation, cross-correlation, L1 norm, Bezier curves, Monte Carlo simulations, and genetic algorithms. It analyzes and subtracts membranes of any shape in cryo-EM, ultimately producing particle stacks and micrographs with membrane signals removed, suitable for subsequent membrane protein analysis.

## Features

* Capable of analyzing biological membranes of any shape, including simple lines and arcs, as well as more complex shapes like S or W curves;
* Accurately locates and subtracts biological membrane signals;
* Utilizes GPU and CUDA acceleration to enhance computational speed;
* Features a user-friendly GUI for ease of use.

## Requirements

* This software requires a GPU and CUDA acceleration, necessitating the installation of CUDA drivers and libraries, as well as the cupy library.
* In some cases, [pyem](https://github.com/asarnow/pyem) is also needed to convert cryoSPARC’s `.cs` files to Relion’s `.star` format for processing.

## Installation

For specific installation methods, please refer to the installation document: [[Installation Document]](./doc/en/Installation_en.md) [[安装文档]](./doc/zh-CN/Installation_zh-CN.md)

## Usage

For detailed usage tutorials, please refer to the documentation: [[Usage Index]](./doc/index.md)

## License

This software is licensed under GPL v3.0.

## Developer

Mr. Zhen Huang at The Zhang Laboratory at Molecular Biophysics and Biochemistry, Yale University.

## Contact

If you have any questions, please contact Zhen Huang at zhen.victor.huang@gmail.com, or Jack Zhang at jack.zhang@yale.edu

# #####################################################################################################################
# PART II:  
