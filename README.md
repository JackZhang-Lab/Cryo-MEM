# PART I:  
# Cryo-EM Membrane Subtraction

A Python-based computational tool for membrane signal analysis and subtraction in cryo-EM.

## Overview

This software utilizes 2D averages and their corresponding alignment information, employing methods such as Radon transformation, cross-correlation, L1 norm, Bezier curves, Monte Carlo simulations, and genetic algorithms. It analyzes and subtracts membranes of any shape in cryo-EM images, ultimately producing particle stacks and micrographs with membrane signals removed. This is particular useful for cryo-EM analysis of membrane proteins embeded in the native lipid bilayers with irregular shapes.

## Features

* Capable of analyzing biological membranes of any shape, including simple lines and arcs, as well as more complex shapes like S or W curves;
* Accurately locates and subtracts biological membrane signals;
* Utilizes GPU and CUDA acceleration to enhance computational speed;
* Features a user-friendly GUI for ease of use.

## Requirements

* This software requires a GPU and CUDA acceleration, necessitating the installation of CUDA drivers and libraries, as well as the cupy library.
* In some cases, [pyem](https://github.com/asarnow/pyem) is also needed to convert cryoSPARC’s `.cs` files to Relion’s `.star` format for processing.

## Installation

For specific installation methods, please refer to the installation document: [[Installation Document]](./Mem_subtraction/doc/en/Installation_en.md) [[安装文档]](./Mem_subtraction/doc/zh-CN/Installation_zh-CN.md)

## Usage

For detailed usage tutorials, please refer to the documentation: [[Usage Index]](./Mem_subtraction/doc/index.md)
https://memxterminator.github.io/wiki/

## License

This software is licensed under GPL v3.0.

## Developer

Mr. Zhen Huang at The Zhang Laboratory at Molecular Biophysics and Biochemistry, Yale University.

## Contact

If you have any questions, please contact Zhen Huang at zhen.victor.huang@gmail.com, or Jack Zhang at jack.zhang@yale.edu

===========================================================================

===========================================================================



# #####################################################################################################################
# PART II:  
# Cryo-EM Membrane Modeling

This is a tool to directly sample and calculate membrane curvature from cryo-EM density maps, as well as build lipids bilayers into your membrane density maps.

## Features

With Cryo-MEM Modeler you can:
- Sample your membrane density with PO4 molecules and convert it into pdb files, where you can check or visualize it in your own way.
- Derive 2D curvature contour maps.
- Build PDB models of lipids out a provided cryo-EM density map.
  
