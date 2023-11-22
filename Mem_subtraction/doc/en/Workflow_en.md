This article primarily introduces the basic concept and process of membrane subtraction. For detailed procedures, please see:
* [1-Obtaining 2D averages with the membrane center at the image center](Tutorials/Obtain-2D-Averages_en.md)
* [2-Removing membrane signals using the Radonfit method](Tutorials/Radonfit_tutorials/index.md)
* [3-Removing membrane signals using the Bezierfit method](Tutorial/Bezierfit_tutorials/index.md)

# Basic Concept

1. Use cryoSPARC to select particles with the biological membrane at the image center and their corresponding 2D averages;
2. For each 2D average, use methods like Radon transformation and Bezier curve fitting to derive the functional expression of that type of biological membrane, which includes information about the membrane center, angle, curvature, etc.;
3. Based on the functional expression of the biological membrane and alignment information in cryoSPARC, perform membrane subtraction on all raw images of that type of biological membrane to obtain a particle stack after membrane subtraction;
4. Place the particle stack after membrane subtraction back into the original micrograph for subsequent analysis.

# Basic Process

## 1 Obtaining 2D Averages with the Membrane Center at the Image Center
The goal is to make the analysis of the membrane in the second step more accurate and easier, as the Radon transformation might fail to analyze properly if the membrane center is not at the image center, especially with significant deviations.

### 1.1 Creating a Series of 2D Average Templates with the Membrane Center at the Image Center

Select a 2D average with a good signal-to-noise ratio and dominant biological membrane signals from cryoSPRAC. Based on it, create a series of 2D average templates where the biological membrane's center is at the image center, the direction is vertical, but with varying curvatures.

<center><img src="../img/kappa-templates-image.gif" alt="Kappa templates"></center>

Use these generated templates as new templates to pick particles again, ensuring the accurate positioning of the particles' centers at the cell membrane center.

### 1.2 Extract Particles & 2D Classification
Extract the aforementioned particles for 2D classification. The obtained 2D averages should have their membrane centers at the image center. Use these 2D averages for subsequent analysis.

## 2 Analyzing 2D Averages to Obtain Corresponding Functional Expressions

### 2.1 Radonfit

This method mainly uses Radon transformation and cross-correlation to fit the 2D averages with **simple lines and arcs as models**, obtaining functional expressions that include the membrane center, angle, curvature, etc. It's suitable for **simpler biological membrane models**, such as viral envelopes.

#### 2.1.1 Radon Analysis Blinking

Analyze each 2D average to determine the parameters for Radon transformation, including crop rate, threshold, and the range of Radon transformation angles. Record these parameters into a JSON file for use in the next step of membrane analysis.

#### 2.1.2 Membrane Analysis
Analyze the 2D averages, obtaining information for each 2D average like membrane center, angle, bilayer distance, monolayer thickness (size of sigma after Gaussian fitting), and curvature. This information corresponds to the functional expression, which can be used to obtain the averaged biological membrane and its mask.

### 2.2 Bezierfit

This method primarily uses Bezier curves combined with Monte Carlo and genetic algorithms to fit the biological membrane in the 2D averages with **more complex irregular curves as models**, resulting in several control points and corresponding functional expressions. It's suitable for **more complex biological membrane models**, like S-shaped, W-shaped, etc., such as mitochondrial membranes.

#### 2.2.1 Using the Monte Carlo Method to Randomly Generate Points in the Membrane Area
For each 2D average, first use a maximum value filter to extract the rough biological membrane area from the 2D average, then randomly generate several points in this area as reference points for the initial fitting of the Bezier curve.

#### 2.2.2 Preliminary Fitting Using Genetic Algorithms
Use genetic algorithms for a preliminary fitting of the points generated in the previous step, obtaining several control points and the corresponding functional expressions. This step allows determining the general direction and shape of the membrane based on the position information of the points from the previous step.

#### 2.2.3 Adjusting Control Points Using Genetic Algorithms to Optimize Fitting Results
Further optimize the control points obtained in the preliminary fitting using genetic algorithms again. By maximizing the value of cross-correlation, adjust the positions of the control points to obtain several optimal control points that fit the shape of the membrane in the 2D average. These optimal control points correspond to a functional expression describing the shape of the membrane.

## 3 Membrane

 Subtraction
Each type of 2D average corresponds to several raw particle stacks. From the previous analysis of the 2D average, we can obtain the functional expression of that type of biological membrane, which gives us the membrane information for each raw particle.

Perform membrane subtraction on each previously extracted raw particle to obtain the corresponding subtracted particle, forming another corresponding particle stack. This stack is the subtracted particle stack.

## 4 Placing the Subtracted Particles Back into the Micrograph
Replace the original particles in the micrographs with the corresponding subtracted particles to obtain new micrographs, where the membrane signals can be weakened or even removed.

## 5 Picking Particles in the New Micrograph Using Membrane Protein Templates
In the new subtracted micrographs, use membrane protein templates to pick particles for subsequent calculations.