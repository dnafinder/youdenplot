[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/youdenplot)

🌐 Overview
youdenplot.m implements a Youden Plot for inter-laboratory variability analysis when two measurements per laboratory are available. It provides a graphical way to decompose the total error of each laboratory into random and systematic components, highlight potentially gross errors, and assess overall agreement among laboratories. The 2025 refactoring preserves the original computational logic while improving clarity, robustness, and documentation in a modern MATLAB style.

⭐ Features
- Computes the Manhattan median of the two measurement series.
- Decomposes total error into random and systematic components.
- Builds a confidence circle based on the random error and Student’s t distribution.
- Draws the 45-degree precision line through the Manhattan median.
- Adds tangents parallel to the 45-degree line to distinguish systematic from gross errors.
- Classifies each point as inside the circle, between tangents, or outside tangents.
- Returns a structured output with key statistics for further analysis.
- Uses distinguishable_colors.m to assign perceptually distinct colors to laboratories.

🛠️ Installation
1. Clone or download the repository:
   https://github.com/dnafinder/youdenplot
2. Make sure the following files are on the MATLAB path:
   - youdenplot.m
   - distinguishable_colors.m
3. Optionally, open and run the function directly in MATLAB Online using the badge above.

▶️ Usage
Basic usage:
    rng(1)
    data = randn(20,2) + [10 10];
    g = 1:size(data,1);
    Y = youdenplot(data, g);

With custom alpha and non-verbose mode:
    Y = youdenplot(data, g, 0.01, 0);

The function opens a figure with the Youden Plot and optionally prints a classification table to the Command Window, depending on the VERBOSE flag.

🔣 Inputs
- DATA  (required)
  - Type: double matrix, size N×2
  - Description: two measurements per laboratory. Column 1 is the first measurement (sample or run 1) and column 2 is the second measurement (sample or run 2).

- G  (optional)
  - Type: numeric vector, length N
  - Description: grouping variable or laboratory index, used for color coding and identification in the plot.
  - Default: 1:N

- ALPHA  (optional)
  - Type: scalar in (0,1)
  - Description: significance level used to build the confidence circle based on the random error and the t distribution.
  - Default: 0.05

- VERBOSE  (optional)
  - Type: scalar, 0 or 1
  - Description: controls whether a summary table is printed to the Command Window.
  - Default: 1 (table displayed)

📤 Outputs
The function returns a struct Y when an output argument is requested.

Fields of Y:
- Y.m
  - 1×2 vector containing the Manhattan median of the two columns of DATA.
- Y.s
  - Standard deviation of the random error component.
- Y.confidence
  - Confidence level expressed as a percentage, (1 − ALPHA) × 100.
- Y.r
  - Radius of the confidence circle centered at the Manhattan median.
- Y.stats
  - N×6 numeric matrix with columns:
    [Total Random Systematic InCircle InTang OutTang]
    where:
    - Total     = total Euclidean distance from the Manhattan median
    - Random    = random error component
    - Systematic= systematic error component
    - InCircle  = logical flag (1 if inside the circle)
    - InTang    = logical flag (1 if between the tangents)
    - OutTang   = logical flag (1 if outside the tangents)

📘 Interpretation
- Points inside the confidence circle:
  - Mainly affected by random error and considered acceptable under the specified confidence level.
- Points outside the circle but between the tangents:
  - Indicate systematic error (consistent bias) while still aligned with the precision line.
- Points close to the 45-degree precision line:
  - Show high precision (small difference between the two measurements).
- Points outside the tangents:
  - Suggest gross errors or out-of-control laboratory performance.

📝 Notes
- The function requires distinguishable_colors.m to generate maximally perceptually distinct colors for plotting different laboratories or groups.
- The refactoring performed in 2025 does not change the original mathematical logic of the function; it only improves readability, robustness, and interface consistency.
- The Youden Plot is particularly useful in inter-laboratory studies, external quality assessment schemes, and method comparison scenarios involving paired measurements.

📚 Citation
If you use this function in scientific work, please cite it as:

Cardillo G. (2014) Youden’s Plot: compute the Youden’s plot for inter-laboratory variability.  
GitHub repository: https://github.com/dnafinder/youdenplot

👤 Author
Author : Giuseppe Cardillo  
Email  : giuseppe.cardillo.75@gmail.com  
GitHub : https://github.com/dnafinder

⚖️ License
This code is distributed under the MIT License. You are free to use, modify, and redistribute it, provided that the original author and source are properly acknowledged.
