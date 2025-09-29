# Lamian Helper Package

This package provides a convenient interface to use **Lamian**, a tool designed to compare expression changes over pseudotime across different conditions. For detailed information on Lamian, please refer to Professor Wenpin Hou's original repository:  
https://github.com/Winnie09/Lamian

---

## Features

Starting from a Seurat object containing **sample** and **gene cluster (meta.data)** information, this package enables the following workflow:  

1. **Harmony integration & normalization**  
2. **Subsetting & renaming samples or metadata**  
3. **Slingshot trajectory inference analysis**  
4. **Filtering of trajectory and features**  
5. **Lamian parameter configuration**  
6. **Visualization options:**  
   - Heatmap  
   - Population plot  
   - FeaturePlot split by sample  

---

## Important Fixes & Modifications

We have incorporated fixes and improvements in the Lamian codebase to ensure robustness and usability:

1. **Fixing matrix non-invertibility issues**  
2. **Adjusting some default parameters and visualization settings** for better performance and clarity  

These fixes are essential to prevent errors during analysis.

---

## Installation

Follow these steps to set up the Lamian Helper Package in your R environment:

1.  **Create a new R project and activate the `renv` environment.**
    *(This step is typically done within your R IDE, such as RStudio.)*

2.  **Clone the repository into your project directory.**
    Open a terminal in your project's root directory and run:
    ```bash
    git clone git@github.com:Ye1203/Lamian_shinyapp.git
    ```

3.  **Restore the package environment using the provided `renv.lock` file.**
    In your R console, run the following commands:
    ```r
    # Install the 'renv' package if you haven't already
    if (!require("renv")) install.packages("renv")
    
    # Initialize and restore the project environment from the lockfile
    renv::init()
    renv::restore()
    ```

4.  **Configure the application path.**
    In the `app.R` file, change the `path_to_shiny` variable to the correct path where the `Lamian_shinyapp` folder are located.

5.  **Launch the Shiny application.**
    In your R console, run:
    ```r
    # Run the application
    shiny::runApp("Lamian_shinyapp")
    ```

---

## Example Usage

The typical workflow includes loading your Seurat object, running normalization, harmony, slingshot, and then using Lamian for comparative pseudotime analysis with the enhanced visualization options provided.

---

## Contact

For any questions or issues, please feel free to open an issue in this repository or contact Bingtian Ye (btye@bu.edu or biangtian@icloud.com).

---
