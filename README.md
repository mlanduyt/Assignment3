# ME700-ASSIGNMENT3

## GETTING STARTED

The code overview of the FEA code provided for Assignment 3 part 3 is located under CodeStructure.md

This file is meant as an addition to the finite element analysis held at: https://github.com/Lejeune-Lab-Graduate-Course-Materials/finite-element-analysis/tree/main.

## File Management

1. Git clone 
    ```bash 
        git clone https://github.com/Lejeune-Lab-Graduate-Course-Materials/finite-element-analysis/
    ```

2. Git clone or copy part1.ipynb into tutorial folder. This will execute the additional example to the assignment. 



## Set Up 

Set up the conda environment and test that the code is functioning. Note: The environemnt name for all parts of Assignemnt 1 refers to the bisection method.  

1. Create a new conda environment and activate it.  
    ```bash 
    conda create --name finite-element-analysis-env python=3.12
    ```
    ```bash
    conda activate finite-element-analysis-env
    ``` 
2. Confirm that the right version of Python and pip setuptools are being used. Python should be version 3.12; the second command will update setuptools if necessary.  
    ```bash
    python --version
    ```
    ```bash
    pip install --upgrade pip setuptools wheel
    ```
3. Make sure you are in the ME700-Assignment3 directory.  
4. Install an editable version of the code. Note that the *pyproject.toml* file is required for this.  
    ```bash
    pip install -e .
    ```
5. Install the finite element analysis package from git.
    ```bash
    pip install git+https://github.com/Lejeune-Lab-Graduate-Course-Materials/finite-element-analysis
    ```