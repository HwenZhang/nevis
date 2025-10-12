The repository contains the MATLAB/Python code for simulating the subglacial hydrology with blisters followed by rapid supraglacial lake drainage events.

![Image text](./figures/Figure_1_schematic.png)

# Contents
The repository contains the essential files for running the model and visualising the results. The directory structure is as follows:
```.
├── README.md
├── figures (figures in the manuscript)
├── results
│   └── {casename}
├── src (essential source files)
├── srcgen (Python scripts to generate matlab scripts to do parameter sweep)
└── srcplot (Python/MATLAB scripts for visualisation)
└── generated_scripts (MATLAB scripts that can be run directly)
```

# How to run
1. Clone the repository to your local machine.
2. Genrate the matlab scripts for a given parameter combination using the Python scripts in the `srcgen` folder. You can modify the parameters in the scripts as needed. Otherwise, you can use the provided scripts in `generated_scripts` folder to run the model with the parameters used in the manuscript.
3. The bash script `nevis_run.sh` can be modified to run the generated matlab scripts. Otherwise, you can run the matlab scripts directly in MATLAB.
4. The results will be saved in the `results/{casename}` folder.
5. Use the scripts in the `srcplot` folder to visualise the results. You can modify the scripts as needed.