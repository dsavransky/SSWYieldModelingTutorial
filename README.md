# Yield Modeling Tutorial for the Sagan Summer Workshop 2024

This repository contains materials for the yield modeling tutorial and associated activities for the 2024 Sagan Summer Workshop.

## Running the Tutorial

These materials can be run entirely via your browser, with no local installation required.  This is recommended for all beginners and anyone who wants to casually explore the code. Alternatively, you can install the materials locally on your own machine, which allows you to run the materials without an internet connection.

### Method 1: Google Colab (No Local Installation)

Navigate to:

https://colab.research.google.com/github/dsavransky/SSWYieldModelingTutorial/blob/main/Notebooks/SSW2024_YieldModelingTutorial1_Setup.ipynb

Execute each of the code cells in this notebook tagged for Colab setup (&#128992;).  

You will be prompted to allow Colab to access your Google account - you must grant all access permissions requested in order for things to work properly. A new directory called 'SSW2024' (or whatever you change the name of the `ssw_dir` variable to be) will be created in your Google drive, and a copy of the repository holding the tutorial code will be cloned there.


> **NOTE**
> The colab setup only needs to be run once.  If you run it again, the final step will produce an error (because the code directory will already exist in your google dirve).

You can now proceed to the tutorial itself, either via the link at the bottom of the setup notebook ('Open the Tutorial') or by navigating directly to: https://colab.research.google.com/github/dsavransky/SSWYieldModelingTutorial/blob/main/Notebooks/SSW2024_YieldModelingTutorial1.ipynb

If you have previously run the setup notebook, you can go directly to the tutorial link without running the setup notebook again.

> **IMPORTANT**
> If you changed the `ssw_dir` variable in the setup notebook, you must update the `YMT_dir` variable in the second code block of the tutorial to match.

The open-ended projects associated with this tutorial can be accessed at: https://colab.research.google.com/github/dsavransky/SSWYieldModelingTutorial/blob/main/Notebooks/SSW2024_YieldModelingGroupProjects.ipynb

### Method 2: Local Installation

If you wish to run this tutorial without an internet connection, you must install the associated materials locally on your own machine. The setup notebook (linked above in the Google Colab section) will do this for you, or you can do it manually, as described below. 

We recommend using a dedicated virtual python environment.  The instructions below assume that you have python (version 3.10 or higher recommended) and pip installed and working on your machine. For help with that, start here: https://wiki.python.org/moin/BeginnersGuide/. We'll assume that python and pip are runnable as `python` and `pip`, respectively, but they might be called `python3` and `pip3` (or something else) on your particular system. These instructions are based on working in a terminal (macOS/Linux) or command prompt/PowerShell (Windows).

1. Download or clone this repository to your computer (https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
2. Create a python virtual environment (we'll call ours `YieldModelingTutorial` but you can replace this with any name you like). In a terminal/command prompt/powershell/etc, navigate to a directory where you want to create the virtual environment and run:
   
   ```python -m venv YieldModelingTutorial```
   
3. Activate the environment. On macOS/Linux (see https://docs.python.org/3/library/venv.html for Windows details):

    ```source YieldModelingTutorial/bin/activate```

4. In the same terminal with the active virtual environment, navigate to the cloned/downloaded repository.  From the top level directory of the repository (the one that contains the file `setup.py`) run:

    ```pip install --user -e .```
    
    This will install all of the python packages required by the tutorial notebooks.
 
5. You will also need a Jypyter environment to execute the notebooks.  We recommend jupyter-lab:

    ```pip install jupyter-lab```


6. Navigate to the `Notebooks` subdirectory of the repository (this should just be `cd Notebooks` from where you were last) and then start JupyterLab by running `jupyter-lab`

7. Skip the first 5 code blocks of the tutorial notebook (the ones marked as for Collab execution only).

8. To stop JupyterLab, type `ctrl+c` in the terminal where it is running and then hit `ctrl+c` again (or type `y` at the prompt). To deactivate the virtual environment just type `deactivate` at the prompt.  Next time you want to run the tutorials again, you just activate the environment again, navigate to the Notebooks directory and run `jupyter-lab`

>**Warning**
>There appears to be an issue (at least on macOS) where if you already have jupyter-lab installed in a system path, it will be executed initially instead of the one you install in your virtual environment.  A simple fix is to deactivate and re-activate the virtual environment after you run the initial pip installation (i.e., between steps 4 and 5).

