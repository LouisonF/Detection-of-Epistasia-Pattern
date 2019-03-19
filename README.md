# SMMB-ACO and Genetic Algorithm

<img src="http://www.ruepremion.fr/images/logos/Logo-universite-nantes-L.jpg" align="center"
     title="University of Nantes logo">

A C++ project made by François Courtin and Louison Fresnais, under the supervision of Christine Sinoquet.

In this project, we have implemented two espistasis dectection methods. The first one is a Genetic Algorithm and the second one is a Stochastic Multiple Markov Blanket algorithm with Ant Colony Optimization (SMMB-ACO).

## Installation

There is a makefile for each method.

1. SMMB-ACO Installation </br>
    -Open a terminal in the SMMB-ACO folder.</br>
    -Make sure that the makefile is here with :
    ```sh
    $ ls
    ```
    -Then run the following commands in the right order:
    make clean remove previous objects
    ```sh
    $ make clean
    ```
    make install create required directories if they are not already created.
    ```sh
    $ make install
    ```
    make compile all sources files and produce an executable file.
    ```sh
    $ make
    ```
1. Genetic Algorithm Installation </br>
    -Open a terminal in the Genetic folder.</br>
    -Make sure that the makefile is here with :
    ```sh
    $ ls
    ```
    -Then run the following commands in the right order:
    make clean remove previous objects
    ```sh
    $ make clean
    ```
    make compile all sources files and produce an executable file.
    ```sh
    $ make
    ```

## Usage

1. SMMB-ACO

    Once you have a correctly compiled program, you can launch SMMB-ACO with the proper python command.</br>
    There is a launch example in the toy_example folder that can help you.

    The first parameter of the launch_smmb.py file is the path to the data folder.</br>
    The second parameter is the name of the dataset that is going to be used in the evaluation script.</br>
    The thrid parameter is the path to the parameters file.</br>
    The fourth parameter is the number of runs of SMMB-ACO per file.</br>
    The fifth parameter is the number of runs in the dataset.</br>
    The sixth parameter is the size of the epistasis pattern.</br>
    ```sh
    python3 launch_smmb.py path_to_dataset_folder dataset_name path_to_parameters_file number_runs number_snps size_epistasis
    ```
2. Genetic

    The launch procedure for the genetic algorithm is exactly the same as for SMMB-ACO. </br>

    The first parameter of the launch_genet.py file is the path to the data folder.</br>
    The second parameter is the name of the dataset that is going to be used in the evaluation script.</br>
    The thrid parameter is the path to the parameters file.</br>
    The fourth parameter is the number of runs of SMMB-ACO per file.</br>
    The fifth parameter is the number of snps in the dataset</br>
    The sixth parameter is the size of the epistasis pattern.</br>
    ```sh
    python3 launch_smmb.py path_to_dataset_folder dataset_name path_to_parameters_file number_runs number_snps size_epistasis
    ```

## Parameter's file edition

Parameter's file are pretty similar for both methods. The edition of a value in the parameters file is simple: </br>
You just have to modifiy the value after the parameter's name. </br>
Lines starting with # are ignored during the parsing of the file.

As an example, the beginning of SMMB-ACO parameter's file:

<img src="https://nsa40.casimages.com/img/2019/02/13/190213081202204633.png" align="center"
     title="Parameters's file parsing example">


## Results

Both algorithm output their results according to the same folder architecture that the provided data architecture.</br>
A result file from SMMB-ACO will return markov-blankets(constitued of SNPs indexes) with sorted signficant p-value scores (accordingto the alpha error risk provided by the user) and other informations such as the G2 score and its reliability, and the number of occurences for a markov blanket.

As an example, a SMMB-ACO output with a generated dataset with the two last SNPs that are causal SNP's (26-27).

<img src="https://nsa40.casimages.com/img/2019/02/13/190213082731794305.png" align="center"
     title="SMMB-ACO results example">
