# BRIDGE

## Overview

BRIDGE is a. 

## Usage

In order to download and start using bridge there are some previous steps to be done, like setting the environment and creating the database.

### Setting the environment
    
First, the user has to clone the git repository to the local machine.

```bash
git clone git@github.com:dmarquezoller/BRIDGE.git
 ```

After copying the repository the environment has to be set up in R so all the libraries are available.

```R
renv::restore() #This command is to be done inside the project
```

After this, your local computer will have all the files and required libraries

### Database creation

In order to use the app, a database is needed. Scripts are provided for the user to be guided through the process.
Firstly, a `.db` file has to be created.

```bash
touch user_database.db
```
Then, after creating the empty database, it has to be filled with tables and annotation files, for that, two scripts are provided that will guide the user through the process.

```bash
python /BRIDGE/Python/db_adding.py
```

For the data adding script this is the set of rules to be followed:

```
############# REQUIREMENTS OF YOUR DATA #############
#                                                   #
#   1. It must be a csv.                            #
#   2. You require AT LEAST 3 identifier columns:   #
#        - gene name -> renamed -> Gene_Name        #
#        - gene id -> renamed -> Gene_ID            #
#        - protein id -> renamed -> Protein ID      #
#      All of the name rules specified must be      #
#      strictly followed.                           #
#   3. All the timeline columns must have a         #
#      single integer in the end specifying         #
#      the replicate                                #
#        -i.e. X6_hpf_1                             #
#   4. There must be no NAs in any identifier       #
#      columns                                      #
#   5. The name of the table must adhere to the     #
#      following structure:                         #
#        - specie_datatype_(optional info)_id       #
#      i.e. zebrafish_rnaseq_investigatorA_1        #                          
#   6. For phosphoproteomics an additional          #
#      identifier is needed: the peptide with the   #
#      mutation, and it has to be called pepG       #
#      i.e. AAAGDEAGGsSR_p1_ac0                     #                
#                                                   #
#                                                   #
#   * If any of these rules is not met for a        #
#     table submited to the database, the app       #
#     will most likely crash                        #                            
#                                                   #
#####################################################
```


```bash
python /BRIDGE/Python/db_adding_annotation.py
```

For the annotation adding scripts this is the set of rules to be followed:

```
############# REQUIREMENTS OF YOUR DATA #############
#                                                   #
#   1. It must be a csv.                            #
#   2. You require the following columns:           #
#          - Gene_ID                                #
#          - Gene_Name                              #
#          - Chromosome                             #
#          - Gene_Start                             #
#          - Gene_End                               #
#          - Gene_Type                              #
#          - Strand                                 #
#   3. The names of the columns have to be          #
#      the precise ones above.                      #
#   4. You can download this annotation table       #
#      with these columns in biomart.               #
#   5. The name of the annotation table uploaded    #
#      has to follow the strict rule:               #
#      <specie>_annotation                          #
#      i.e. zebrafish annotation                    #
#   6. The name of the specie has to be exactly     #
#      the same as the one used for the data        #
#      table.                                       #
#                                                   #
#                                                   #
#   * If any of these rules is not met for a        #
#     table submited to the database, the app       #
#     will most likely crash                        #                            
#                                                   #
#####################################################
```

Both this scripts can be executed as many times as needed. 

After all this, the user will have the usable database.

### App execution