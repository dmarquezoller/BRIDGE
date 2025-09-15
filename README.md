# BRIDGE

## Overview

BRIDGE is a user-friendly app that enables scientists to **explore, analyze and integrate multi-omics datasets** (proteomics, phosphoproteomics and RNA-seq) interactively, privately and without the need of programming skills. It supports both **individual** and **integrative** analysis of datasets and generates interactive visualizations such as heatmaps, volcano plots, and time-course. BRIDGE is especially powerful for identifying shared biological signals across different omics layers. 

## Installation

In order to download and start using bridge there are some previous steps to be done, like setting the environment and creating the database.

### Setting up the environment

First, the user has to clone the git repository to the local machine.

```bash
git clone https://github.com/dmarquezoller/BRIDGE
 ```

After copying the repository the environment has to be set up in R so all the libraries are available.

```R
renv::restore() #This command is to be done in R after opening the BRIDGE project
```

After this, your local computer will have all the files and required libraries

### Database creation

In order to use the app, a database is needed. Scripts are provided for the user to be guided through the process.
Firstly, a `.db` file has to be created.

```bash
touch user_database.db
```
Then, after creating the empty database, it has to be filled with tables and annotation files, for that, two scripts are provided that will guide the user through the process.

Both scripts assume certain homogeneity in the data. For a correct functioning of the scripts and the app, that is why we put together this set of rulese to be followed both in the manual curation of the data prior to the database creation and in the creation of the database itself.

```bash
python /BRIDGE/Python/db_adding.py
```

```bash
python /BRIDGE/Python/db_adding_annotation.py
```

Both this scripts can be executed as many times as needed.

After all this, the user will have the usable database.

###  Requirements for Your Data

Before submitting tables to the database, please ensure your data follows **all** of the rules below.  
If any rule is not met, the app will most likely crash.  

---

#### <u>1. File format</u>
- The file **must be a CSV**.

#### 2. Identifier columns
- You must provide **at least 3 identifier columns**, with exact names:
  - `Gene_Name` (from gene name)
  - `Gene_ID` (from gene id)
  - `Protein_ID` (from protein id)  

All naming rules must be followed **strictly**.

#### 3. Value columns
- All value (measurement) columns must end with an integer specifying the replicate, **preceded by an underscore** (`_`).  
- **No additional underscores** are allowed in the column name.  
- For extra separation, use other symbols instead.  

Example:  
`X6.hpf_1`

#### 4. Missing values
- **No `NA`s** are allowed in any identifier columns.

#### 5. Table name
- The name of the table must follow this structure: 
`<species>_<datatype>_<optional info>_<id>`

Example:
`zebrafish_proteomics_test_1`

#### 6. Phosphoproteomics data
- An additional identifier is required: **the peptide with the mutation**, named `pepG`.  

Example:
`AAAGDEAGGsSR_p1_ac0`

#### 7. Processed data
- If you are adding processed data, ensure:
  - It is an object of class **`SummarizedExperiment`**.  
  - It was generated using the **same columns** as the raw table (for matching cache keys and tables).

###  Important Notes
- Read these rules carefully and verify that your database meets all requirements.  
- Otherwise, the app may fail to run properly.  

#### Required R packages:
- `storr`
- `DBI`
- `RSQLite`


## Usage

Now, the user already has all the files and libraries as well as the database, so it only remains the execution of the app.

```bash
Rscript /BRIDGE/app.R user_database.db
```

After the execution the app will be open, copy the `url` to your browser and start using BRIDGE!

## File structure

As said before, the code has been heavily modularized to ease the editing, debugging and improvement of the app.
This also allows the user to further locally customize the app with more pipelines or plots without the need of understanding and editing the whole code
but rather just changing the corresponding files.

Here is a diagram showing all the different code files and their hierarchy. Moreover, you can also see which functions are declared in which files and a brief description of what is inside each file.
Following this diagram, the user should find in a fairly easy manner the part of the code they are interested in.

![Code Hierarchy Diagram](./CODE_DIAGRAM_FINAL.png)