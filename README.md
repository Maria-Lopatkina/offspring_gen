# Population Data Generation from STR Allele Frequencies with PI, PP, and LR Calculations

## Motivation
Understanding the genetic structure and ethnic variability of populations is essential for elucidating 
genetic-demographic processes and advancing DNA identification technologies. Standard forensic markers 
often demonstrate reduced effectiveness in populations with high levels of inbreeding, such as those in 
Dagestan and Siberia.

This project aims to address this limitation by developing a comprehensive test system 
that can offer enhanced performance in these challenging forensic scenarios. By generating and analyzing 
population data, we seek to improve the accuracy of genetic examinations and support more reliable 
identification processes for individuals in these regions.

## Aim and tasks

The **goal** of this project is to develop and implement software tools for generating genetic data and 
calculating key metrics to demonstrate the limitations of current STR marker sets and frequencies. We aimed 
to show that existing frequency databases and marker panels may not always provide accurate or informative 
results. By generating simulated data and calculating **Probability of Paternity** (PP), **Parentage Index**
(PI), and **Likelihood Ratio** (LR) based on real population frequencies, we aim to illustrate that more 
precise and population-specific data lead to more accurate assessments of familial relationships.

Our tasks include:
* **Data Generation**: Developing programs to generate simulated genetic data using realistic allele frequency 
distributions.
* **Metric Calculation**: Creating software to compute PP, PI, and LR from the generated data to evaluate 
parent-child relationships and other familial connections.
* **False Positives Detection**: Implementing tools to identify false positives in "child-potential father" 
comparisons.
* **Sibling and Grandparent Analysis**: Developing programs to calculate LR for "child-potential sibling" and 
"child-potential grandparent" pairs.

## Conclusion
This project developed new tools for generating genetic data and calculating key metrics such as Probability 
of Paternity (PP), Parentage Index (PI), and Likelihood Ratio (LR). Our findings demonstrated that current 
forensic autosomal marker sets may lack informativeness, particularly in highly inbred populations. 
By using real population frequencies, our tools showed improved accuracy in determining familial 
relationships and identifying false positives. This highlights the need for population-specific data to 
enhance forensic genetic testing and improve the reliability of relationship assessments.

## Installation
To set up and run the Python scripts for this project, follow these steps:

### 1. Clone the Repository:

```commandline
git clone https://github.com/Maria-Lopatkina/offspring_gen.git
```

### 2. Navigate to the Project Directory:

```commandline
cd offspring_gen
```

### 3. Create a Virtual Environment (Optional but recommended):

```commandline
python -m venv venv
```

### 4. Activate the Virtual Environment:

+ On Windows:

```commandline
venv\Scripts\activate
```

+ On macOS/Linux:

```commandline
source venv/bin/activate
```

### 5. Install Dependencies:
Ensure you have pip installed. Then run:

```commandline
pip install -r requirements.txt
```

### 6. Configure the Config File:

* Locate the configuration file **config_file.yaml**, open the file and fill in the necessary 
details according to the project's requirements.

### Requirements

* python==3.9.12
* pandas==2.1.0
* numpy==1.25.2
* PyYAML==6.0.1
* et-xmlfile==1.1.0
* openpyxl==3.1.2
* six==1.16.0

## Usage

### 1. Script for genotype data simulation

This script generates a table with the genotypes of individuals based on STR allele frequencies 
and saves it as an .xlsx file.

### Input
```commandline

```

### Run script

```commandline
python3 genotype_simutation.py
```

### Output

```commandline 
```
