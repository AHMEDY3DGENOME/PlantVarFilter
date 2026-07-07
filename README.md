# PlantOmicsGwas: An Integrated GWAS and Genomic Prediction Pipeline for Plant Genomes

## Developers & Contributors

| Developer                                 | Expertise                     |
|-------------------------------------------|-------------------------------|
| **Ahmed Yassin**                          | Computational Biologist       |
| **Falak Sher Khan**                       | Computational Biologist       |
| **PlantOmicsGwas Software (Affiliation)** | Ye-Lab (PKU-IAAS)        |
 | **PlantOmicsGwas Status** | Early Demonstration Release (Under Active Development) |

------------------------------------------------------------------------------------------------------------------------
## Developed By: 

- Ahmed Yassin, Computational Biologist || PhD Candidate
- Falak Sher Khan, PhD, computational biologist 
The Peking university of Advanced Agricultural Science (PKU-IAAS), China


> **⚠️ Early Demonstration Release (v0.2.4)**
>
> **PlantOmicsGwas v0.1.0** is the **first public demonstration release** of the software and is intended exclusively for evaluation, testing, and feedback from researchers, collaborators, and the plant genomics community.
>
> The platform is currently under active development. Although the core analytical workflow is functional, we are continuously improving computational performance, large-scale genomic data processing, workflow efficiency, the graphical user interface, and the overall user experience.
>
> Following the publication of our accompanying research manuscript, we will release the **complete production version** of PlantOmicsGwas, including comprehensive documentation, finalized software licensing, detailed installation guides, tutorials, and the full implementation of all planned analytical modules.
>
> We sincerely appreciate feedback from researchers and early adopters, as it will play an important role in the continued refinement and development of the platform.


## Quick Installation (Linux)

PlantOmicsGwas depends on several external bioinformatics tools in addition to the Python package. Please follow the steps below to install all required components.

### Step 1 — Install Miniforge (Conda)

Download and install Miniforge:

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

bash Miniforge3-Linux-x86_64.sh

source ~/.bashrc
```

---

### Step 2 — Create the PlantOmicsGwas Environment

Create a dedicated Conda environment and install all required bioinformatics dependencies.

```bash
mamba create -n plantomicsgwas \
-c conda-forge \
-c bioconda \
python=3.11 \
samtools \
bcftools \
bowtie2 \
minimap2 \
plink \
-y
```

---

### Step 3 — Activate the Environment

```bash
conda activate plantomicsgwas
```

---

### Step 4 — Install PlantOmicsGwas

Install the latest release from **PyPI** together with all optional modules.

```bash
pip install --upgrade pip

pip install plantomicsgwas
```

---

### Step 5 — Launch PlantOmicsGwas

Start the graphical user interface.

```bash
plantomicsgwas-gui
```

Alternatively, launch the command-line interface.

```bash
plantomicsgwas --help
```

---

##  Verify Installation

Verify that all required bioinformatics tools are correctly installed.

```bash
samtools --version
bcftools --version
bowtie2 --version
minimap2 --version
plink --version
```

If all commands return their version numbers without errors, PlantOmicsGwas has been installed successfully and is ready for use.
````
