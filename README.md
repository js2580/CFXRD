# CFXRD

This library provides Python scripts for X-ray diffraction (XRD) analysis in Polyacrylonitrile (PAN)-based Carbon Fibre Reinforced Polymer, particularly at the semi-crystalline structure of the fibre. This includes Radial intergration analysis (Fibre orientation determination), Azimuthal integration (Lattice spacing determination) and finally lattice strain development under in-situ loading analysis. 

## References: 
J. Srisuriyachot, Carbon fibre lattice strain mapping via microfocus synchrotron X-ray diffraction of a reinforced composite, *Carbon*, 200, pp.347-360, DOI: https://doi.org/10.1016/j.carbon.2022.08.041

# Clone git and create virtual environment

```
git clone https://github.com/js2580/CFXRD.git
```

## Setup and Installation

### Install a Virtual Environment

```
pip install virtualenv
python -m venv CFXRD
```

### Activate environment

Windows:
```
./CFXRD/Script/activate
```
Linux: 
```
source CFXRD/bin/activate
```

### Install necessary libraries and version

```
pip install -r scripts/requirements.txt
```

Run scripts with Spyder

[Set up Python Virtual environment in Spyder](https://medium.com/analytics-vidhya/5-steps-setup-python-virtual-environment-in-spyder-ide-da151bafa337)

# Data reduction (Caking XRD integration)

[Download DAWN science](https://dawnsci.org/)

[Processing steps](https://lightform-group.github.io/wiki/tutorials/sxrd-caking)
