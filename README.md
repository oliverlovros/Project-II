#  Van der Waals gas of particle

 Van der Waals gas of particle is a Fortran/Python program that allow us to study the interaction of hard sppheres subjected to a Leenard-Jones potential


<p align="center">
  <img src="https://user-images.githubusercontent.com/66941005/155822626-9a3d667d-cf97-44cb-b0a4-29d1b485c6d4.gif" alt="animated" />
</p>

<div align="center">
 
[![GitHub forks](https://img.shields.io/github/forks/Eines-Informatiques-Avancades/Project-II)](https://github.com/Eines-Informatiques-Avancades/Project-II/network)
[![GitHub issues](https://img.shields.io/github/issues/Eines-Informatiques-Avancades/Project-II)](https://github.com/Eines-Informatiques-Avancades/Project-II/issues)
![GitHub pull requests](https://img.shields.io/github/issues-pr/Eines-Informatiques-Avancades/Project-II)
[![GitHub stars](https://img.shields.io/github/stars/Eines-Informatiques-Avancades/Project-II)](https://github.com/Eines-Informatiques-Avancades/Project-II/stargazers)
![GitHub repo size](https://img.shields.io/github/repo-size/Eines-Informatiques-Avancades/Project-II)
[![GitHub license](https://img.shields.io/github/license/Eines-Informatiques-Avancades/Project-II)](https://github.com/Eines-Informatiques-Avancades/Project-II)
![GitHub language count](https://img.shields.io/github/languages/count/Eines-Informatiques-Avancades/Project-II)

![Fortran](https://img.shields.io/badge/Fortran-%23734F96.svg?style=for-the-badge&logo=fortran&logoColor=white)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
</div>


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#Pre-requirements">Pre-requirements</a></li>
    <li><a href="#Installation">Installation</a></li>
    <li><a href="#Distribution-of-tasks">Distribution of tasks</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- Pre-requirements -->
### Pre-requirements 📋

_Que cosas necesitas para instalar el software y como instalarlas_

```
sudo apt install python3.8
sudo apt-get install gfortran
```
You also need/want to install VMD software
For installation following the steps below: Download it from (http://www.ks.uiuc.edu/Research/vmd/) Then in the terminal type:
```
cd Downloads
  tar -zxvf vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz
```
This will extract the folder to your downloads folder. You can then just
```
cd vmd-1.9.1

./configure LINUXAMD64

cd src

sudo make install
```

<!-- Installation -->
## Installation 🔧

Use  to install the program.

Ex:
```bash
pip install foobar
```
<!-- Usage -->
## Usage ⚙️
In order to run this program, the necessary files are found in the [Working_area](https://github.com/Eines-Informatiques-Avancades/Project-II/tree/master/Working_Area) directory.

This directory contains the following files:

- **Makefile**: Compiles and executes the whole program. It is also designed to manage all the generated files for a more user-friendly experience.
- **parameters.txt**: Data file that contains all the parameters related to the system of study (Number of particles, geometry of the lattice, density, mass,...), the data related to the simulation (Initial temperature, initial distribution, thermostat, integration method,...) and a final section where the names of the output data files are defined.
- **gofr_params.txt**: Data file that contains the necessary parameters in order to run correctly the radial distribution function (g(r)) calculation. **Important**: This file contains information such as the number of particles, the number of configurations, the density,... so it is very important that these parameters coincide with those defined in the parameters.txt file.

To start using the program, the following command has to be used:
```
make all
```
This will compile and execute the program and all the statistic calculations will also be performed. It is **important** to point out that in order to perform the whole program, the user is asked a couple of questions related to the statistical study, so keep an eye on this, otherwise the program will not reach completion.

Once the process has been completed, the unecessary files are cleaned and the output files are sent to their corresponding directories (Results -> Data and Figures).

To visualize the trajectory of the system, VMD software is needed.

<!-- DISTRIBUTION OF TASKS -->
## Distribution of tasks ✒️ 
Project coordinator: Àlex Teruel

- Main program (Fortran program): Working together
- Initial state (module): Oliver Loveday
- Boundary conditions (module): Adrià Calzada
- Forces (module): Daniel Conde
- Integration (module): Àlex Teruel
- Radial distribution function (Fortran program): Àlex Teruel
- Statistics (python program): Daniel Conde
- Visualisation of Results (makefile + gnuplot + results analysis): Joint work

The joint work tasks will be carried out (to a greater extent) by those members who are more advanced in their corresponding tasks.

## Wiki 📖

Puedes encontrar mucho más de cómo utilizar este proyecto en nuestra [Wiki](https://github.com/Eines-Informatiques-Avancades/Project-II/wiki)

<!-- CONTRIBUTING -->
## Contributing 🖇️
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.




<!-- LICENSE -->
## License
[REPOSITORIO](https:...)


