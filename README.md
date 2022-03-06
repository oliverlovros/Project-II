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

In order to run the program it is necessary to have the following requirements installed:

The program is written in Fortran language so a compiler like [gfortran](https://gcc.gnu.org/wiki/GFortran) or equivalent is needed.

Python programs are used to perform statistical analysis of the obtained data. Therefore, Python version 3 or higher is required. Also, for python programs to work properly the following libraries are needed:

- [Numpy](https://numpy.org)
- [Matplotlib](https://matplotlib.org)
- [Sys](https://docs.python.org/3/library/sys.html)

To be able to use these libraries, it is recommended to work with anaconda environment.

In case the user wishes to see the trajectory generated in the simulation program, the use of the [VMD](https://www.ks.uiuc.edu/Research/vmd/) program is recommended.

<!-- Installation -->
## Installation 🔧

Use  to install the program.

Ex:
```bash
pip install foobar
```
<!-- Usage -->
## Usage ⚙️
In the Working_area directory
```
#make usage
make modulos
make all

#to visualize...
vmd (explain)
python
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
```
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


