# Pathlines of Kármán vortex street
This is a small MATLAB(R) project for the lecture "[Introduction to Numerical Fluid Mechanics](https://www.uni-ulm.de/index.php?id=100164
)".
The aim is the calculation and visualisation of pathlines inside the Kármán vortex street.

<img src="https://snipboard.io/mGSb0d.jpg" width="40%"/>

## Numerical methods for solving
The following solve methods are available in [`pathline.m`](pathline.m):
* Explicit Euler `explicitEuler`
* Implicit Euler `implicitEuler`
* Improved explicit Euler `betterEuler`
* Euler-Heun `eulerHeun`
* Crank-Nicolson `crankNicolson`
* Runge-Kutta 4. order (3/8-rule) `rungeKutta4Newton`
* Runge-Kutta 4. order (classic) `rungeKutta4Classic`

## How to run
Choose one of the following two ways to run the scripts.
### Run [`postprocessing_mini.m`](postprocessing_mini.m)
Running [`postprocessing_mini.m`](postprocessing_mini.m) with no arguments will calculate the pathlines with the following default values:
```
N        = 1000
first    = 100
last     = N-1
subSteps = 100  % number of substeps for the numerical solve method
method   = 'explicitEuler'
```
It is possible to change the default values with named arguments, the order of appearance is irrelevant:
```
postprocessing_mini('method','implicitEuler', ...
                    'first',1, ...
                    'N',300, ...
                    'subSteps',150)
```
The return value of [`postprocessing_mini.m`](postprocessing_mini.m) is optional and will return the mean processing times for the calculation of the pathlines in each time step.

### Run [`evaluatePerfomance.m`](evaluatePerfomance.m)
Running [`evaluatePerfomance.m`](evaluatePerfomance.m) will compare the runtime between the different numerical methods over a changing number of `subSteps`.

<img src="https://snipboard.io/kIsciD.jpg" width="45%"/>
<img src="https://snipboard.io/sDLb19.jpg" width="45%"/>
