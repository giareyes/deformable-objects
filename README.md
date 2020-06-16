# deformable-objects
(brief description)

## Instructions
To compile the code, run make within 2d/projects/simulation_code. After compiling, all code should be run from within the 2D directory.

### Pre-made Commands
There are some pre-made commands you can run. For example, executing

```{r}
./makeBunSim
```

will take the bunny.svg file from the data folder, create a triangle mesh, execute quasistatics on the mesh 60 times in order to create a reduction basis, then run a reduced motion simulation on this code. Within the 2D directory, there are three executables like this: ./makeBunSim, ./makeCircleSim, and ./makeCragSim.

### Manually Create Simulation
To run code without using these pre-made commands, there are three steps. If a .poly, .node, and .ele file already exist for the shape you want to model, you can skip to step 3.

#### Step 1
First, a poly file type must be created from the initial svg file. This can be done by running the command

```{r}
./bin/svgToPoly ./data/[yourImageName].svg ./data/[yourImageName].poly
```

#### Step 2
After a poly file is created, the specifications for a triangle mesh can be made. This is done by running the command

```{r}
./bin/triangle ./data/[yourImageName].poly
```

#### Step 3
After the .poly, .node, and .ele files are created, the simulation can be run. There are two different kinds of simulations that can be run: a quasistatics simulation or a motion simulation. In addition, both of these simulations can be reduced or unreduced.

To run an unreduced sim, run the command

```{r}
./bin/sim ./data/[yourImageName].1 [SQUASH/STRETCH/LSHEAR/RSHEAR/MOTION] [-m]
```

To run a reduced sim, run the command

```{r}
./bin/sim ./data/[yourImageName].1 [SQUASH/STRETCH/LSHEAR/RSHEAR/MOTION] -q [int of your choice] [-m]
```

To create a basis, add the flag -n to any command 

## Concepts Behind the Code
(describe StVK and Stable Neo-Hookean approaches, describe tensor algebraic approach)

Refer to writeup.pdf within the folder writeup for a more in-depth description of the theory behind this work.
