# deformable-objects
(brief description)

## Instructions
To compile the code, run make within 2d/projects/simulation_code. After compiling, all code should be run from within the 2D directory.

There are some pre-made commands you can run. For example, executing

```{r}
./makeBunSim
```

will take the bunny.svg file from the data folder, create a triangle mesh, execute quasistatics on the mesh 60 times in order to create a reduction basis, then run a reduced motion simulation on this code. Within the 2D directory, there are three executables like this: ./makeBunSim, ./makeCircleSim, and ./makeCragSim.

To run code without using these pre-made commands, there are three steps. If a .poly, .node, and .ele file already exist for the shape you want to model, you can skip to step 3.

First, a poly file type must be created from the initial svg file. This can be done by running the command

```{r}
./bin/svgToPoly ./data/[yourImageName].svg ./data/[yourImageName].poly
```

After a poly file is created, the specifications for a triangle mesh can be made. This is done by running the command

```{r}
./bin/triangle ./data/[yourImageName].poly
```

After the .poly, .node, and .ele files are created, the simulation can be run. There are two different kinds of simulations that can be run: a quasistatics simulation or a motion simulation. In addition, both of these simulations can be reduced or unreduced.

To run an unreduced motion sim, run the command

```{r}
./bin/sim ./data/[yourImageName].1 MOTION -b [-m]
```

To run the reduced motion sim, run the command

```{r}
./bin/sim ./data/[yourImageName].1 MOTION -b -q [int of your choice] [-m]
```

To run an unreduced quasistatics simulation, run the command

```{r}
./bin/sim ./data/[yourImageName].1 [SQUASH/STRETCH/LSHEAR/RSHEAR] [-m]
```

To run a reduced quasistatics simulation, run the command

...

## Concepts Behind the Code
(describe StVK and Stable Neo-Hookean approaches, describe tensor algebraic approach)

Refer to writeup.pdf within the folder writeup for a more in-depth description of the theory behind this work.
