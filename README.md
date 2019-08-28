
# XFEM_Fracture2D

### Description

This is a Matlab program that can be used to solve fracture problems involving arbitrary multiple crack propagations in a 2D linear-elastic solid based on the principle of minimum potential energy. The extended finite element method is used to discretise the solid continuum considering cracks as discontinuities in the displacement field. To this end, a strong discontinuity enrichment and a square-root singular crack tip enrichment are used to describe each crack. Several crack growth criteria are available to determine the evolution of cracks over time; apart from the classic maximum tension (or hoop-stress) criterion, the minimum total energy criterion and the local symmetry criterion are implemented implicitly with respect to the discrete time-stepping. 

### Key features 

* *Fast:* The stiffness matrix and the force vector (i.e. the equations' system) and the enrichment tracking data structures are updated at each time step only with respect to the changes in the fracture topology. This ultimately results in the major part of the computational expense in the solution to the linear system of equations rather than in the post-processing of the solution or in the assembly and updating of the equations. As Matlab offers fast and robust direct solvers, the computational times are reasonably fast.

* *Robust.* Suitable for multiple crack propagations with intersections. Furthermore, the stress intensity factors are computed robustly via the interaction integral approach (with the inclusion of the terms to account for crack surface pressure, residual stresses or strains). The minimum total energy criterion and the principle of local symmetry are implemented implicitly in time. The energy release rates are computed based on the stiffness derivative approach using algebraic differentiation (rather than finite differencing of the potential energy). On the other hand, the crack growth direction based on the local symmetry criterion is determined such that the local mode-II stress intensity factor vanishes; the change in a crack tip kink angle is approximated using the ratio of the crack tip stress intensity factors.

* *Easy to run.* Each job has its own input files which are independent form those of all other jobs. The code especially lends itself to running parametric studies. Various results can be saved relating to the fracture geometry, fracture mechanics parameters, and the elastic fields in the solid domain. Extensive visualisation library is available for plotting results.

### Instructions

1. Get started by running the demo to showcase some of the capabilities of the program and to determine if it can be useful for you. At the Matlab's command line enter:

```Matlab
>> RUN_JOBS.m
```

This will execute a series of jobs located inside the *jobs directory* `./JOBS_LIBRARY/`. These jobs do not take very long to execute (around 5 minutes in total).

2. Subsequently, you can pick one of the jobs inside `./JOBS_LIBRARY/` by defining the job title:

```Matlab
>> job_title = 'several_cracks/edge/vertical_tension'
```

3. Then you can open all the relevant scripts for this job as follows:

```Matlab
>> open_job
```

The following input scripts for the *job* will be open in the Matlab's editor:

1. `JOB_MAIN.m`: This is the job's main script. It is called when executing `RUN_JOB` (or `RUN_JOBS`) and acts like a wrapper. Notably, it can serve as a convenient interface to run parametric studies and to save intermediate simulation results.
2. `Input_Scope.m`: This defines the scope of the simulation. From which crack growth criteria to use, to what to compute and what results to show via plots and/or movies. To put it simply, the script is a bunch of "switches" that tell the program what the user wants to be done.
3. `Input_Material.m`: Defines the material's elastic properties in different regions or layers (called "phases") of the computational domain. Moreover, it defines the fracture toughness of the material (assumed to be constant in all material phases). 
4. `Input_Crack.m`: Defines the initial crack geometry.
5. `Input_BC.m`: Defines boundary conditions, such as displacements, tractions, crack surface pressure (assumed to be constant in all cracks), body loads (e.g. gravity, pre-stress or pre-strain).
6. `Mesh_make.m`: In-house structured mesh generator for rectangular domains using either linear triangle or bilinear quadrilateral elements. It is possible to mesh horizontal layers using different mesh sizes.
7. `Mesh_read.m`: Gmsh based mesh reader for version-1 mesh files. Of course you can use your own mesh reader provided the output variables are of the correct format (see later).
8. `Mesh_file.m`: Specifies the mesh input file (.msh). At the moment, only Gmsh mesh files of version-1 are allowed.

### Mesh_file.m

A mesh file needs to be able to output the following data or variables:

* `mNdCrd`: Node coordinates, size = `[nNdStd, 2]`
* `mLNodS`: Element connectivities, size = `[nElemn,nLNodS]`
* `vElPhz`: Element material phase (or region) ID's, size = `[nElemn,1]`
* `cBCNod`: cell of boundary nodes, cell size = `{nBound,1}`, cell element size = `[nBnNod,2]`

Example mesh files are located in `./JOBS_LIBRARY/`. Gmsh version-1 file format is described [here](http://www.manpagez.com/info/gmsh/gmsh-2.4.0/gmsh_60.php).

### Additional notes

* global variables are defined in `.\Routines_AuxInput\Declare_Global.m`
* External libraries are `.\Other_Libs\distmesh` and `.\Other_Libs\mesh2d`

### References

Two external meshing libraries are used for the local mesh refinement and remeshing at the crack tip during crack propagation or prior to a crack intersection with another crack or with a boundary of the domain. Specifically, these libraries, which are located in `.\Other_Libs\`, are the following:

* [*mesh2d*](https://people.sc.fsu.edu/~jburkardt/m_src/mesh2d/mesh2d.html) by Darren Engwirda
* [*distmesh*](http://persson.berkeley.edu/distmesh/) by Per-Olof Persson and Gilbert Strang.  

### Issues and Support

For support or questions please email [sutula.danas@gmail.com](mailto:sutula.danas@gmail.com).

### Authors

Danas Sutula, University of Luxembourg, Luxembourg.


If  you  find  this  code  useful, we  kindly  ask  that  you  consider  citing  us.

* [Minimum energy multiple crack propagation](http://hdl.handle.net/10993/29414)
