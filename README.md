# Position-Based Dynamics

Group members: Jingshan Feng (jfa99), Yuke Wang (ywa390)

## Unity Setup

To begin, kindly clone the repository, open it as a Unity project in Unity Hub, and select editor version `2021.3.16f1`, with which this project was developed.

## Repository Structure - Key Directories and Files
- **`Assets/Scripts/`: Holds scripts for both scenes in their own sub-directories.**
    - `ClothSimulation/`: Scripts for the cloth simulation scene.
      - `Constraints.cs`: Contains structures for constraints and the mesh.
      - `PBDModel.cs`: The behavior model for the cloth in which vertices' positions are initialized and updated.
    - `FluidSimulation/`: Scripts for the fluid simulation scene.
        - `HelperModules.cs`: All helper classes used by `PBFDemo`, including a smoothing kernel, a spatial hashing module, a fluid solver definition, etc.
        - `PBFDemo.cs`: Defines the actual fluid model attaching to an object, giving controls from the inspector.
- **`Demos/`: Contains demo videos and gifs for two scenes in their respective sub-directories.**
    - `ClothSimulationDemo/`: Demo videos and gifs for the cloth simulation scene.
    - `FluidSimulationDemo/`: Demo videos are stored directly in this directory, while gifs are stored in a sub-directory.
        - `Gif/`: Gifs for the fluid scene, generated from the demo videos and used in this README for a quick demonstration.
- **`Assets/Scenes/`: Contains the cloth and fluid scenes.**
    - `ClothSimulation.unity`: The cloth simulation scene.
    - `FluidSimulation.unity`: The (base) fluid simulation scene, with fluid particles initially bound to a box at the corner and released upon playing the scene.
        - Additionally, a `FluidIncompressibility.unity` scene was duplicated from the base fluid scene with modified parameters to specifically demonstrate the relationship between the incompressibility and how many iterations we solve it for. 

## Overview

Position-Based Dynamics was first proposed by Müller et al. $^1$ as it provides multiple benefits over traditional methods. Specifically, it is more stable, computationally efficient, less laggy, and would not overshoot or drift. 

In this project, we created two scenes with this idea in Unity. 
- In the first scene, we created a plane to simulate a horizontal cloth with fixed four corners. We applied wind force and gravity to it to observe the deformation of the original plane.
- In the second scene, we created multiple particles to simulate a fluid body, primarily following Macklin et al.'s later research $^2$ for the simulation algorithm. We also partially used Akinci etl.'s SPH research $^3$ for boundary handling and the smoothing kernels (the Poly 6 and the Spiky kernels) from Müller et al. $^4$.

Both scenes involve creating a Unity playground, defining the physical objects, as well as various fine-tuning to identify the best constraint parameters to govern their behaviors. Please kindly refer to our full project report for more details on the intuition, implementation, and results.

## Demos

### The Cloth Simulation Scene

The cloth is positioned horizontally at the initial position. After pressing the play button, the plane will deform like cloth under the force of gravity and wind.

In this demo, 2 parameters can be changed.
- **Wind Force**: By changing the Main property of the WindZone object, we can adjust the wind force in the environment. The demo shows the situation when the wind force was changed from 1 to 5, then to 0. The adjustable range of wind power is 0~5.
    ![image](Demos/ClothSimulationDemo/ChangeWindForce.gif)
- **Shrink Coefficient**: By changing the shrinkCoefficient property of the cloth object, the tightness of clothes can be adjusted. The smaller the shrinkCoefficient is, the looser the cloth is. The demo shows the situation that the shrink coefficient change from 10 to 5. The adjustable range of wind power is 1~10.
    ![image](Demos/ClothSimulationDemo/ChangeShrinkCoefficient.gif)

### The Fluid Simulation Scene

In the fluid scene, please allow me to demonstrate a "baseline" version below, where we tuned viscosity, mass, time step, and other parameters to achieve a stable simulation. A few of these parameters were taken from the original paper $^2$ directly; however, the majority of them were tuned by us to achieve the best results, as we didn't follow the paper's exact implementation.

![image](Demos/FluidSimulationDemo/Gif/Baseline.gif)

Additionally, we present two key comparisons, one for the effect of iterations on how the incompressibility is maintained, and the other for the effect of viscosity on the fluid's behavior.

- **Iterations and Incompressibility**: The density constraint keeps fluid particles from moving too close to each other and maintains the incompressibility of the fluid. As we solve for constraints in iterations in PBD, solving for fewer iterations here means a "less enforced" density constraint, making the fluid somewhat compressible, and vice versa. 

    The demo below shows the difference between 40 iterations (left) and 2 iterations (right), where the fluid is more compressible with 2 iterations. We chose 10 iterations for our baseline as it provides a good balance between incompressibility and efficiency.

    ![image](Demos/FluidSimulationDemo/Gif/IncompressibilityComparison.gif)

- **Viscosity and Fluid Behavior**: As defined in the simulation loop, after solving for constraints, a viscosity term is added to the velocity of each particle to simulate the fluid's viscosity, essentially "smoothing" out the difference in velocities between neighboring particles and making the fluid more cohesive. In contrast, particles tend to undergo more independent motion when the viscosity is set to 0.

    The two demos below show the difference between 0 viscosity (first) and 0.75 viscosity (second). The fluid is more cohesive and sticky with the higher viscosity, though inevitably slightly more fake-looking, as some particles lose their individuality and become more "glued" to their neighbors. We selected 0.2 for our baseline for its balance.

    ![image](Demos/FluidSimulationDemo/Gif/0Viscosity.gif)

    ![image](Demos/FluidSimulationDemo/Gif/0.75Viscosity.gif)

## Work Breakdown

Jingshan Feng: Create the cloth simulation scene. Implement the position-based dynamic algorithm to a plane, recalculate triangles, vertices’ positions, and normals of mesh to make the plane deform like cloth. 

Yuke Wang: Create the fluid simulation scene. Implement the position-based algorithm, including a Poly 6 and a Spiky kernel, a spatial hashing module, a solver module with gravity, the density constraint, and the viscosity term, and a boundary handling module.

## Reference

1. Müller, M., Heidelberger, B., Hennix, M., & Ratcliff, J. (2007). Position based dynamics. *Journal of Visual Communication and Image Representation*, 18(2), 109-118.

    **[Link to the above paper in PDF](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf)**

2. Macklin, M., & Müller, M. (2013). Position based fluids. *ACM Transactions on Graphics (TOG)*, 32(4), 1-12.

    **[Link to the above paper in PDF](https://dl.acm.org/doi/pdf/10.1145/2461912.2461984)**

3. Akinci, N., Ihmsen, M., Akinci, G., Solenthaler, B., & Teschner, M. (2012). Versatile rigid-fluid coupling for incompressible SPH. *ACM Transactions on Graphics (TOG)*, 31(4), 1-8.

    **[Link to the above paper in PDF](https://dl.acm.org/doi/pdf/10.1145/2185520.2185558)**

4. Müller, M., Charypar, D., & Gross, M. (2003, July). Particle-based fluid simulation for interactive applications. In *Proceedings of the 2003 ACM SIGGRAPH/Eurographics symposium on Computer animation* (pp. 154-159).

    **[Link to the above paper in PDF](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=1739fd145ef1d327ab301cacc017af2a87f33086)**
