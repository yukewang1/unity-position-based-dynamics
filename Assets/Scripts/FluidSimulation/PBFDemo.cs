using UnityEngine;
using PBFHelperModules;


public class PBFDemo : MonoBehaviour {
    public float dt;  // Should be same as how often FixedUpdate is called, 0.02f
    public float radius;  // Radius of fluid particles (visualization)
    public float viscosity;
    public float density;
    public float massMultiplier;
    public int iterations;
    // Defines the initial releasing box for fluid particles; Y is height
    public Vector3 fluidBoundsMin;
    public Vector3 fluidBoundsMax;
    public Material fluidMaterial;
    public Vector3 envBoundsMin;
    public Vector3 envBoundsMax;
    public int boundaryLayers;
    private float spatialHashCellSize;
    private float spatialHashRadius; 
    private Fluid fluid;
    private GameObject[] fluidSpheres;

    void Start() {
        this.spatialHashCellSize = 4 * radius;
        this.spatialHashRadius = spatialHashCellSize * spatialHashCellSize;
        initializeBoundaryAndFluid();
    }

    void FixedUpdate() {
        fluid.solve(dt);
        // Update the actual visual fluid particles
        for (int i = 0; i < fluidSpheres.Length; i++) {
            fluidSpheres[i].transform.position = fluid.positions[i];
        }
    }

    private void initializeBoundaryAndFluid() {
        BoundingBox envInnerBounds = new BoundingBox(envBoundsMin, envBoundsMax);
        float boundaryThickness = radius * boundaryLayers * 2;
        Vector3 boundaryDiff = new Vector3(boundaryThickness, boundaryThickness, boundaryThickness);
        BoundingBox envOuterBounds = new BoundingBox(
            envInnerBounds.min - boundaryDiff,
            envInnerBounds.max + boundaryDiff
        );

        Boundary boundary = new Boundary(
            envInnerBounds, envOuterBounds, radius, density, spatialHashCellSize, spatialHashRadius
        );

        // The spheres are actual game objects
        // fluid.positions are the positions of particles, used to update the spheres
        BoundingBox fluidBounds = new BoundingBox(fluidBoundsMin, fluidBoundsMax);
        fluid = new Fluid(
            radius, fluidBounds, boundary, density, iterations, viscosity, massMultiplier
        );

        fluidSpheres = new GameObject[fluid.positions.Length];
        for (int i = 0; i < fluid.positions.Length; i++) {
            Vector3 particle = fluid.positions[i];
            GameObject sphere = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            sphere.transform.parent = transform;
            sphere.transform.position = particle;
            sphere.transform.localScale = new Vector3(
                fluid.diameter, fluid.diameter, fluid.diameter
            );
            sphere.GetComponent<Renderer>().material = fluidMaterial;
            sphere.GetComponent<Collider>().enabled = false;
            fluidSpheres[i] = sphere;
        }
    }
}
