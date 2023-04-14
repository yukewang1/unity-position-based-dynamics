using System.Collections;
using System.Collections.Generic;
using UnityEngine;


namespace PBFHelperModules {
    public class SmoothingKernel {
        private float h;  // radius
        private float poly6Multiplier;
        private float spikyMultiplier;
        public float wZero;
        public SmoothingKernel(float radius) {
            this.h = radius;
            this.poly6Multiplier = 315f / (64f * Mathf.PI * Mathf.Pow(h, 9));
            this.spikyMultiplier = -45f / (Mathf.PI * Mathf.Pow(h, 6));
            this.wZero = w(Vector3.zero);
        }

        // The Poly 6 smoothing kernel
        public float w(Vector3 distanceVector) {
            // r is the distance vector
            float r = distanceVector.magnitude;
            if (r > h) { return 0f; }
            return poly6Multiplier * Mathf.Pow(h * h - r * r, 3);
        }

        // Gradient of the spiky kernel
        public Vector3 gradW(Vector3 distanceVector) {
            float r = distanceVector.magnitude;
            if (r > h || r == 0f) { return Vector3.zero; }
            return spikyMultiplier * Mathf.Pow(h - r, 2) * distanceVector.normalized;
        }
    }


    public struct indexedVector3 {
        public int index;
        public Vector3 position;
        public bool isBoundary;
        public indexedVector3(int index, Vector3 position, bool isBoundary = false) {
            this.index = index;
            this.position = position;
            this.isBoundary = isBoundary;
        }
    }

    public class SpatialHashing {
        // The key is for referencing the given query particles array
        // The value refers to stored particles and include information of 
        // whether a particle is a boundary particle
        public Dictionary<int, List<indexedVector3>> neighbors;
        private Dictionary<Vector3Int, List<indexedVector3>> grid;
        private float cellSize;
        private float radius;
        private float radius2;

        public SpatialHashing(float cellSize, float radius) {
            this.cellSize = cellSize;
            this.radius = radius;
            this.radius2 = radius * radius;
            this.grid = new Dictionary<Vector3Int, List<indexedVector3>>();
            this.neighbors = new Dictionary<int, List<indexedVector3>>();
        }

        public void addParticles(
            IEnumerable<Vector3> fluidParticles, IEnumerable<Vector3> boundaryParticles
        ) {
            grid.Clear();
            addParticlesHelper(fluidParticles, false);
            addParticlesHelper(boundaryParticles, true);
        }

        private void addParticlesHelper(IEnumerable<Vector3> particles, bool isBoundary) {
            int i = 0;
            foreach (Vector3 particle in particles) {
                Vector3Int cell = getCell(particle);
                if (!grid.ContainsKey(cell)) {
                    grid.Add(cell, new List<indexedVector3>());
                }
                grid[cell].Add(new indexedVector3(i, particle, isBoundary));
                i++;
            }
        }

        public void getNeighbors(IEnumerable<Vector3> particles, bool boundarySearch) {
            neighbors.Clear();
            int i = 0;
            // Add to neighbors for index to index look up
            foreach (Vector3 p in particles) {
                neighbors.Add(i, new List<indexedVector3>());
                Vector3Int cell = getCell(p);
                // Search surrounding 27 cells
                for (int x = -1; x <= 1; x++) {
                    for (int y = -1; y <= 1; y++) {
                        for (int z = -1; z <= 1; z++) {
                            Vector3Int searchCell = new Vector3Int(cell.x + x, cell.y + y, cell.z + z);
                            if (!grid.ContainsKey(searchCell)) { continue; }
                            // Loop through neighbors, add index to neighbors if within radius
                            foreach (indexedVector3 neighbor in grid[searchCell]) {
                                // boundarySearch means whether the incoming particles are boundary 
                                // particles. We need it to skip the "neighbor" that is the particle
                                // we are seraching its neighbors for itself.
                                if (boundarySearch && neighbor.index == i) { continue; } 
                                if (!boundarySearch && !neighbor.isBoundary && neighbor.index == i) {
                                    continue;
                                }
                                if ((p - neighbor.position).sqrMagnitude <= radius2) {
                                    neighbors[i].Add(neighbor);
                                }
                            }
                        }
                    }
                }
                i++;
            }
        }

        private Vector3Int getCell(Vector3 particle) {
            return new Vector3Int(
                (int)(particle.x / cellSize),
                (int)(particle.y / cellSize),
                (int)(particle.z / cellSize)
            );
        }
    }


    public class Boundary {
        public List<Vector3> boundaryParticles;
        public List<float> psi;
        public float radius;
        public float diameter;
        public float density;
        public float spatialHashCellSize;
        public float spatialHashRadius;
        public SmoothingKernel kernel;
        public Boundary(
            BoundingBox innerBox, BoundingBox outerBox, float radius, float density, 
            float spatialHashCellSize, float spatialHashRadius
        ) {
            this.radius = radius;
            this.diameter = radius * 2;
            this.density = density;
            this.spatialHashCellSize = spatialHashCellSize;
            this.spatialHashRadius = spatialHashRadius;
            this.kernel = new SmoothingKernel(spatialHashRadius);
            this.boundaryParticles = new List<Vector3>();
            this.psi = new List<float>();

            Vector3 numParticles = new Vector3(
                (int)(outerBox.width() / diameter),
                (int)(outerBox.height() / diameter),
                (int)(outerBox.depth() / diameter)
            );

            for (int x = 0; x < numParticles.x; x++) {
                for (int y = 0; y < numParticles.y; y++) {
                    for (int z = 0; z < numParticles.z; z++) {
                        Vector3 particle = new Vector3(
                            outerBox.min.x + x * diameter + radius,
                            outerBox.min.y + y * diameter + radius,
                            outerBox.min.z + z * diameter + radius
                        );
                        if (!innerBox.contains(particle)) {
                            boundaryParticles.Add(particle);
                        }
                    }
                }
            }

            computePsi();
        }

        public void computePsi() {
            SpatialHashing neighborHash = new SpatialHashing(spatialHashCellSize, spatialHashRadius);
            neighborHash.addParticles(new List<Vector3>(), boundaryParticles);
            neighborHash.getNeighbors(boundaryParticles, true);

            for (int i = 0; i < boundaryParticles.Count; i++) {
                List<indexedVector3> neighbors = neighborHash.neighbors[i];
                float accumulativeDensity = kernel.wZero;

                foreach (indexedVector3 neighbor in neighbors) {
                    accumulativeDensity += kernel.w(
                        boundaryParticles[neighbor.index] - boundaryParticles[i]
                    );
                }
                // Combining volume
                psi.Add(density / accumulativeDensity);
            }
        }
    }

    public class Fluid {
        public Vector3[] positions; // fluid particles, not the sphere
        public Vector3[] predictedPositions;
        public Vector3[] velocities;
        public float[] densities;
        public float[] lambdas;  // the Lagrange multiplier
        public float radius;
        public float diameter;
        public BoundingBox fluidBounds;
        public Boundary boundary;
        public SpatialHashing neighborHash;  // Spatial hashing for fluid particles
        public float density;
        public float mass;
        public float viscosity;
        public int iterations;
        // Some initialization to speed up the simulation
        private float initialDensity;

        public Fluid(
            float radius, BoundingBox fluidBounds, Boundary boundary, float density, int iterations,
            float viscosity, float massMultiplier
        ) {
            this.viscosity = viscosity;
            this.radius = radius;
            this.diameter = radius * 2;
            this.iterations = iterations;
            this.mass = massMultiplier * Mathf.Pow(radius, 3) * density;
            this.fluidBounds = fluidBounds;
            this.boundary = boundary;
            this.density = density;
            initParticles();
            this.neighborHash = new SpatialHashing(
                boundary.spatialHashCellSize, boundary.spatialHashRadius
            );
            this.initialDensity = boundary.kernel.wZero * mass;
        }

        public void solve(float dt) {
            if (dt <= 0) { return; }
            Vector3 gravityVelocityDelta = new Vector3(0, -9.81f, 0) * dt;
            for (int i = 0; i < positions.Length; i++) {
                velocities[i] += gravityVelocityDelta;
                predictedPositions[i] = positions[i] + velocities[i] * dt;
            }
            neighborHash.addParticles(predictedPositions, boundary.boundaryParticles);
            neighborHash.getNeighbors(predictedPositions, false);
            Dictionary<int, List<indexedVector3>> neighbors = neighborHash.neighbors;
            solveIncompressibility(neighbors);
            for (int i = 0; i < positions.Length; i++) {
                velocities[i] = (predictedPositions[i] - positions[i]) / dt;
            }
            solveViscosity();
            for (int i = 0; i < positions.Length; i++) {
                positions[i] = predictedPositions[i];
            }
        }

        private void solveIncompressibility(Dictionary<int, List<indexedVector3>> neighbors) {
            for (int i = 0; i < iterations; i++) {
                for (int particleIdx = 0; particleIdx < predictedPositions.Length; particleIdx++) {
                    computeDensities(particleIdx, neighbors[particleIdx]);
                    computeLambdas(particleIdx, neighbors[particleIdx]);
                }
                // DEBUG
                // string arrayAsString = string.Join(", ", densities);
                // Debug.Log(i + ": densities [" + arrayAsString + "]");
                
                // Cannot combine loops because it uses neighbors' lambdas
                for (int particleIdx = 0; particleIdx < predictedPositions.Length; particleIdx++) {
                    shiftPositions(particleIdx, neighbors[particleIdx]);
                }
            }
        }

        private void computeDensities(int idx, List<indexedVector3> neighbors) {
            densities[idx] = initialDensity;
            foreach (indexedVector3 neighbor in neighbors) {
                int neighborIdx = neighbor.index;
                if (neighbor.isBoundary) {
                    densities[idx] += boundary.kernel.w(
                        boundary.boundaryParticles[neighborIdx] - predictedPositions[idx]
                    ) * boundary.psi[neighborIdx];
                } else {
                    densities[idx] += boundary.kernel.w(
                        predictedPositions[neighborIdx] - predictedPositions[idx]
                    ) * mass;
                }
            }
        }

        private void computeLambdas(int idx, List<indexedVector3> neighbors) {
            float c = densities[idx] / density - 1;
            if (c <= 0) {
                lambdas[idx] = 0.0f;
                return;
            }

            float accumulatedSqrMag = 0.0f;
            Vector3 particle = predictedPositions[idx];
            Vector3 accumulatedGrad = Vector3.zero;
            foreach (indexedVector3 neighbor in neighbors) {
                int neighborIdx = neighbor.index;
                if (neighbor.isBoundary) {
                    Vector3 boundaryNeighbor = boundary.boundaryParticles[neighborIdx];
                    Vector3 gradient = -boundary.kernel.gradW(particle - boundaryNeighbor) * 
                                    (boundary.psi[neighborIdx] / density);
                    accumulatedGrad -= gradient;
                    accumulatedSqrMag += gradient.sqrMagnitude;
                } else {
                    Vector3 fluidNeighbor = predictedPositions[neighborIdx];
                    Vector3 gradient = -boundary.kernel.gradW(particle - fluidNeighbor) * 
                                    (mass / density);
                    accumulatedGrad -= gradient;
                    accumulatedSqrMag += gradient.sqrMagnitude;
                }
            }
            accumulatedSqrMag += accumulatedGrad.sqrMagnitude;
            lambdas[idx] = -c / (accumulatedSqrMag + 1e-6f);
        }

        private void shiftPositions(int idx, List<indexedVector3> neighbors) {
            Vector3 deltaPosition = Vector3.zero;
            foreach (indexedVector3 neighbor in neighbors) {
                int neighborIdx = neighbor.index;
                if (neighbor.isBoundary) {
                    deltaPosition += boundary.kernel.gradW(
                        predictedPositions[idx] - boundary.boundaryParticles[neighborIdx]
                    ) * lambdas[idx] * boundary.psi[neighborIdx] / density;
                } else {
                    deltaPosition += boundary.kernel.gradW(
                        predictedPositions[idx] - predictedPositions[neighborIdx]
                    ) * (lambdas[idx] + lambdas[neighborIdx]) * mass / density;
                }
            }
            predictedPositions[idx] += deltaPosition;
        }

        // The XSPH viscosity
        private void solveViscosity() {
            SmoothingKernel kernel = boundary.kernel;
            float viscosityTimesMass = viscosity * mass;
            // Spatial hash contains predicted positions
            // This particle's predicted position = p, its neighbor = q
            for (int i = 0; i < predictedPositions.Length; i++) {
                Vector3 p = predictedPositions[i];
                List<indexedVector3> neighbors = neighborHash.neighbors[i];
                for (int j = 0; j < neighbors.Count; j++) {
                    if (neighbors[j].isBoundary) { continue; }
                    // neighborIndex is fluid particle index now
                    int neighborIndex = neighbors[j].index;
                    Vector3 q = predictedPositions[neighborIndex];
                    float scaleFactor = kernel.w(p - q) * viscosityTimesMass / densities[neighborIndex];
                    velocities[i] += (velocities[neighborIndex] - velocities[i]) * scaleFactor;
                }
            }
        }

        private void initParticles() {
            Vector3 numParticles = new Vector3(
                (int)(fluidBounds.width() / diameter),
                (int)(fluidBounds.height() / diameter),
                (int)(fluidBounds.depth() / diameter)
            );
            int numParticlesTotal = (int)(numParticles.x * numParticles.y * numParticles.z);
            positions = new Vector3[numParticlesTotal];
            predictedPositions = new Vector3[numParticlesTotal];
            velocities = new Vector3[numParticlesTotal];
            densities = new float[numParticlesTotal];
            lambdas = new float[numParticlesTotal];

            for (int x = 0; x < numParticles.x; x++) {
                for (int y = 0; y < numParticles.y; y++) {
                    for (int z = 0; z < numParticles.z; z++) {
                        int idx = (int)(x * numParticles.y * numParticles.z + y * numParticles.z + z);
                        Vector3 particle = new Vector3(
                            fluidBounds.min.x + x * diameter + radius,
                            fluidBounds.min.y + y * diameter + radius,
                            fluidBounds.min.z + z * diameter + radius
                        );
                        positions[idx] = particle;
                    }
                }
            }
        }
    }

    public class BoundingBox {
        public Vector3 min;
        public Vector3 max;

        public BoundingBox(Vector3 min, Vector3 max) {
            this.min = min;
            this.max = max;
            // Debug.Log("Width: " + width() + " Height: " + height() + " Depth: " + depth());
        }

        public float width() { return max.x - min.x; }
        public float height() { return max.y - min.y; }
        public float depth() { return max.z - min.z; }

        public bool contains(Vector3 point) {
            return point.x >= min.x && point.x <= max.x &&
                point.y >= min.y && point.y <= max.y &&
                point.z >= min.z && point.z <= max.z;
        }
    }
}
