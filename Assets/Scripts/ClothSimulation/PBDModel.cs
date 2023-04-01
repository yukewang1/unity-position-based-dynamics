using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using System.Threading.Tasks;
using Unity.Jobs;
using Unity.Burst;
using Unity.Collections;

namespace ClothSimulation
{
    public class PBDModel : MonoBehaviour
    {
        // Density of the cloth.
        public float density = 0.2f;

        // The iteration number to solve constraints.
        public float solverIterations = 4;

        // The coefficient for the cloth to compress.
        public float compressCoefficient = 1f;

        // The coefficient for the cloth to stretch.
        public float stretchCoefficient = 1f;

        // The coefficient for the cloth to bend.
        public float bendCoefficient = 0.2f;

        // The damping coefficient for the cloth.
        public float dampCoefficient = 0.1f;

        // Unit time
        public float deltaT = 0.01f;

        // Coefficient for gravity
        public Vector3 G;

        public WindZone windZone;

        private Vector3 windForce;

        private Mesh mesh;

        // The total number of vertices;
        private int n;

        // Triangle vertices array: every three elements consist a triangle.
        private int[] triangles;
        // Vertice array: store the edge and the triangle information that the vertice involves.
        private Vertice[] vertices;
        // Vertice masses
        private float[] masses;
        // Vertice velocities
        private Vector3[] velocities;
        // Vertice current position
        private Vector3[] verticesPositions;
        // Vertice next position
        private Vector3[] normals;
        // Remain duration time
        private float durationTime = 0;

        // Key is the index of vertice 1 add the index of vertice 2.
        // For example, vertice1 is 0, vertice2 is 121. Key is 000121.
        // Value is the corresponding Edge object.
        // Is designed to verify that no edge is shared by more than two triangles.
        Hashtable edgeTables = new Hashtable();

        private List<DistanceConstraintInfo> distanceConstraints = new List<DistanceConstraintInfo>();

        private List<BendConstraintInfo> bendConstraints = new List<BendConstraintInfo>();

        void Awake()
        {
            // Resize the mesh from the original 121 vertices to 289.
            Mesh originalMesh = GetComponent<MeshFilter>().mesh;
            int originalVertexCount = originalMesh.vertices.Length;
            Vector3[] oldVertices = originalMesh.vertices;
            int l = (int)System.Math.Sqrt((double)originalVertexCount);
            Mesh newMesh = new Mesh();
            int newVertexCount = 289;
            int lNew = (int)System.Math.Sqrt((double)newVertexCount);
            Vector3[] newVertices = new Vector3[newVertexCount];
            int[] newTriangles = new int[(int)System.Math.Pow((double)(lNew - 1), 2) * 6];

            // Construct new positions of mesh vertices.
            float dx = (oldVertices[l - 1].x - oldVertices[0].x) / (lNew - 1);
            float dz = (oldVertices[originalVertexCount - l].z - oldVertices[0].z) / (lNew - 1);

            for (int i = 0; i < lNew; i++)
            {
                for (int j = 0; j < lNew; j++)
                {
                    newVertices[i * lNew + j] = new Vector3(oldVertices[0].x + dx * j, oldVertices[0].y, oldVertices[0].z + dz * i);
                }
            }

            // Construct new triangles of mesh.
            int tmp = 0;
            for (int i = 0; i < lNew - 1; i++)
            {
                for (int j = 0; j < lNew - 1; j++)
                {
                    newTriangles[tmp] = i * lNew + j;
                    newTriangles[tmp + 1] = (i + 1) * lNew + j + 1;
                    newTriangles[tmp + 2] = (i + 1) * lNew + j;
                    newTriangles[tmp + 3] = i * lNew + j;
                    newTriangles[tmp + 4] = i * lNew + j + 1;
                    newTriangles[tmp + 5] = (i + 1) * lNew + j + 1;
                    tmp = tmp + 6;
                }
            }

            // Construct new UV of mesh.
            Vector2[] newUV = new Vector2[newVertexCount];
            float a = 1.0f / (float)(lNew - 1);
            for (int i = 0; i < lNew; i++)
            {
                for (int j = 0; j < lNew; j++)
                {
                    newUV[i * lNew + j] = new Vector2(a * j, a * i);
                }
            }

            // Construct new normals of mesh.
            Vector3[] newNormals = new Vector3[newVertexCount];
            // Update normals
            for (int j = 0; j < newTriangles.Length; j += 3)
            {
                int vertice0 = newTriangles[j];
                int vertice1 = newTriangles[j + 1];
                int vertice2 = newTriangles[j + 2];
                Vector3 p0 = newVertices[vertice0];
                Vector3 p1 = newVertices[vertice1];
                Vector3 p2 = newVertices[vertice2];
                Vector3 n = Vector3.Normalize(Vector3.Cross(p1 - p0, p2 - p0));
                newNormals[vertice0] += n;
                newNormals[vertice1] += n;
                newNormals[vertice2] += n;
            }
            for (int j = 0; j < newVertices.Length; j++)
            {
                newNormals[j] = Vector3.Normalize(newNormals[j]) * -1.0f;
            }

            newMesh.vertices = newVertices;
            newMesh.triangles = newTriangles;
            newMesh.normals = newNormals;
            newMesh.uv = newUV;
            newMesh.RecalculateBounds();
            GetComponent<MeshFilter>().mesh = newMesh;
        }

        // Start is called before the first frame update
        // 1. Init mass, velocity and position of vertices.
        // 2. Build connection among vertices, edges and triangles.
        void Start()
        {
            mesh = GetComponent<MeshFilter>().mesh;
            // Init origin vertice position
            verticesPositions = mesh.vertices;
            n = verticesPositions.Length;
            vertices = new Vertice[n];
            velocities = new Vector3[n];
            normals = new Vector3[n];
            masses = new float[n];
            G = new Vector3(0, -9.81f, 0);

            // Init veclocities -> 0
            for (int i = 0; i < vertices.Length; i++)
            {
                vertices[i] = new Vertice();
                vertices[i].edges = new List<Edge>();
                vertices[i].triangles = new List<int>();
                velocities[i] = new Vector3(0f, 0f, 0f);
                normals[i] = new Vector3(0f, 0f, 0f);
                masses[i] = 0f;
                verticesPositions[i] = verticesPositions[i];
            }

            // Build vertices map with edges and edges map with triangles.
            triangles = mesh.GetTriangles(0);
            for (int i = 0; i < triangles.Length; i += 3)
            {
                int triangleIndex = i / 3;
                int vertice1 = triangles[i];
                int vertice2 = triangles[i + 1];
                int vertice3 = triangles[i + 2];
                int edge1 = getEdgeIndex(vertice1, vertice2);
                int edge2 = getEdgeIndex(vertice2, vertice3);
                int edge3 = getEdgeIndex(vertice1, vertice3);

                // Build DistanceConstraints
                if (!addEdgeCount(edge1))
                {
                    vertices[vertice1].edges.Add(new Edge(vertice1, vertice2));
                    vertices[vertice2].edges.Add(new Edge(vertice1, vertice2));
                    float distance = Vector3.Distance(verticesPositions[vertice1], verticesPositions[vertice2]);
                    distanceConstraints.Add(new DistanceConstraintInfo()
                    {
                        fixedLength = distance,
                        index0 = vertice1,
                        index1 = vertice2
                    });
                }
                if (!addEdgeCount(edge2))
                {
                    vertices[vertice2].edges.Add(new Edge(vertice2, vertice3));
                    vertices[vertice3].edges.Add(new Edge(vertice2, vertice3));
                    float distance = Vector3.Distance(verticesPositions[vertice2], verticesPositions[vertice3]);
                    distanceConstraints.Add(new DistanceConstraintInfo()
                    {
                        fixedLength = distance,
                        index0 = vertice2,
                        index1 = vertice3
                    });
                }
                if (!addEdgeCount(edge3))
                {
                    vertices[vertice1].edges.Add(new Edge(vertice1, vertice3));
                    vertices[vertice3].edges.Add(new Edge(vertice1, vertice3));
                    float distance = Vector3.Distance(verticesPositions[vertice1], verticesPositions[vertice3]);
                    distanceConstraints.Add(new DistanceConstraintInfo()
                    {
                        fixedLength = distance,
                        index0 = vertice1,
                        index1 = vertice3
                    });
                }
                vertices[vertice1].triangles.Add(triangleIndex);
                vertices[vertice2].triangles.Add(triangleIndex);
                vertices[vertice3].triangles.Add(triangleIndex);

                // Init the mass of each vertice.
                // The mass of the triangle will be evenly distributed to the mass of its vertices.
                Vector3 p0 = verticesPositions[vertice1] - verticesPositions[vertice2];
                Vector3 p1 = verticesPositions[vertice1] - verticesPositions[vertice3];
                float tiangleArea = Vector3.Cross(p0, p1).magnitude * 0.5f;
                float triangleMass = tiangleArea * density;
                float contributedMass = triangleMass / 3;
                masses[vertice1] += contributedMass;
                masses[vertice2] += contributedMass;
                masses[vertice3] += contributedMass;

            }

            // Init BendConstraints
            foreach (int key in edgeTables.Keys)
            {
                if ((int)edgeTables[key] != 2)
                {
                    continue;
                }
                int index0 = key / 1000;
                int index1 = key % 1000;
                List<int> triangles0 = vertices[index0].triangles;
                List<int> triangles1 = vertices[index1].triangles;
                List<int> commonTriangle = triangles0.Intersect(triangles1).ToList();
                List<int> triangleIndex = new List<int>();
                foreach (int i in commonTriangle)
                {
                    if (triangles[i*3] != index0 && triangles[i*3] != index1)
                    {
                        triangleIndex.Add(triangles[i * 3]);
                    }
                    if(triangles[i * 3 + 1] != index0 && triangles[i * 3 + 1] != index1)
                    {
                        triangleIndex.Add(triangles[i * 3 + 1]);
                    }
                    if(triangles[i * 3 + 2] != index0 && triangles[i * 3 + 2] != index1)
                    {
                        triangleIndex.Add(triangles[i * 3 + 2]);
                    }
                }
                if (triangleIndex.Count != 2)
                {
                    Debug.Log("Edge has more than two triangles.");
                }
                int index2 = Mathf.Min(triangleIndex[0], triangleIndex[1]);
                int index3 = Mathf.Max(triangleIndex[0], triangleIndex[1]);
                Vector3 p0 = verticesPositions[index0];
                Vector3 p1 = verticesPositions[index1];
                Vector3 p2 = verticesPositions[index2];
                Vector3 p3 = verticesPositions[index3];
                Vector3 n1 = Vector3.Normalize(Vector3.Cross(p1 - p0, p2 - p0));
                Vector3 n2 = Vector3.Normalize(Vector3.Cross(p1 - p0, p3 - p0));
                float fixedAngle = (float)System.Math.Acos((double)Vector3.Dot(n1, n2));
                bendConstraints.Add(new BendConstraintInfo()
                {
                    fixedAngle = fixedAngle,
                    index0 = index0,
                    index1 = index1,
                    index2 = index2,
                    index3 = index3
                });
            }
        }

        // Update is called once per frame
        void Update()
        {
            NativeArray<float> _masses = new NativeArray<float>(vertices.Length, Allocator.Persistent);
            _masses.CopyFrom(masses);
            NativeArray<Vector3> _velocities = new NativeArray<Vector3>(vertices.Length, Allocator.Persistent);
            _velocities.CopyFrom(velocities);
            NativeArray<Vector3> _verticesPositions = new NativeArray<Vector3>(vertices.Length, Allocator.Persistent);
            _verticesPositions.CopyFrom(verticesPositions);
            NativeArray<Vector3> _normals = new NativeArray<Vector3>(vertices.Length, Allocator.Persistent);
            _normals.CopyFrom(normals);
            verticesPositions = mesh.vertices;
            durationTime += Time.deltaTime;
            windForce = windZone.transform.forward * windZone.windMain;
            for (int i = 0; i < vertices.Length; i++)
            {
                normals[i] = new Vector3(0f, 0f, 0f);
                verticesPositions[i] = verticesPositions[i];
            }
            for (int i = 0; i < Mathf.FloorToInt(durationTime / deltaT); i++)
            {
                // Update normals
                for (int j = 0; j < triangles.Length; j += 3)
                {
                    int vertice0 = triangles[j];
                    int vertice1 = triangles[j + 1];
                    int vertice2 = triangles[j + 2];
                    Vector3 p0 = verticesPositions[vertice0];
                    Vector3 p1 = verticesPositions[vertice1];
                    Vector3 p2 = verticesPositions[vertice2];
                    Vector3 n = Vector3.Normalize(Vector3.Cross(p1 - p0, p2 - p0));
                    normals[vertice0] += n;
                    normals[vertice1] += n;
                    normals[vertice2] += n;
                }
                for (int j = 0; j < vertices.Length; j++)
                {
                    normals[j] = Vector3.Normalize(normals[j] * -1.0f);
                }
                NativeArray<Vector3> currrentPositions = new NativeArray<Vector3>(vertices.Length, Allocator.Persistent);
                
                // Predict positions
                prePredictPositions(ref currrentPositions, _masses, _velocities, _verticesPositions, _normals);
                // Solve constraints
                solveConstraints(ref currrentPositions);
                for (int j = 0; j < vertices.Length; j++)
                {
                    velocities[j] = (currrentPositions[j] - verticesPositions[j]) / deltaT;
                    verticesPositions[j] = currrentPositions[j];
                }
                mesh.vertices = verticesPositions;
                currrentPositions.Dispose();
            }
            durationTime = durationTime % deltaT;
            _masses.Dispose();
            _velocities.Dispose();
            _verticesPositions.Dispose();
            _normals.Dispose();
        }

        void LateUpdate()
        {
            mesh.SetVertices(verticesPositions);
            for (int i = 0; i < vertices.Length; i++)
            {
                normals[i] = new Vector3(0f, 0f, 0f);
            }
            for (int j = 0; j < triangles.Length; j += 3)
            {
                int vertice0 = triangles[j];
                int vertice1 = triangles[j + 1];
                int vertice2 = triangles[j + 2];
                Vector3 p0 = verticesPositions[vertice0];
                Vector3 p1 = verticesPositions[vertice1];
                Vector3 p2 = verticesPositions[vertice2];
                Vector3 n = Vector3.Normalize(Vector3.Cross(p1 - p0, p2 - p0));
                normals[vertice0] += n;
                normals[vertice1] += n;
                normals[vertice2] += n;
            }
            for (int j = 0; j < vertices.Length; j++)
            {
                normals[j] = Vector3.Normalize(normals[j]) * -1.0f;
            }
            mesh.SetNormals(normals);
            mesh.RecalculateBounds();
        }

        public void prePredictPositions(ref NativeArray<Vector3> currentPositions, NativeArray<float> _masses, NativeArray<Vector3> _velocities, NativeArray<Vector3> _verticesPositions, NativeArray<Vector3> _normals)
        {
            var jobData = new PredictPosition() { 
                currentPositions = currentPositions,
                windForce = windForce,
                dampCoefficient = dampCoefficient,
                masses = _masses,
                velocities = _velocities,
                verticesPositions = _verticesPositions,
                normals = _normals,
                deltaT = deltaT,
                G = G,
            };
            var handle = jobData.Schedule(vertices.Length, 32);
            handle.Complete();
        }

        public void solveConstraints(ref NativeArray<Vector3> currentPositions)
        {
            Vector3[] deltaPositions = new Vector3[vertices.Length];
            for (int i = 0; i < solverIterations; i++)
            {
                float di = 1.0f / (solverIterations - i);
                double compressIt = 1 - System.Math.Pow(1 - compressCoefficient, 1.0f / solverIterations);
                double stretchIt = 1 - System.Math.Pow(1 - stretchCoefficient, 1.0f / solverIterations);
                double bendIt = 1 - System.Math.Pow(1 - bendCoefficient, 1.0f / solverIterations);
                for (int j = 0; j < vertices.Length; j++)
                {
                    deltaPositions[j] = new Vector3(0f, 0f, 0f);
                }

                // Solve distance constraints
                for (int j = 0; j < distanceConstraints.Count; j++)
                {
                    DistanceConstraintInfo cur = distanceConstraints[j];
                    Vector3 p0 = currentPositions[cur.index0];
                    Vector3 p1 = currentPositions[cur.index1];
                    float m0 = masses[cur.index0];
                    float m1 = masses[cur.index1];
                    Vector3 direction = Vector3.Normalize(p1 - p0);
                    float length = Vector3.Distance(p1, p0);
                    Vector3 deltaP = new Vector3(0f, 0f, 0f);
                    if (length > cur.fixedLength)
                    {
                        deltaP = (float)stretchIt * direction * (length - cur.fixedLength);
                    } else
                    {
                        deltaP = (float)compressIt * direction * (length - cur.fixedLength);
                    }
                    deltaPositions[cur.index0] += deltaP * di * m1 / (m0 + m1);
                    deltaPositions[cur.index1] -= deltaP * di * m0 / (m0 + m1);
                }

                // Solve bending constraints
                for (int j = 0; j < bendConstraints.Count; j++)
                {
                    BendConstraintInfo cur = bendConstraints[j];
                    Vector3 vertice1 = currentPositions[cur.index0];
                    Vector3 vertice2 = currentPositions[cur.index1];
                    Vector3 vertice3 = currentPositions[cur.index2];
                    Vector3 vertice4 = currentPositions[cur.index3];
                    Vector3 p2 = vertice2 - vertice1;
                    Vector3 p3 = vertice3 - vertice1;
                    Vector3 p4 = vertice4 - vertice1;
                    Vector3 n1 = Vector3.Normalize(Vector3.Cross(p2, p3));
                    Vector3 n2 = Vector3.Normalize(Vector3.Cross(p2, p4));
                    float d = Vector3.Dot(n1, n2);
                    Vector3 q3 = (Vector3.Cross(p2, n2) + Vector3.Cross(n1, p2) * d) / Vector3.Cross(p2, p3).magnitude;
                    Vector3 q4 = (Vector3.Cross(p2, n1) + Vector3.Cross(n2, p2) * d) / Vector3.Cross(p2, p4).magnitude;
                    Vector3 q2 = -(Vector3.Cross(p3, n2) + Vector3.Cross(n1, p3) * d) / Vector3.Cross(p2, p3).magnitude
                        - (Vector3.Cross(p4, n1) + Vector3.Cross(n2, p4) * d) / Vector3.Cross(p2, p4).magnitude;
                    Vector3 q1 = -q2 - q3 - q4;
                    float w1 = 1 / masses[cur.index0];
                    float w2 = 1 / masses[cur.index1];
                    float w3 = 1 / masses[cur.index2];
                    float w4 = 1 / masses[cur.index3];
                    float sum = w1 * q1.magnitude * q1.magnitude;
                    sum += w2 * q2.magnitude * q2.magnitude;
                    sum += w3 * q3.magnitude * q3.magnitude;
                    sum += w4 * q4.magnitude * q4.magnitude;
                    sum = System.Math.Max(0.01f, sum);
                    double a = -System.Math.Sqrt(1 - d * d) * (System.Math.Acos((double)d) - cur.fixedAngle) / (double)sum;
                    
                    if (!double.IsNaN(a))
                    {
                        Vector3 deltaP1 = (float)bendIt * w1 * q1 * di * (float)a / sum;
                        Vector3 deltaP2 = (float)bendIt * w2 * q2 * di * (float)a / sum;
                        Vector3 deltaP3 = (float)bendIt * w3 * q3 * di * (float)a / sum;
                        Vector3 deltaP4 = (float)bendIt * w4 * q4 * di * (float)a / sum;
                        deltaPositions[cur.index0] += deltaP1;
                        deltaPositions[cur.index1] += deltaP2;
                        deltaPositions[cur.index2] += deltaP3;
                        deltaPositions[cur.index3] += deltaP4;
                    }
                }

                int k = (int)System.Math.Sqrt((double)n);
                // The position of four vertices of the plane won't change.
                deltaPositions[0] = verticesPositions[0] - currentPositions[0];
                deltaPositions[k - 1] = verticesPositions[k - 1] - currentPositions[k - 1];
                deltaPositions[n - k] = verticesPositions[n - k] - currentPositions[n - k];
                deltaPositions[n - 1] = verticesPositions[n - 1] - currentPositions[n - 1];

                for (int j = 0; j < vertices.Length; j++)
                {
                    currentPositions[j] = currentPositions[j] + deltaPositions[j];
                }
            }
        }

        // Give the edge unique index.
        public int getEdgeIndex(int vertice1, int vertice2)
        {
            if (vertice1 < vertice2)
            {
                return vertice1 * 1000 + vertice2;
            } else
            {
                return vertice2 * 1000 + vertice1;
            }
        }

        // Add the edge to hash table and accelerate the count.
        public bool addEdgeCount(int index)
        {
            bool isExist = false;
            if (edgeTables.Contains(index))
            {
                edgeTables[index] = (int)edgeTables[index] + 1;
                isExist = true;
            }
            else
            {
                edgeTables.Add(index, 1);
            }
            if ((int)edgeTables[index] > 2)
            {
                Debug.Log("More than two triangles share the same Edge.");
            }
            return isExist;
        }
    }
}
