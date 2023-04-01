using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Jobs;
using Unity.Burst;
using Unity.Collections;

namespace ClothSimulation
{
    public class Vertice
    {
        public List<Edge> edges;
        public List<int> triangles;

        public Vertice()
        {
            edges = null;
            triangles = null;
        }
    }

    public class Edge
    {
        public int vertice1;
        public int vertice2;

        public Edge(int vertice1, int vertice2)
        {
            this.vertice1 = Mathf.Min(vertice1, vertice2);
            this.vertice2 = Mathf.Max(vertice1, vertice2);
        }
    }

    public class Triangles
    {
        public int vertice1;
        public int vertice2;
        public int vertice3;

        public Triangles(int vertice1, int vertice2, int vertice3)
        {
            int[] array = new int[] { vertice1, vertice2, vertice3 };
            System.Array.Sort(array);
            this.vertice1 = array[0];
            this.vertice2 = array[1];
            this.vertice3 = array[2];
        }
    }

    public struct DistanceConstraintInfo
    {
        public int index0;
        public int index1;
        public float fixedLength;
    }

    public struct BendConstraintInfo
    {
        public int index0;
        public int index1;
        public int index2;
        public int index3;
        public float fixedAngle;
    }

    struct PredictPosition : IJobParallelFor
    {
        public NativeArray<Vector3> currentPositions;
        [ReadOnly]
        public Vector3 windForce;
        [ReadOnly]
        public float dampCoefficient;
        [ReadOnly]
        public NativeArray<float> masses;
        // Vertice velocities
        [ReadOnly]
        public NativeArray<Vector3> velocities;
        // Vertice current position
        [ReadOnly]
        public NativeArray<Vector3> verticesPositions;
        // Vertice next position
        [ReadOnly]
        public NativeArray<Vector3> normals;
        // Unit time
        [ReadOnly]
        public float deltaT;

        // Coefficient for gravity
        [ReadOnly]
        public Vector3 G;
        public void Execute(int i)
        {
            currentPositions[i] = new Vector3(0f, 0f, 0f);
            Vector3 p = verticesPositions[i];
            Vector3 v = velocities[i];
            float m = masses[i];
            Vector3 normal = normals[i];
            Vector3 v_next = Vector3.Dot(windForce, normal) * normal * deltaT / m + G * deltaT + v;
            v_next = v_next * Mathf.Max(1 - dampCoefficient * deltaT / m, 0);
            Vector3 p_next = p + v_next * deltaT;
            currentPositions[i] = p_next;
        }
    }
}
