using System.Collections;
using System.Collections.Generic;
using UnityEngine;

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
}
