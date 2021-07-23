using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Diffusion2D_Library
{
    /// <summary>
    /// Define a vector structre for real-valued data
    /// </summary>
    [Serializable]
    public class RVector : ICloneable
    {
        // Fields
        private readonly int ndim;
        private double[] vector;

        // Constructors
        public RVector(int ndim)
        {
            this.ndim = ndim;
            vector = new double[ndim];
            for (int i = 0; i < ndim; i++)
            {
                vector[i] = 0.0;
            }
        }
        public RVector(double[] cv)
        {
            ndim = cv.Length;
            vector = cv;
        }
        //Indexers
        public double this[int i]
        {
            get
            {
                if (i < 0 || i > ndim)
                {
                    throw new Exception("i is out of range!");
                }
                return vector[i];
            }
            set { vector[i] = value; }
        }
        // Accessors 
        public int GetRVectorSize
        {
            get { return ndim; }
        }
        public RVector Section(int si, int ei)
        {
            if (ei > ndim || si > ndim) { throw new Exception("indices are out of range!"); }
            if (ei < si) { int hold = si; si = ei; ei = hold; }
            int length = ei - si;
            RVector rv = new(length);
            for (int i = si; i < ei; i++) { rv[i - si] = vector[i]; }

            return rv;
        }
        // Override Methods
        public override string ToString()
        {
            string str = "{";
            for (int i = 0; i < ndim - 1; i++)
            {
                str += vector[i] + ", ";
            }
            str += vector[ndim - 1] + "}";
            return str;
        }
        public override bool Equals(object obj)
        {
            return (obj is RVector vector1) && Equals(vector1);
        }
        public bool Equals(RVector cv)
        {
            return vector == cv.vector;
        }
        public override int GetHashCode()
        {
            return vector.GetHashCode();
        }
        public static bool operator ==(RVector v1, RVector v2)
        {
            return v1.Equals(v2);
        }
        public static bool operator !=(RVector v1, RVector v2)
        {
            return !v1.Equals(v2);
        }
        public static RVector operator +(RVector cv)
        {
            return cv;
        }
        public static RVector operator +(RVector v1, RVector v2)
        {
            RVector result = new(v1.GetRVectorSize);
            for (int i = 0; i < v1.GetRVectorSize; i++)
            {
                result[i] = v1[i] + v2[i];
            }
            return result;
        }
        public static RVector operator +(RVector cv, double d)
        {
            RVector result = new(cv.GetRVectorSize);
            for (int i = 0; i < cv.GetRVectorSize; i++)
            {
                result[i] = cv[i] + d;
            }
            return result;
        }
        public static RVector operator +(double d, RVector cv)
        {
            RVector result = new(cv.GetRVectorSize);
            for (int i = 0; i < cv.GetRVectorSize; i++)
            {
                result[i] = cv[i] + d;
            }
            return result;
        }
        public static RVector operator -(RVector cv)
        {
            double[] result = new double[cv.GetRVectorSize];
            for (int i = 0; i < cv.GetRVectorSize; i++)
            {
                result[i] = -cv[i];
            }
            return new RVector(result);
        }
        public static RVector operator -(RVector v1, RVector v2)
        {
            RVector result = new(v1.GetRVectorSize);
            for (int i = 0; i < v1.GetRVectorSize; i++)
            {
                result[i] = v1[i] - v2[i];
            }
            return result;
        }
        public static RVector operator -(RVector cv, double d)
        {
            RVector result = new(cv.GetRVectorSize);
            for (int i = 0; i < cv.GetRVectorSize; i++)
            {
                result[i] = cv[i] - d;
            }
            return result;
        }
        public static RVector operator -(double d, RVector cv)
        {
            RVector result = new(cv.GetRVectorSize);
            for (int i = 0; i < cv.GetRVectorSize; i++)
            {
                result[i] = d - cv[i];
            }
            return result;
        }
        public static RVector operator *(RVector cv, double d)
        {
            RVector result = new(cv.GetRVectorSize);
            for (int i = 0; i < cv.GetRVectorSize; i++)
            {
                result[i] = cv[i] * d;
            }
            return result;
        }
        public static RVector operator *(double d, RVector cv)
        {
            RVector result = new(cv.GetRVectorSize);
            for (int i = 0; i < cv.GetRVectorSize; i++)
            {
                result[i] = d * cv[i];
            }
            return result;
        }
        public static RVector Product(RVector v1, RVector v2)
        {
            RVector result = new(v1.GetRVectorSize);
            for (int i = 0; i < v1.GetRVectorSize; i++)
            {
                result[i] = v1[i] * v2[i];
            }
            return result;
        }
        public static RVector operator /(RVector cv, double d)
        {
            RVector result = new(cv.GetRVectorSize);
            for (int i = 0; i < cv.GetRVectorSize; i++)
            {
                result[i] = cv[i] / d;
            }
            return result;
        }
        public static RVector operator /(double d, RVector cv)
        {
            RVector result = new(cv.GetRVectorSize);
            for (int i = 0; i < cv.GetRVectorSize; i++)
            {
                result[i] = d / cv[i];
            }
            return result;
        }
        // Makes a clone copy of a complex vector
        public RVector Clone()
        {
            RVector cv = new(vector);
            cv.vector = (double[])vector.Clone();
            return cv;
        }
        object ICloneable.Clone()
        {
            return Clone();
        }
        // Methods
        // Calculates the dot product of a complex vector
        public static double DotProduct(RVector v1, RVector v2)
        {
            double result = 0.0;
            for (int i = 0; i < v1.GetRVectorSize; i++)
            {
                result += v1[i] * v2[i];
            }
            return result;
        }
        // Calculates the norm of a complex vector
        public double GetNorm()
        {
            double result = 0.0;
            for (int i = 0; i < GetRVectorSize; i++)
            {
                result += vector[i] * vector[i];
            }
            return Math.Sqrt(result);
        }
        // Calculates the square of the norm of a complex vector
        public double GetNormSquare()
        {
            double result = 0.0;
            for (int i = 0; i < GetRVectorSize; i++)
            {
                result += vector[i] * vector[i];
            }
            return result;
        }
        // Normalizes a complex vector
        public void Normalize()
        {
            double norm = GetNorm();
            if (norm == 0)
            {
                throw new Exception("Normalized a vector with norm of zero!");
            }
            for (int i = 0; i < GetRVectorSize; i++)
            {
                vector[i] /= norm;
            }
        }
        // Calculates the unit vector of a complex vector
        public RVector GetUnitVector()
        {
            RVector result = new(vector);
            result.Normalize();
            return result;
        }
        // Swaps entries in a complex vector
        public RVector SwapRVectorEntries(int m, int n)
        {
            double temp = vector[m];
            vector[m] = vector[n];
            vector[n] = temp;
            return new RVector(vector);
        }
        // Calculates the cross product between two complex vectors 
        public static RVector CrossProduct(RVector v1, RVector v2)
        {
            if (v1.GetRVectorSize != 3)
            {
                throw new Exception("Vector v1 must be 3 dimensional!");
            }
            RVector result = new(3);
            result[0] = v1[1] * v2[2] - v1[2] * v2[1];
            result[1] = v1[2] * v2[0] - v1[0] * v2[2];
            result[2] = v1[0] * v2[1] - v1[1] * v2[0];
            return result;
        }
        public static RVector ConvertToRVector(RMatrix A)
        {
            int nrows = A.GetnRows;
            int ncols = A.GetnCols;

            if (nrows != 1 && ncols != 1) { throw new Exception("A must contain either 1 column or 1 row! "); }
            RVector x;
            if (nrows == 1)
            {
                x = new RVector(ncols);
                for (int i = 0; i < ncols; i++) { x[i] = A[0, i]; }
            }
            else
            {
                x = new RVector(nrows);
                for (int i = 0; i < nrows; i++) { x[i] = A[i, 0]; }
            }
            return x;
        }

    }
}
