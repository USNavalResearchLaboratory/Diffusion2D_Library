using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Diffusion2D_Library
{
    /// <summary>
    /// Defines a matrix class that can store real-valued data
    /// </summary>
    public class RMatrix : ICloneable
    {
        // Fields
        private readonly int nRows;
        private readonly int nCols;
        private double[,] matrix;

        // Constructors
        public RMatrix(int nRows, int nCols)
        {
            this.nRows = nRows;
            this.nCols = nCols;
            matrix = new double[nRows, nCols];
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    matrix[i, j] = 0.0;
                }
            }
        }
        public RMatrix(double[,] matrix)
        {
            nRows = matrix.GetLength(0);
            nCols = matrix.GetLength(1);
            this.matrix = matrix;
        }
        public RMatrix(RVector rv)
        {
            nRows = 1;
            nCols = rv.GetRVectorSize;
            matrix = new double[1, rv.GetRVectorSize];
            for (int i = 0; i < rv.GetRVectorSize; i++)
            {
                matrix[0, i] = rv[i];
            }
        }
        public RMatrix IdentityMatrix()
        {
            RMatrix m = new(nRows, nCols);
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    if (i == j)
                    {
                        m[i, j] = 1.0;
                    }
                }
            }
            return m;
        }
        public RMatrix(double[] rv)
        {
            nRows = rv.Length;
            nCols = 1;
            matrix = new double[nRows, nCols];
            for (int i = 0; i < nRows; i++)
            {
                matrix[i, 0] = rv[i];
            }
        }

        // Accessors
        public int GetnRows
        { get { return nRows; } }
        public int GetnCols
        { get { return nCols; } }
        // Indexers
        public double this[int m, int n]
        {
            get
            {
                if (m < 0 || m > nRows)
                {
                    throw new Exception("m-th row is out of range!");
                }
                if (n < 0 || n > nCols)
                {
                    throw new Exception("n-th col is out of range!");
                }
                return matrix[m, n];
            }
            set { matrix[m, n] = value; }
        }
        // Override Methods
        public override string ToString()
        {
            string strMatrix = "(";
            for (int i = 0; i < nRows; i++)
            {
                string str = "";
                for (int j = 0; j < nCols - 1; j++)
                {
                    str += matrix[i, j].ToString() + ", ";
                }
                str += matrix[i, nCols - 1].ToString();
                if (i != nRows - 1 && i == 0)
                    strMatrix += str + "\n";
                else if (i != nRows - 1 && i != 0)
                    strMatrix += " " + str + "\n";
                else
                    strMatrix += " " + str + ")";
            }
            return strMatrix;
        }
        public override bool Equals(object obj)
        {
            return (obj is RMatrix matrix1) && Equals(matrix1);
        }
        public bool Equals(RMatrix cm)
        {
            return matrix == cm.matrix;
        }
        public override int GetHashCode()
        {
            return matrix.GetHashCode();
        }
        public static bool operator ==(RMatrix cm1, RMatrix cm2)
        {
            return cm1.Equals(cm2);
        }
        public static bool operator !=(RMatrix cm1, RMatrix cm2)
        {
            return !cm1.Equals(cm2);
        }
        public static bool operator >=(RMatrix rm1, RMatrix rm2)
        {
            int nrows1 = rm1.GetnRows;
            int ncols1 = rm1.GetnCols;
            int nrows2 = rm2.GetnRows;
            int ncols2 = rm2.GetnCols;

            bool gt = true;
            if (nrows1! != nrows2 || ncols1 != ncols2) { throw new Exception("Matrices do not have the same dimensions!"); }

            for (int i = 0; i < nrows1; i++)
            {
                for (int j = 0; j < ncols1; j++)
                {
                    if (rm1[i, j] < rm2[i, j]) { gt = false; }
                }
            }
            return gt;
        }
        public static bool operator <=(RMatrix rm1, RMatrix rm2)
        {
            int nrows1 = rm1.GetnRows;
            int ncols1 = rm1.GetnCols;
            int nrows2 = rm2.GetnRows;
            int ncols2 = rm2.GetnCols;

            bool lt = true;
            if (nrows1! != nrows2 || ncols1 != ncols2) { throw new Exception("Matrices do not have the same dimensions!"); }
            for (int i = 0; i < nrows1; i++)
            {
                for (int j = 0; j < ncols1; j++)
                {
                    if (rm1[i, j] > rm2[i, j]) { lt = false; }
                }
            }

            return lt;
        }
        public static RMatrix operator +(RMatrix rm)
        {
            return rm;
        }
        public static RMatrix operator +(RMatrix rm1, RMatrix rm2)
        {
            if (!RMatrix.CompareDimension(rm1, rm2))
            {
                throw new Exception("The dimensions of 2 matrices must be the same!");
            }
            RMatrix result = new(rm1.GetnRows, rm1.GetnCols);
            for (int i = 0; i < rm1.GetnRows; i++)
            {
                for (int j = 0; j < rm1.GetnCols; j++)
                {
                    result[i, j] = rm1[i, j] + rm2[i, j];
                }
            }
            return result;
        }
        public static RMatrix operator +(RMatrix rm, double cn)
        {
            RMatrix result = new(rm.GetnRows, rm.GetnCols);
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i, j] = rm[i, j] + cn;
                }
            }
            return result;
        }
        public static RMatrix operator +(double cn, RMatrix rm)
        {
            RMatrix result = new(rm.GetnRows, rm.GetnCols);
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i, j] = rm[i, j] + cn;
                }
            }
            return result;
        }
        public static RMatrix operator +(RVector rv, RMatrix rm)
        {
            if (rv.GetRVectorSize != rm.GetnRows)
            {
                throw new Exception("The dimensions of vector and matrix must be the same!");
            }
            int nvals = rv.GetRVectorSize;
            RMatrix result = new(nvals, 1);
            for (int i = 0; i < nvals; i++)
            {
                result[i, 0] = rv[i] + rm[i, 0];
            }
            return result;
        }
        public static RMatrix operator -(RMatrix rm)
        {
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    rm[i, j] = -rm[i, j];
                }
            }
            return rm;
        }
        public static RMatrix operator -(RMatrix rm1, RMatrix rm2)
        {
            if (!RMatrix.CompareDimension(rm1, rm2))
            {
                throw new Exception("The dimensions of two matrices must be the same!");
            }
            RMatrix result = new(rm1.GetnRows, rm1.GetnCols);
            for (int i = 0; i < rm1.GetnRows; i++)
            {
                for (int j = 0; j < rm1.GetnCols; j++)
                {
                    result[i, j] = rm1[i, j] - rm2[i, j];
                }
            }
            return result;
        }
        public static RMatrix operator -(RMatrix rm, double cn)
        {
            RMatrix result = new(rm.GetnRows, rm.GetnCols);
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i, j] = rm[i, j] - cn;
                }
            }
            return result;
        }
        public static RMatrix operator -(double cn, RMatrix rm)
        {
            RMatrix result = new(rm.GetnRows, rm.GetnCols);
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i, j] = cn - rm[i, j];
                }
            }
            return result;
        }
        public static RMatrix operator *(RMatrix rm, double cn)
        {
            RMatrix result = new(rm.GetnRows, rm.GetnCols);
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i, j] = rm[i, j] * cn;
                }
            }
            return result;
        }
        public static RMatrix operator *(double cn, RMatrix rm)
        {
            RMatrix result = new(rm.GetnRows, rm.GetnCols);
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i, j] = rm[i, j] * cn;
                }
            }
            return result;
        }
        public static RMatrix operator *(RMatrix m1, RMatrix m2)
        {
            if (m1.GetnCols != m2.GetnRows)
            {
                throw new Exception("# columns of the matrix 1 must = # columns of the matrix 2");
            }

            double ctmp;
            int nrows1 = m1.GetnRows;
            int ncols2 = m2.GetnCols;
            RMatrix result = new(nrows1, ncols2);

            //for (int i = 0; i < cm1.GetnRows; i++)
            //{
            //    for (int j = 0; j < cm2.GetnCols; j++)
            //    {
            //        ctmp = result[i, j];
            //        for (int k = 0; k < cm2.GetnRows; k++)
            //        {
            //            ctmp += cm1[i, k] * cm2[k, j];
            //        }
            //        result[i, j] = ctmp;
            //    }
            //}
            for (int i = 0; i < nrows1; i++)
            {
                for (int j = 0; j < ncols2; j++)
                {
                    ctmp = 0.0;
                    RVector row = m1.GetRowVector(i);
                    RVector col = m2.GetColVector(j);

                    for (int k = 0; k < m1.GetnCols; k++)
                    {
                        ctmp += row[k] * col[k];
                    }
                    result[i, j] = ctmp;
                }
            }
            return result;
        }
        public static RMatrix operator *(RMatrix rm, RVector rv)
        {
            int num_rows = rm.GetnRows;
            int num_cols = rv.GetRVectorSize;
            RMatrix Result = new(num_rows, num_cols);

            for (int i = 0; i < num_rows; i++)
            {
                for (int j = 0; j < num_cols; j++)
                {
                    Result[i, j] = rm[i, 0] * rv[j];
                }
            }
            return Result;
        }
        public static RMatrix operator /(RMatrix rm, double cn)
        {
            RMatrix result = new(rm.GetnRows, rm.GetnCols);
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i, j] = rm[i, j] / cn;
                }
            }
            return result;
        }
        public static RMatrix operator /(double cn, RMatrix rm)
        {
            RMatrix result = new(rm.GetnRows, rm.GetnCols);
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i, j] = rm[i, j] / cn;
                }
            }
            return result;
        }

        // Methods
        // Checks for a square matrix where #rows = #cols
        public bool IsSquared()
        {
            if (nRows == nCols)
                return true;
            else
                return false;
        }
        // Compares the dimension of two real matrices
        public static bool CompareDimension(RMatrix cm1, RMatrix cm2)
        {
            if (cm1.GetnRows == cm2.GetnRows && cm1.GetnCols == cm2.GetnCols)
                return true;
            else
                return false;
        }
        // Makes a clone copy of a real matrix
        public RMatrix Clone()
        {
            RMatrix cm = new(matrix);
            cm.matrix = (double[,])matrix.Clone();
            return cm;
        }
        object ICloneable.Clone()
        {
            return Clone();
        }

        // Sets up a call to calculate the transpose of a real matrix
        public RMatrix GetTranspose()
        {
            RMatrix ct = this;
            if (IsSquared())
            {
                ct.Transpose();
                return ct;
            }
            else
            {
                RMatrix nt = Transpose(ct);
                return nt;
            }
        }

        // Calculates the transpose of a real matrix
        public void Transpose()
        {
            RMatrix cm = new(nCols, nRows);
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    cm[j, i] = matrix[i, j];
                }
            }
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    matrix[i, j] = cm[i, j];
                }
            }

        }
        public static RMatrix Transpose(RMatrix cm)
        {
            RMatrix nm = new(cm.nCols, cm.nRows);
            for (int i = 0; i < cm.nRows; i++)
            {
                for (int j = 0; j < cm.nCols; j++)
                {
                    nm[j, i] = cm.matrix[i, j];
                }
            }
            return nm;
        }
        // Calculates the trace of a real matrix
        public double GetTrace()
        {
            double sum_of_diag = 0.0;
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    if (i == j)
                        sum_of_diag += matrix[i, j];
                }
            }
            return sum_of_diag;
        }
        // Extracts a row vector from a real matrix at specified row
        public RVector GetRowVector(int m)
        {
            if (m < 0 || m > nRows)
            {
                throw new Exception("m-th row is out of range!");
            }
            RVector RowVector = new(nCols);
            for (int i = 0; i < nCols; i++)
            {
                RowVector[i] = matrix[m, i];
            }
            return RowVector;
        }
        // Extracts a column vector from a real matrix at specified column
        public RVector GetColVector(int m)
        {
            if (m < 0 || m > nCols)
            {
                throw new Exception("n-th col is out of range!");
            }
            RVector ColCVector = new(nRows);
            for (int i = 0; i < nRows; i++)
            {
                ColCVector[i] = matrix[i, m];
            }
            return ColCVector;
        }
        // Swaps specificed real matrix row with another row
        public RMatrix SwapMatrixRow(int m, int n)
        {
            double ctemp = 0.0;
            for (int i = 0; i < nCols; i++)
            {
                ctemp = matrix[m, i];
                matrix[m, i] = matrix[n, i];
                matrix[n, i] = ctemp;
            }
            return new RMatrix(matrix);
        }
        // Swaps specificed real matrix column with another column
        public RMatrix SwapMatrixColumn(int m, int n)
        {
            double ctemp = 0.0;
            for (int i = 0; i < nRows; i++)
            {
                ctemp = matrix[i, m];
                matrix[i, m] = matrix[i, n];
                matrix[i, n] = ctemp;
            }
            return new RMatrix(matrix);
        }
        // Calculates the transform of a real matrix
        public static RVector RTransform(RMatrix rm, RVector rv)
        {
            RVector result = new(rv.GetRVectorSize);
            if (!rm.IsSquared())
            {
                throw new Exception("The matrix must be squared!");
            }
            if (rm.GetnCols != rv.GetRVectorSize)
            {
                throw new Exception("Vector size must = # rows in matrix");
            }
            for (int i = 0; i < rm.GetnRows; i++)
            {
                result[i] = 0.0;
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i] += rm[i, j] * rv[j];
                }
            }
            return result;
        }
        public static RVector RTransform(RVector rv, RMatrix rm)
        {
            RVector result = new(rv.GetRVectorSize);
            if (!rm.IsSquared())
            {
                throw new Exception("The matrix must be squared!");
            }
            if (rm.GetnRows != rv.GetRVectorSize)
            {
                throw new Exception("Vector size must = # rows in matrix");
            }
            for (int i = 0; i < rm.GetnRows; i++)
            {
                result[i] = 0.0;
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    result[i] += rv[j] * rm[j, i];
                }
            }
            return result;
        }
        public static RMatrix RTransform(RVector rv1, RVector rv2)
        {
            if (rv1.GetRVectorSize != rv2.GetRVectorSize)
            {
                throw new Exception("The vectors must have the same size!");
            }
            RMatrix result = new(rv1.GetRVectorSize, rv1.GetRVectorSize);
            for (int i = 0; i < rv1.GetRVectorSize; i++)
            {
                for (int j = 0; j < rv1.GetRVectorSize; j++)
                {
                    result[j, i] = rv1[i] * rv2[j];
                }
            }
            return result;
        }
        // Calculates the determinant of a real matrix
        public static double Determinant(RMatrix rm)
        {
            double result = 0.0;
            if (!rm.IsSquared())
            {
                throw new Exception("The matrix must be squared!");
            }
            if (rm.GetnRows == 1)
                result = rm[0, 0];
            else
            {
                for (int i = 0; i < rm.GetnRows; i++)
                {
                    double mDeterm = Determinant(Minor(rm, 0, i));
                    result += Math.Pow(-1, i) * rm[0, i] * mDeterm;
                }
            }
            return result;
        }
        // Calculates the minor of a real matrix at specified row and column
        public static RMatrix Minor(RMatrix rm, int row, int col)
        {
            RMatrix cmm = new(rm.GetnRows - 1, rm.GetnCols - 1);
            int ii = 0, jj = 0;
            for (int i = 0; i < rm.GetnRows; i++)
            {
                if (i == row)
                    continue;
                jj = 0;
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    if (j == col)
                        continue;
                    cmm[ii, jj] = rm[i, j];
                    jj++;
                }
                ii++;
            }
            return cmm;
        }
        // Calculates the adjoint of a real matrix
        public static RMatrix Adjoint(RMatrix rm)
        {
            if (!rm.IsSquared())
            {
                throw new Exception("The matrix must be squared!");
            }
            RMatrix ma = new(rm.GetnRows, rm.GetnCols);
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    ma[i, j] = Math.Pow(-1, i + j) * (Determinant(Minor(rm, i, j)));
                }
            }
            return ma.GetTranspose();
        }

        // Calculates the inverse of a real matrix
        public static RMatrix Inverse(RMatrix rm)
        {
            double Dm = Determinant(rm);
            if (Dm == 0)
            {
                throw new Exception("Cannot inverse a matrix with 0 determinant!");
            }
            return (Adjoint(rm) / Dm);
        }
        public static RMatrix InverseAccurate(RMatrix rm)
        {
            double Dm = Determinant(rm);
            if (Dm == 0)
            {
                throw new Exception("Cannot inverse a matrix with 0 determinant!");
            }
            int nrows = rm.GetnRows;
            int ncols = rm.GetnCols;
            RMatrix InvMat = new(nrows, ncols);
            bool rowsum = false;

            RMatrix D0 = new(nrows, ncols), D0B, I = new(nrows, ncols), C;
            //double start_pow = 0.1;
            //double delta_pow = 0.01;
            int counter = 0;
            while (!rowsum)
            {
                double themax = Max(rm);
                double apow = Math.Log10(Math.Abs(themax)) + 5;
                double the_pow = apow; // start_pow - counter * delta_pow;
                double multiplier = Math.Pow(10, -the_pow);
                D0 = multiplier * Inverse(rm);
                D0B = D0 * rm;
                I = rm.IdentityMatrix();
                C = I - D0B;

                RVector row_sum = new(nrows);
                for (int i = 0; i < nrows; i++)
                {
                    for (int j = 0; j < ncols; j++)
                    {
                        row_sum[i] += C[i, j];
                    }
                }

                int rowcounter = 0;
                for (int i = 0; i < nrows; i++)
                {
                    if (row_sum[i] < 1.0) { rowcounter++; }
                }
                if (rowcounter == nrows - 1) { rowsum = true; }
                counter++;
            }

            if (rowsum == true)
            {
                RMatrix Dnew = new(nrows, ncols);
                for (int c = 1; c < 30; c++)
                {
                    Dnew = (2 * I - D0 * rm) * D0;
                    if (Dnew == D0) { break; }
                    else { D0 = Dnew; }
                }
                InvMat = Dnew;
            }

            return InvMat;
        }
        // Replaces the n-th row of a real matrix with contents of a real vector
        public RMatrix ReplaceRow(RVector rv, int m)
        {
            if (m < 0 || m > nRows)
            {
                throw new Exception("m-th row is out of range!");
            }
            if (rv.GetRVectorSize != nCols)
            {
                throw new Exception("Vector size is out of range!");
            }
            for (int i = 0; i < nCols; i++)
            {
                matrix[m, i] = rv[i];
            }
            return new RMatrix(matrix);
        }
        // Replaces the n-th column of a real matrix with contents of a real vector
        public RMatrix ReplaceCol(RVector rv, int n)
        {
            if (n < 0 || n > nCols)
            {
                throw new Exception("n-th col is out of range!");
            }
            if (rv.GetRVectorSize != nRows)
            {
                throw new Exception("Vector size is out of range!");
            }
            for (int i = 0; i < nRows; i++)
            {
                matrix[i, n] = rv[i];
            }
            return new RMatrix(matrix);
        }
        public static (int rval, int cval) IndexOf(RMatrix rm, double c)
        {
            int ifound = -1, jfound = -1;
            for (int i = 0; i < rm.GetnRows; i++)
            {
                for (int j = 0; j < rm.GetnCols; j++)
                {
                    if (rm[i, j] == c)
                    {
                        ifound = i;
                        jfound = j;
                    }
                }
            }

            return (ifound, jfound);
        }
        public static RMatrix GetDiagonal(RMatrix rm)
        {
            int m = rm.GetnRows;
            int n = rm.GetnCols;
            RMatrix Result = new(m, n);

            for (int i = 0; i <= m - 1; i++)
            {
                for (int j = 0; j <= n - 1; j++)
                {
                    if (i == j)
                    {
                        Result[i, j] = rm[i, j];
                    }
                }
            }

            return Result;
        }
        public static void Cholesky_Decomposition(RMatrix rm)
        {
            int n = rm.GetnRows;
            RMatrix lower = new(n, n);

            // Decomposing a matrix
            // into Lower Triangular
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double sum = 0;

                    // summation for diagonals
                    if (j == i)
                    {
                        for (int k = 0; k < j; k++) { sum += (int)Math.Pow(lower[j, k], 2); }
                        lower[j, j] = Math.Sqrt(rm[j, j] - sum);
                    }

                    else
                    {
                        // Evaluating L(i, j)
                        // using L(j, j)
                        for (int k = 0; k < j; k++) { sum += (lower[i, k] * lower[j, k]); }
                        lower[i, j] = (rm[i, j] - sum) / lower[j, j];
                    }
                }
            }
        }
        public static double Max(RMatrix rm)
        {
            int nrows = rm.GetnRows;
            int ncols = rm.GetnCols;

            double max = 0.0;
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    if (Math.Abs(rm[i, j]) > max) { max = rm[i, j]; }
                }
            }
            return max;
        }
        public static double MinDiagonal(RMatrix rm)
        {
            int nrows = rm.GetnRows;
            int ncols = rm.GetnCols;

            double min = Math.Abs(rm[0, 0]);
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    if (i == j && Math.Abs(rm[i, j]) < min) { min = rm[i, j]; }
                }
            }
            return min;
        }
        public static RVector DotProduct(RMatrix rm, RVector rv)
        {
            int nrows = rm.GetnRows;
            int ncols1 = rm.GetnCols;
            int nrows2 = rv.GetRVectorSize;
            RVector result = new(nrows2);

            for (int i = 0; i < nrows; i++)
            {
                double col_sum = 0;
                for (int j = 0; j < ncols1; j++)
                {
                    double p = rm[i, j] * rv[j];
                    col_sum += p;
                }
                result[i] = col_sum;
            }
            return result;
        }
        public static RMatrix GaussSeidel(RMatrix A, RVector x0, RMatrix b)
        {
            double tol = 1.0e-8;
            bool stop = false;
            int iter = 0;
            double scale = Math.Pow(10, -6);

            int N = 300;
            int nrows = A.GetnRows; //2; // 
            int ncols = A.GetnCols; //2; // 
            RMatrix p = new(nrows, 1);
            RVector x = new(nrows);
            RVector y = new(nrows); ;
            double sum = 0.0;

            //A = new RMatrix(2, 2);
            //A[0, 0] = 16;
            //A[0, 1] = 3;
            ////A[0, 2] = 2;
            //A[1, 0] = 7;
            //A[1, 1] = -11;
            ////A[1, 2] = 1;
            ////A[2, 0] = 1;
            ////A[2, 1] = 1;
            ////A[2, 2] = 3;

            //b = new RMatrix(3, 1);
            //b[0, 0] = 11;
            //b[1, 0] = 13;
            ////b[2, 0] = 3;

            //x[0] = 0;
            //x[1] = 0;
            ////x[2] = 0;

            for (int m = 0; m < nrows; m++) { y[m] = 0.0; }

            while (!stop && iter < N)
            {
                for (int i = 0; i < nrows; i++)
                {
                    sum = b[i, 0];
                    for (int j = 0; j < ncols; j++)
                    {
                        if (j != i) { sum -= (A[i, j] * x[j]); }
                    }
                    x[i] = sum / A[i, i];
                }

                int check_sum = 0;
                for (int k = 0; k < nrows; k++)
                {
                    double check_term = Math.Abs((x[k] - y[k]) / x[k]);
                    if (check_term < tol) { check_sum += 1; }

                    if (check_sum >= nrows) { stop = true; }
                    else { y[k] = x[k]; }

                }
                iter += 1;
            }


            //RMatrix D = GetDiagonal(A);
            //RMatrix E = A - D;
            //RMatrix inVD = Inverse(D);
            //RMatrix intx0 = new RMatrix(nrows, 1);
            //for (int i = 0; i < nrows; i++)
            //{
            //    intx0[i, 0] = scale * x0[i];
            //}

            //p = (inVD * b) - (inVD * (E * intx0));
            for (int i = 0; i < nrows; i++)
            {
                p[i, 0] = x[i];
            }
            return p;
        }
        //public static RMatrix GetScaledAMatrix(RMatrix A)
        //{
        //    int nrows = A.GetnRows;
        //    int ncols = A.GetnCols;
        //    RMatrix A_prime = new RMatrix(nrows, ncols);

        //    double denom1 = 0.0;
        //    double denom2 = 0.0;
        //    double denom3 = 0.0;
        //    for (int i = 0; i < nrows; i++)
        //    {
        //        for (int j = 0; j < ncols; j++)
        //        {
        //            denom1 = Math.Sqrt(Math.Abs(A[i, i]));
        //            denom2 = Math.Sqrt(Math.Abs(A[j, j]));
        //            denom3 = 1.0 / (denom1 * denom2);
        //            if (!double.IsInfinity(denom3)) { A_prime[i, j] = denom3 * A[i, j]; }
        //            else { throw new Exception("Infinity detected!"); }

        //        }
        //    }

        //    return A_prime;
        //}
        //public static RMatrix GetScaledAMatrix(RMatrix A)
        //{
        //    int nrows = A.GetnRows;
        //    int ncols = A.GetnCols;
        //    RMatrix A_prime = new RMatrix(nrows, ncols);

        //    double tr_scale = Math.Sqrt(A.GetTrace()); //A.GetTrace(); //  
        //    double denom1 = 0.0, denom2 = 0.0;
        //    for (int i = 0; i < nrows; i++)
        //    {                
        //        for (int j = 0; j < ncols; j++)
        //        {
        //            denom1 = Math.Sqrt(Math.Abs(A[i, i]));
        //            denom2 = Math.Sqrt(Math.Abs(A[i, j]));
        //            A_prime[i, j] = A[i, j] / (denom1 * denom2);
        //        }
        //    }

        //    return A_prime;
        //}
        public static double L2Norm(RMatrix A)
        {
            int nrows = A.GetnRows;
            int ncols = A.GetnCols;

            double sum = 0.0;
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++) { sum += Math.Pow(A[i, j], 2); }
            }
            double value = Math.Sqrt(sum);
            return value;
        }
        public static double L2_Error_Norm(RMatrix A, RMatrix B, double h)
        {
            int nrows1 = A.GetnRows;
            int ncols1 = A.GetnCols;
            int nrows2 = B.GetnRows;
            int ncols2 = B.GetnCols;

            if (nrows1! != nrows2 || ncols1 != ncols2) { throw new Exception("Matrices do not have the same dimensions!"); }

            double sum = 0.0;
            for (int i = 0; i < nrows1; i++)
            {
                for (int j = 0; j < ncols1; j++) { sum += Math.Pow(h, 2) * Math.Pow(A[i, j] - B[i, j], 2); }
            }

            double value = Math.Sqrt(sum);
            return value;
        }

        public static RVector Cholesky(RMatrix A, RVector b)
        {
            int nrows = A.GetnRows;
            int ncols = A.GetnCols;
            RVector x = new(b.GetRVectorSize);
            RVector y = new(b.GetRVectorSize);

            RMatrix L = new(nrows, ncols);
            RMatrix U = new(nrows, ncols);

            L[0, 0] = Math.Sqrt(A[0, 0]);
            if (nrows == 1 && ncols == 1)
            {
                y[0] = b[0] / L[0, 0];
                x[0] = y[0] / L[0, 0];
            }
            else
            {
                for (int i = 0; i < nrows; i++) { L[i, 0] = A[i, 0] / L[0, 0]; }

                for (int i = 1; i < nrows-1; i++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < i; k++) { sum += L[i, k] * L[i, k]; }
                    L[i, i] = Math.Sqrt(A[i, i] - sum);

                    for (int j = i + 1; j < ncols; j++)
                    {
                        double sum2 = 0.0;
                        for (int k = 0; k < i; k++) { sum2 += L[j, k] * L[i, k]; }
                        L[j, i] = (A[j, i] - sum2) / L[i, i];
                    }
                }

                double sum3 = 0.0;
                for (int k = 0; k < ncols; k++) { sum3 += L[nrows - 1, k] * L[nrows - 1, k]; }
                L[nrows - 1, ncols - 1] = Math.Sqrt(A[nrows - 1, ncols - 1] - sum3);

                y[0] = b[0] / L[0, 0];

                for (int i = 1; i < nrows; i++)
                {
                    double sum4 = 0.0;
                    for (int k = 0; k < i; k++) { sum4 += L[i, k] * y[k]; }
                    y[i] = (b[i] - sum4) / L[i, i];
                }

                x[nrows - 1] = y[nrows - 1] / L[nrows - 1, ncols - 1];
                for (int i = nrows-2; i > -1; i--)
                {
                    double sum5 = 0.0;
                    for (int k = i + 1; k < nrows; k++) { sum5 += L[k, i] * x[k]; }
                    x[i] = (y[i] - sum5) / L[i, i];
                }
            }

            return x;
        }
    }
}
