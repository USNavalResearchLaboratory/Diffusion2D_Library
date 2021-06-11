using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Diffusion2D_Library
{
    /// <summary>
    /// Defines a tri-diagonal, banded, matrix class for solving martix equations
    /// </summary>
    public class TridiagonalMatrix
    {
        // Fields
        private readonly int nRows;
        private readonly int nCols;
        private readonly RVector main;
        private readonly RVector lower;
        private readonly RVector upper;

        // Constructors
        public TridiagonalMatrix(int nRows, int nCols)
        {
            if (!IsSquared()) { throw new Exception("This matrix needs to be square!"); }

            this.nRows = nRows;
            this.nCols = nCols;
            main = new RVector(nRows);
            lower = new RVector(nRows - 1);
            upper = new RVector(nRows - 1);

            for (int i = 0; i < nRows; i++)
            {
                main[i] = 1.0;
            }
            for (int i = 0; i < nRows - 1; i++)
            {
                lower[i] = 1.0;
                upper[i] = 1.0;
            }
        }
        public TridiagonalMatrix(RVector main, RVector upper, RVector lower)
        {
            this.main = main;
            if (upper.GetRVectorSize == main.GetRVectorSize - 1)
            {
                this.upper = upper;
            }
            if (lower.GetRVectorSize == main.GetRVectorSize - 1)
            {
                this.lower = lower;
            }

        }
        public TridiagonalMatrix(int N, double mainDiagonal, double upperVec, double lowerVec)
        {
            nRows = N;
            nCols = N;

            main = new RVector(N);
            for (int i = 0; i < main.GetRVectorSize; i++) { main[i] = mainDiagonal; }
            upper = new RVector(N - 1);
            for (int i = 0; i < upper.GetRVectorSize; i++) { upper[i] = upperVec; }
            lower = new RVector(N - 1);
            for (int i = 0; i < lower.GetRVectorSize; i++) { lower[i] = lowerVec; }
        }

        // Accessors
        public int GetnRows
        { get { return nRows; } }
        public int GetnCols
        { get { return nCols; } }

        // Indexers
        public double this[int col_idx, int row_idx]
        {
            get
            {
                if (col_idx < 0 || col_idx > nRows)
                {
                    throw new Exception("m-th row is out of range!");
                }
                if (row_idx < 0 || row_idx > nCols)
                {
                    throw new Exception("n-th col is out of range!");
                }
                double out_value;
                if (row_idx == col_idx) { out_value = main[col_idx]; }
                else if (row_idx == col_idx - 1) { out_value = lower[col_idx - 1]; }
                else if (row_idx == col_idx + 1) { out_value = upper[col_idx]; }
                else
                {
                    out_value = 0.0;
                }
                return out_value;
            }
            set
            {
                if (row_idx == col_idx) { main[col_idx] = value; }
                else if (row_idx == col_idx - 1) { lower[col_idx - 1] = value; }
                else if (row_idx == col_idx + 1) { upper[col_idx] = value; }
                //else
                //{
                //    double avalue = 0.0;
                //    //Console.WriteLine("Trying to set a value outside of the tridiagonal band.");
                //}
            }
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
        public RVector GetRow(int row_idx)
        {
            if (row_idx < 0 || row_idx > nRows - 1)
            {
                throw new Exception("r-th row is out of range!");
            }
            RVector RowRVector = new(nCols);
            if (row_idx == 0)
            {
                RowRVector[row_idx] = main[row_idx];
                RowRVector[row_idx + 1] = upper[row_idx];
            }
            else if (row_idx == main.GetRVectorSize - 1)
            {
                RowRVector[row_idx - 1] = lower[row_idx - 1];
                RowRVector[row_idx] = main[row_idx];
            }
            else
            {
                RowRVector[row_idx - 1] = lower[row_idx - 1];
                RowRVector[row_idx] = main[row_idx];
                RowRVector[row_idx + 1] = upper[row_idx];
            }

            return RowRVector;
        }
        public RVector GetCol(int c)
        {
            if (c < 0 || c > nCols - 1)
            {
                throw new Exception("c-th column is out of range!");
            }
            RVector ColRVector = new(nRows);

            if (c == 0)
            {
                ColRVector[c] = main[c];
                ColRVector[c + 1] = lower[c];
            }
            else if (c == main.GetRVectorSize - 1)
            {
                ColRVector[c - 1] = upper[c - 1];
                ColRVector[c] = main[c];
            }
            else
            {
                ColRVector[c - 1] = upper[c - 1];
                ColRVector[c] = main[c];
                ColRVector[c + 1] = lower[c - 1];
            }
            return ColRVector;
        }
        public RVector Dot(RVector rv)
        {
            if (GetnCols != rv.GetRVectorSize)
            {
                throw new Exception("# columns of the matrix must = # rows of the vector");
            }

            double ctmp;
            int nrows1 = GetnRows;
            int ncols2 = rv.GetRVectorSize;
            int n = rv.GetRVectorSize;
            RVector result = new(n);

            for (int i = 0; i < n; i++)
            {
                RVector row = GetRow(i);
                if (i == 0) { ctmp = (row[i] * rv[i]) + (row[i + 1] * rv[i + 1]); }
                else if (i == n - 1)
                {
                    ctmp = (row[i - 1] * rv[i - 1]) + (row[i] * rv[i]);
                }
                else
                {
                    double t_term = row[i - 1] * rv[i - 1];
                    double c_term = row[i] * rv[i];
                    double l_term = row[i + 1] * rv[i + 1];
                    ctmp = t_term + c_term + l_term;
                }

                result[i] = ctmp;
            }
            return result;
        }

        // Public methods for solving linear algebraic expressions
        public static RVector GaussSeidel_Old(TridiagonalMatrix A, RVector xold, RVector b)
        {
            int n = xold.GetRVectorSize;
            RVector xnew = new(n);

            for (int i = 0; i < n; i++)
            {
                double Lval = 0.0;
                for (int j = 0; j < i - 1; j++)
                {
                    Lval += A[i, j] * xnew[j];
                }
                double Rval = 0.0;
                for (int k = i + 1; k < n; k++)
                {
                    Rval += A[i, k] * xold[k];
                }
                xnew[i] = (1.0 / A[i, i]) * (b[i] - Lval - Rval);
            }

            return xnew;
        }
        public static RVector Jacobi(TridiagonalMatrix A, RVector XO, RVector b)
        {
            int n = XO.GetRVectorSize;
            RVector x = new(n);
            double tol = 1.0e-5;
            int N = 100;
            double sum;
            int outK;
            for (int k = 0; k < N; k++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (j == 0) { sum = (-A[j, j + 1] * XO[j + 1]) + b[j]; }
                    else if (j == n - 1) { sum = (-A[j, j - 1] * XO[j - 1]) + b[j]; }
                    else { sum = (-A[j, j - 1] * XO[j - 1]) + (-A[j, j + 1] * XO[j + 1]) + b[j]; }

                    x[j] = sum / A[j, j];
                }
                RVector normx = x - XO;
                double normx_val = normx.GetNorm();

                if (normx_val < tol) { outK = k; break; }

                for (int i = 0; i < n; i++) { XO[i] = x[i]; }

            }
            return x;
        }
        public static RVector GaussSeidel(TridiagonalMatrix A, RVector b) //RVector XO, 
        {
            int n = b.GetRVectorSize;
            RVector x = new(n);
            RVector XO = new(n);
            for (int i = 0; i < n; i++) { XO[i] = b[i]; }
            double tol = 1.0e-5;
            int N = 100;
            double sum;
            int outK;
            double Ajjm1, Ajjp1, Ajj;

            for (int k = 0; k < N; k++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (j == 0) { Ajjp1 = A[j, j + 1]; sum = b[j] - (Ajjp1 * XO[j + 1]); }
                    else if (j == n - 1) { Ajjm1 = A[j, j - 1]; sum = b[j] - (Ajjm1 * x[j - 1]); }
                    else
                    {
                        Ajjm1 = A[j, j - 1];
                        Ajjp1 = A[j, j + 1];
                        sum = b[j] - (Ajjm1 * x[j - 1]) - (Ajjp1 * XO[j + 1]);
                    }
                    Ajj = A[j, j];
                    x[j] = (1.0 / Ajj) * (sum);
                }
                RVector normx = x - XO;
                double normx_val = normx.GetNorm();

                if (normx_val < tol) { outK = k; break; }

                for (int i = 0; i < n; i++) { XO[i] = x[i]; }

            }
            //x[0] = b[0];
            //x[n - 1] = b[n - 1];
            return x;
        }
        public static RVector SOR(TridiagonalMatrix A, RVector b)
        {
            double omega = 1.24;
            int n = b.GetRVectorSize;
            RVector x = new(n);
            RVector XO = new(n);
            for (int i = 0; i < n; i++) { XO[i] = 0.0; }
            double tol = 1.0e-5;
            int N = 100;
            double sum;
            int outK;
            double Ajjm1, Ajjp1, Ajj;
            for (int k = 0; k < N; k++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (j == 0) { Ajjp1 = A[j, j + 1]; sum = b[j] - (Ajjp1 * XO[j + 1]); }
                    else if (j == n - 1) { Ajjm1 = A[j, j - 1]; sum = b[j] - (Ajjm1 * x[j - 1]); }
                    else
                    {
                        Ajjm1 = A[j, j - 1];
                        Ajjp1 = A[j, j + 1];
                        sum = b[j] - (Ajjm1 * x[j - 1]) - (Ajjp1 * XO[j + 1]);
                    }
                    Ajj = A[j, j];
                    x[j] = ((1 - omega) * XO[j]) + (omega / Ajj) * (sum);
                }
                RVector normx = x - XO;
                double normx_val = normx.GetNorm();

                if (normx_val < tol) { outK = k; break; }

                for (int i = 0; i < n; i++) { XO[i] = x[i]; }
            }
            return x;
        }
        public static RVector Thomas_Algorithm(TridiagonalMatrix A, RVector b)
        {
            int n = b.GetRVectorSize;
            RVector x = new(n);
            RVector c = new(n);
            RVector d = new(n);

            for (int i = 0; i < n - 1; i++)
            {
                if (i == 0) { c[i] = A[i, i + 1] / A[i, i]; }
                else { c[i] = A[i, i + 1] / (A[i, i] - (A[i, i - 1] * c[i - 1])); }
            }

            for (int i = 0; i < n; i++)
            {
                if (i == 0) { d[i] = b[i] / A[i, i]; }
                else { d[i] = (b[i] - A[i, i - 1] * d[i - 1]) / (A[i, i] - (A[i, i - 1] * c[i - 1])); }
            }

            for (int i = n - 1; i > -1; i--)
            {
                if (i == n - 1) { x[i] = d[i]; }
                else { x[i] = d[i] - (c[i] * x[i + 1]); }
            }
            return x;
        }

        // Override Methods
        public override string ToString()
        {
            int lcount = 0;
            int ucount = 0;

            string strMatrix = "{ ";
            Parallel.For(0, nRows, i =>
            {
                string str = "";
                for (int j = 0; j < nCols; j++)
                {
                    // Main diagonal
                    if (i == j)
                    {
                        if (i < nRows - 1 && j < nCols - 1)
                        {
                            str += main[i].ToString() + ", ";
                        }
                        else
                        {
                            str += main[i].ToString() + " }";
                        }

                    }
                    else
                    {
                        if (j == 0)
                        {
                            if (j == i - 1)
                            {
                                str += "{ " + lower[lcount].ToString() + ", "; lcount += 1;
                            }
                            else
                            {
                                str += "{ 0, ";
                            }

                        }
                        else if (j == nCols - 1)
                        {
                            if (j == i + 1)
                            {
                                str += upper[ucount].ToString() + " }" + "\n"; ucount += 1;
                            }
                            else
                            {
                                str += "0 }" + "\n";
                            }

                        }
                        else
                        {
                            if (j == i - 1)
                            {
                                str += lower[lcount].ToString() + ", "; lcount += 1;
                            }
                            else if (j == i + 1)
                            {
                                str += upper[ucount].ToString() + ", "; ucount += 1;
                            }
                            else
                            {
                                str += "0, ";
                            }

                        }
                    }

                    //if (i == j-1) { str += lower[lcount].ToString() + ", "; lcount += 1; }
                    //if (i == j+1 && j < nCols-1) { str += upper[ucount].ToString() + ", "; ucount += 1; }
                    //else if ((i == j + 1) && (j == nCols - 1)) { str += upper[ucount].ToString() + "\n "; ucount += 1; }

                    //if ((i != j-1) && (i !=j) && (i != j+1) && j < nCols - 1)  { str += "0, "; }                    
                    //if ((j == nCols - 1) && (j != i)  && (j != i+1) ) { str += "0}\n"; }
                }
                strMatrix += str;

                //if (i != nRows - 1 && i == 0)
                //    strMatrix += str + "\n";
                //else if (i != nRows - 1 && i != 0)
                //    strMatrix += " " + str + "\n";
                //else
                //    strMatrix += " " + str + ")";
            });
            return strMatrix;
        }
    }
}
