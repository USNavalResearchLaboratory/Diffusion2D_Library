using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Diffusion2D_Library
{
    class BoundaryAndInitialConditionsMethods
    {
        public static RVector Boundary_2D_Constant(double t, RVector SlidingSide, double FixedValue)
        {
            int n = SlidingSide.GetRVectorSize;
            RVector BC = new(n);
            //double y = FixedValue;

            for (int i = 0; i < n; i++)
            {
                //double x = SlidingSide[i];
                BC[i] = 55.55; // Math.Exp(x * y * t);
            }

            return BC;
        }
        public static RVector Boundary_2D_Zero(double t, RVector SlidingSide, double FixedValue)
        {
            int n = SlidingSide.GetRVectorSize;
            RVector BC = new(n);
            //double x = FixedValue;

            for (int i = 0; i < n; i++)
            {
                //double y = SlidingSide[i];
                BC[i] = 0.0; // Math.Exp(x * y * t);
            }

            return BC;
        }
        public static RMatrix InitialCondition_2D(RMatrix xposition, RMatrix yposition)
        {
            int nrows = xposition.GetnRows;
            int ncols = xposition.GetnCols;


            RMatrix composition_field = new(nrows, ncols);

            //for (int i = 0; i < nrows; i++)
            //{
            //    for (int j = 0; j < ncols; j++)
            //    {
            //        //if (i >= nx_start && i <= nx_end && j >= ny_start && j <= ny_end) { composition_field[i, j] = 0.0; }
            //        //else { composition_field[i, j] = 0.0; }
            //        double x = xposition[i, j];
            //        double y = yposition[i, j];
            //        composition_field[i, j] = 0.0; // Math.Exp(x * y * 0.0);
            //    }
            //}

            return composition_field;
        }
        public static RMatrix SourceTerm_2D(RMatrix xposition, RMatrix yposition, double time, RMatrix D, RMatrix composition)
        {
            //double t = time;
            int nrows = xposition.GetnRows;
            int ncols = xposition.GetnCols;

            RMatrix C_field = new(nrows, ncols);

            //for (int i = 1; i < nrows - 1; i++)
            //{
            //    for (int j = 1; j < ncols - 1; j++)
            //    {
            //        double x = xposition[i, j];
            //        double y = yposition[i, j];


            //        double term1 = Math.Exp(x * y * t);
            //        double term2 = x * y;
            //        double term3 = Math.Pow(y * t, 2);
            //        double term4 = Math.Pow(x * t, 2);

            //        C_field[i, j] = 0.0; // (term2 - (D[i, j] * (term3 + term4))) * term1;
            //    }
            //}
            return C_field;
        }
    }
}
