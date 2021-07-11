using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Diffusion2D_Library
{
    class BoundaryAndInitialConditionsMethods
    {
        /// <summary>
        /// Method for a boundary condition with a specified constant value
        /// </summary>
        /// <param name="n">Size of the vector of entries</param>
        /// <param name="constant_value">Constant value to be specified on the boundary</param>
        /// <returns>A vector of the specified size with the specified constant value in each entry</returns>
        public static RVector Boundary_2D_Constant(int n, double constant_value)
        {
            RVector BC = new(n);
            for (int i = 0; i < n; i++) { BC[i] = constant_value; }
            return BC;
        }
        /// <summary>
        /// Method for a boundary condition with a zero value
        /// </summary>
        /// <param name="n">Size of the vector to be returned</param>
        /// <returns>A vector of the specified size with a zero in each entry</returns>
        public static RVector Boundary_2D_Zero(int n)
        {   
            RVector BC = new(n);         
            return BC;
        }
        /// <summary>
        /// Method to create an initial composition matrix with 0 in all entried
        /// </summary>
        /// <param name="nrows">Number of rows in the matrix</param>
        /// <param name="ncols">Number of columns in the matrix</param>
        /// /// <param name="value">Constant value to set the initial condition</param>
        /// <returns>A matrix with a constant value in all of its entries</returns>
        public static RMatrix InitialCondition_2D_Zero(int nrows, int ncols, double value)
        {
            RMatrix composition_field = new(nrows, ncols);
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    composition_field[i, j] = value;
                }
            }
            return composition_field;
        }
        /// <summary>
        /// A method to create a source term that returns a specified vlaue in all of its entries
        /// </summary>
        /// <param name="nrows">Number of rows in the matrix</param>
        /// <param name="ncols">Numbe of columns in the matrix</param>
        /// <returns>A matrix with the specified value in all of its rows and columns</returns>
        public static RMatrix SourceTerm_2D(int nrows, int ncols, double value)
        {
            RMatrix composition_field = new(nrows, ncols);
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    composition_field[i, j] = value;
                }
            }
            return composition_field;
        }
    }
}
