using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static Diffusion2D_Library.BoundaryCondition;

namespace Diffusion2D_Library
{
    /// <summary>
    /// Delegate for methods for boundary condition functions
    /// </summary>
    /// <param name="time"></param>
    /// <returns>A vector storing the composition on a boundary of a 2D field</returns>
    //public delegate RVector BoundaryCondition_Del(double time, int n);
    public delegate RVector Del_BC_Method(double t, RVector SlidingSide, double FixedValue);
    /// <summary>
    /// Delegate for assigning the initial composition field in a 2D region
    /// </summary>
    /// <param name="c"></param>
    /// <returns>A matrix of composition values</returns>
    public delegate RMatrix Del_IC_xy(RMatrix x, RMatrix y);
    /// <summary>
    /// Delegate for calculating sources/sinks for reactions occurring in a 2D composition field
    /// </summary>
    /// <param name="position"></param>
    /// <param name="time"></param>
    /// <param name="composition"></param>
    /// <returns>A matrix of composition values</returns>
    public delegate RMatrix Del_Source_xy(RMatrix xposition, RMatrix yposition, double time, RMatrix DiffCoeff, RMatrix composition);
    /// <summary>
    /// Delegate for calculating sources/sinks for reactions occurring in a 2D composition field with diffusivity as a matrix
    /// </summary>
    /// <param name="xposition"></param>
    /// <param name="yposition"></param>
    /// <param name="time"></param>
    /// <param name="DiffCoeff"></param>
    /// <param name="composition"></param>
    /// <returns></returns>
    public delegate RMatrix SourceTerm_MatrixD_Del(RMatrix xposition, RMatrix yposition, double time, RMatrix DiffCoeff, RMatrix composition);
    public delegate RVector Del_BC_Constant(int n, double v);
    public delegate RMatrix Del_IC_Constant(int nr, int nc, double v);
    public delegate RMatrix Del_Source_Constant(int nr, int nc, double v);
    public class DiffusionSimulator2D
    {
        public struct CompositionField2D
        {
            public RMatrix InitialComposition;
            public RMatrix FinalComposition;
            public RMatrix xposition_matrix;
            public RMatrix yposition_matrix;
            public RMatrix DiffusionCoefficient;
            public CompositionField2D(int nx, int ny)
            {
                InitialComposition = new(nx, ny);
                FinalComposition = new(nx, ny);
                xposition_matrix = new(nx, ny);
                yposition_matrix = new(nx, ny);
                DiffusionCoefficient = new(nx, ny);
            }
        }
        

        /// <summary>
        /// Enumeration that specifies the type of boundary condition to be applied to a boundary
        /// </summary>
        /// 
        public enum BoundingBox
        {
            left,
            right,
            top,
            bottom
        }
        public enum Mode { Quiet, Verbose }
        public struct BoundaryWithFunction
        {
            public BoundingBox BoundaryLocation;
            public ABoundaryCondition TypeBC;
            public Del_BC_Method BoundaryFunction;
            public RVector PositionVaries;
            public double PositionFixed;
        }
        public struct BoundaryWithVector
        {
            public BoundingBox BoundaryLocation;
            public ABoundaryCondition TypeBC;
            public RVector FunctionValues;
        }
    }
}
