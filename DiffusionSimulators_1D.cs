using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Diffusion2D_Library
{
    /// <summary>
    /// Solves the parabolic partial differential equation: ∂u/∂t-(∂^2 u)/(∂x^2 )=f(x,t) in one dimension
    /// </summary>
    public class DiffusionSimulators_1D
    {
        /// <summary>
        /// Method delegate for boundary conditions on the composition field
        /// </summary>
        /// <param name="time"></param>
        /// <returns></returns>
        public delegate double BoundaryCondition_Del(double time);
        /// <summary>
        /// Method delegate for the initial condition of the composition field
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public delegate double InitialCondition_Del(double x);
        /// <summary>
        /// Method delegate to handle a function for a source term
        /// </summary>
        /// <param name="x">Position</param>
        /// <param name="t">Time</param>
        /// <returns></returns>
        public delegate RVector SourceTerm_Del(RVector x, double t, RVector c);
        private enum BoundaryConditions
        {
            Dirichlet,
            Neumann,
            Mixed,
            Unknown
        }

        // Fields
        private readonly double D;
        private readonly double dt;
        private readonly double dx;
        private RVector cinitial;
        private RVector cfinal;
        private readonly RVector position;
        private BoundaryConditions bc_type;

        // Properties
        public RVector C_Initial
        {
            get { return cinitial; }
            set
            {
                if (value.GetRVectorSize > 0) { cinitial = value; }
            }
        }
        public RVector C_Final
        {
            get { return cfinal; }
            set { if (value.GetRVectorSize > 0) { cfinal = value; } }
        }
        public string Boundary_Conditions
        {
            get { return ConvertEnumBCToString(bc_type); }
            set { bc_type = ConvertStringToEnumBC(value); }
        }

        // Constructors
        public DiffusionSimulators_1D(double D, double dx, int nx, double dt, int nt, string Boundary_Conditions)
        {
            this.D = D;
            this.dx = dx;
            this.dt = dt;
            cinitial = new(nx);
            cfinal = new(nx);
            position = new(nx);
            for (int i = 0; i < nx; i++) { position[i] = i * dx; }
            this.Boundary_Conditions = Boundary_Conditions;

        }
        public DiffusionSimulators_1D(double[] coeffs, int[] n, string Boundary_Conditions)
        {
            if (coeffs.Length >= 3)
            {
                D = coeffs[0];
                dx = coeffs[1];
                dt = coeffs[2];
            }
            cinitial = new(n[1]);
            cfinal = new(n[1]);
            position = new(n[0]);
            this.Boundary_Conditions = Boundary_Conditions;
        }

        // Public methods for solving the 1D diffusion equation
        /// <summary>
        /// Method for solving the 1D diffusion equation using the explicit Forward-Euler algorithm
        /// </summary>
        /// <param name="Lbc">Function for the left-side boundary condition for composition</param>
        /// <param name="Rbc">Function for the right-side boundary condition for composition</param>
        /// <param name="I0">Function for the initial condition for composition</param>
        /// <param name="n_steps">Number of time-steps</param>
        public void OneD_FE(BoundaryCondition_Del Lbc, BoundaryCondition_Del Rbc, InitialCondition_Del I0, SourceTerm_Del g, int n_steps)
        {
            int n = C_Initial.GetRVectorSize;
            RVector xold = new(n);
            RVector xnew = new(n);
            RVector b = new(n);

            // Define the A matrix
            double nu = D * dt / Math.Pow(dx, 2);
            double off_d_val = nu;
            double diag_val = 1 - (2 * nu);
            TridiagonalMatrix A = new(n, diag_val, off_d_val, off_d_val);
            switch (bc_type)
            {
                case BoundaryConditions.Dirichlet:
                    A[0, 0] = 1.0;
                    A[0, 1] = 0.0;
                    A[n - 1, n - 1] = 1.0;
                    A[n - 1, n - 2] = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        if (i == 0) { C_Initial[i] = Lbc(0.0) + I0(i * dx); }
                        else if (i == n - 1) { C_Initial[i] = Rbc(0.0) + I0(i * dx); }
                        else { C_Initial[i] = I0(i * dx); }

                        xold[i] = C_Initial[i];
                    }
                    break;
                case BoundaryConditions.Neumann:
                    A[0, 0] = -2 * nu;
                    A[0, 1] = 2 * nu;
                    A[n - 1, n - 1] = -2 * nu;
                    A[n - 1, n - 2] = 2 * nu;
                    for (int i = 0; i < n; i++) { C_Initial[i] = I0(i * dx); xold[i] = C_Initial[i]; }
                    break;
                case BoundaryConditions.Mixed:
                    break;
                case BoundaryConditions.Unknown:
                    break;
                default:
                    break;
            }

            // Time evolution
            for (int t = 0; t < n_steps; t++)
            {
                switch (bc_type)
                {
                    case BoundaryConditions.Dirichlet:
                        xnew = A.Dot(xold) + g(position, t * dt, xold);
                        xnew[0] = nu * Lbc(t * dt);
                        xnew[n - 1] = nu * Rbc(t * dt);
                        break;
                    case BoundaryConditions.Neumann:
                        xnew = A.Dot(xold) + g(position, t * dt, xold);
                        xnew[0] = -2 * nu * Lbc(t * dt);
                        xnew[n - 1] = 2 * nu * Rbc(t * dt);
                        break;
                    case BoundaryConditions.Mixed:
                        xnew = A.Dot(xold) + b;
                        break;
                }

                xold = xnew;
            }
            cfinal = xnew;
        }
        /// <summary>
        /// Method for solving the 1D diffusion equation using the implicit Backward-Euler algorithm
        /// </summary>
        /// <param name="Lbc">Function for the left-side boundary condition for composition</param>
        /// <param name="Rbc">Function for the right-side boundary condition for composition</param>
        /// <param name="I0">Function for the initial condition for composition</param>
        /// <param name="n_steps">Number of time-steps</param>
        public void OneD_BE(BoundaryCondition_Del Lbc, BoundaryCondition_Del Rbc, InitialCondition_Del I0, SourceTerm_Del g, int n_steps)
        {
            int n = C_Initial.GetRVectorSize;
            RVector b = new(n);
            RVector xnew = new(n);

            // Define the A matrix
            double nu = D * dt / Math.Pow(dx, 2);
            double off_d_val = -nu;
            double diag_val = 1 + (2 * nu);
            TridiagonalMatrix A = new(n, diag_val, off_d_val, off_d_val);
            switch (bc_type)
            {
                case BoundaryConditions.Dirichlet:
                    A[0, 0] = 1.0;
                    A[0, 1] = 0.0;
                    A[n - 1, n - 1] = 1.0;
                    A[n - 1, n - 2] = 0.0;
                    for (int i = 0; i < n; i++)
                    {
                        if (i == 0) { C_Initial[i] = Lbc(0.0) + I0(i * dx); }
                        else if (i == n - 1) { C_Initial[i] = Rbc(0.0) + I0(i * dx); }
                        else { C_Initial[i] = I0(i * dx); }
                    }

                    b = C_Initial;  // + g(position, 0.0, C_Initial)
                    break;
                case BoundaryConditions.Neumann:
                    A[0, 1] = -2 * nu;
                    A[n - 1, n - 2] = -2 * nu;
                    for (int i = 0; i < n; i++) { C_Initial[i] = I0(i * dx); b[i] = C_Initial[i]; }
                    break;
                case BoundaryConditions.Mixed:
                    break;
            }

            // Time evolution
            for (int t = 0; t < n_steps; t++)
            {
                //xreturn = TridiagonalMatrix.Jacobi_Solver(A, x_init_guess, b);
                //xnew = TridiagonalMatrix.GaussSeidel(A, b); //x_init_guess, 
                xnew = TridiagonalMatrix.SOR(A, b);  //+ g(position, t * dt, b)
                                                     //xnew = TridiagonalMatrix.Thomas_Algorithm(A, b);
                for (int i = 1; i < n - 1; i++) { b[i] = xnew[i]; }

                switch (bc_type)
                {
                    case BoundaryConditions.Dirichlet:
                        b[0] = Lbc(t * dt);
                        b[n - 1] = Rbc(t * dt);
                        break;
                    case BoundaryConditions.Neumann:
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                }

            }
            for (int i = 0; i < n; i++) { cfinal[i] = b[i]; }
        }
        /// <summary>
        /// Method for solving the 1D diffusion equation using the implicit Crank-Nicoloson algorithm
        /// </summary>
        /// <param name="Lbc">Function for the left-side boundary condition for composition</param>
        /// <param name="Rbc">Function for the right-side boundary condition for composition</param>
        /// <param name="I0">Function for the initial condition for composition</param>
        /// <param name="n_steps">Number of time-steps</param>
        public void OneD_CN(BoundaryCondition_Del Lbc, BoundaryCondition_Del Rbc, InitialCondition_Del I0, SourceTerm_Del g, int n_steps)
        {
            int n = C_Initial.GetRVectorSize;
            RVector xnew = new(n);
            RVector xold = new(n);
            RVector b = new(n);
            RVector bj = new(n);
            RVector bj1 = new(n);

            // Setup the A matrix
            double nu = D * dt / Math.Pow(dx, 2);
            double off_d_val, diag_val;
            off_d_val = -nu / 2;
            diag_val = 1.0 + nu;
            TridiagonalMatrix A = new(n, diag_val, off_d_val, off_d_val);
            switch (bc_type)
            {
                case BoundaryConditions.Dirichlet:
                    A[0, 0] = 1.0;
                    A[0, 1] = 0.0;
                    A[n - 1, n - 1] = 1.0;
                    A[n - 1, n - 2] = 0.0;
                    break;
                case BoundaryConditions.Neumann:
                    A[0, 1] = -nu;
                    A[n - 1, n - 2] = -nu;
                    break;
                case BoundaryConditions.Mixed:
                    break;
            }

            // Setup the B matrix
            off_d_val = nu / 2.0;
            diag_val = 1.0 - nu;
            TridiagonalMatrix B = new(n, diag_val, off_d_val, off_d_val);
            switch (bc_type)
            {
                case BoundaryConditions.Dirichlet:
                    B[0, 0] = 1.0;
                    B[0, 1] = 0.0;
                    B[n - 1, n - 1] = 1.0;
                    B[n - 1, n - 2] = 0.0;

                    // Define the initial x vector      
                    for (int i = 0; i < n; i++)
                    {
                        if (i == 0) { C_Initial[i] = Lbc(0.0) + I0(i * dx); }
                        else if (i == n - 1) { C_Initial[i] = Rbc(0.0) + I0(i * dx); }
                        else { C_Initial[i] = I0(i * dx); }
                    }
                    xold = C_Initial; // + g(position, 0.0, C_Initial)
                    break;
                case BoundaryConditions.Neumann:
                    B[0, 1] = nu;
                    B[n - 1, n - 2] = nu;
                    for (int i = 0; i < n; i++) { C_Initial[i] = I0(i * dx); xold[i] = C_Initial[i]; }
                    break;
                case BoundaryConditions.Mixed:
                    break;
            }

            // Time evolution
            for (int t = 0; t < n_steps; t++)
            {
                b = B.Dot(xold); //+ g(position, t * dt, xold)

                switch (bc_type)
                {
                    case BoundaryConditions.Dirichlet:
                        bj1[0] = Lbc((t + 1) * dt);
                        bj1[n - 1] = Rbc((t + 1) * dt);
                        bj[0] = Lbc(t * dt);
                        bj[n - 1] = Rbc(t * dt);
                        break;
                    case BoundaryConditions.Neumann:
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                }
                //xnew = TridiagonalMatrix.Jacobi_Solver(A, x_init_guess, b);
                //xnew = TridiagonalMatrix.GaussSeidel(A, x_init_guess, b);
                //xnew = TridiagonalMatrix.SOR(A, b + bj1 + bj);
                xnew = TridiagonalMatrix.Thomas_Algorithm(A, b + bj1 + bj);
                for (int i = 1; i < n - 1; i++) { xold[i] = xnew[i]; }
            }
            cfinal = xnew;
        }

        /// <summary>
        /// Method for outputting composition and position data to a csv file
        /// </summary>
        /// <param name="of">Filename</param>
        /// <param name="x">Position</param>
        /// <param name="c">Composition</param>
        public static void FileWriteData_CSV(string of, RVector x, RVector c)
        {
            string check_dir = Path.GetDirectoryName(of);
            string owd;

            if (check_dir == null)
            {
                string cwd = Directory.GetCurrentDirectory();
                owd = Path.Combine(cwd, of);
            }
            else
            {
                owd = of;
            }
            if (File.Exists(owd)) { File.Delete(owd); }

            FileStream fS = new(owd, FileMode.OpenOrCreate);
            StreamWriter sW = new(fS);

            string header = "x,c";
            sW.WriteLine(header);
            int nvals = x.GetRVectorSize;
            if (nvals > 0)
            {
                for (int i = 0; i < nvals; i++)
                {
                    string line = x[i].ToString() + "," + c[i].ToString();
                    sW.WriteLine(line);
                }
            }
            else
            {
                throw new Exception("No data available to write to the file!");
            }

            sW.Close();
            fS.Close();

        }

        // Private methods
        /// <summary>
        /// Converts the enumerated boundary condition to a string
        /// </summary>
        /// <param name="bc">variable of type BoundaryCondition</param>
        /// <returns></returns>
        private static string ConvertEnumBCToString(BoundaryConditions bc)
        {
            string result = bc switch
            {
                BoundaryConditions.Dirichlet => "Dirichlet",
                BoundaryConditions.Neumann => "Neumann",
                BoundaryConditions.Mixed => "Mixed",
                BoundaryConditions.Unknown => "Unknown",
                _ => "Unknown",
            };
            return result;
        }
        /// <summary>
        /// Method to convert a string value to one of the enumerated boundary conditions
        /// </summary>
        /// <param name="value">string variable to be converted to the boundary condition enum</param>
        /// <returns></returns>
        private static BoundaryConditions ConvertStringToEnumBC(string value)
        {
            BoundaryConditions bc = value switch
            {
                "Dirichlet" => BoundaryConditions.Dirichlet,
                "Neumann" => BoundaryConditions.Neumann,
                "Mixed" => BoundaryConditions.Mixed,
                _ => BoundaryConditions.Unknown,
            };
            return bc;
        }
    }
}
