using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Diffusion2D_Library
{
    /// <summary>
    /// This class simulates diffusion in a 2D region with a non-constant diffusion coefficient 
    /// </summary>
   public class DiffusionSimulators_2D_MatrixD
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
        /// Delegate for methods for boundary condition functions
        /// </summary>
        /// <param name="time"></param>
        /// <returns>A vector storing the composition on a boundary of a 2D field</returns>
        //public delegate RVector BoundaryCondition_Del(double time, int n);
        public delegate RVector BoundaryCondition_Del(double t, RVector SlidingSide, double FixedValue);
        /// <summary>
        /// Delegate for assigning the initial composition field in a 2D region
        /// </summary>
        /// <param name="c"></param>
        /// <returns>A matrix of composition values</returns>
        public delegate RMatrix InitialCondition_Del(RMatrix x, RMatrix y);
        /// <summary>
        /// Delegate for calculating sources/sinks for reactions occurring in a 2D composition field
        /// </summary>
        /// <param name="position"></param>
        /// <param name="time"></param>
        /// <param name="composition"></param>
        /// <returns>A matrix of composition values</returns>
        public delegate RMatrix SourceTerm_Del(RMatrix xposition, RMatrix yposition, double time, RMatrix DiffCoeff, RMatrix composition);
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
        /// <summary>
        /// Enumeration that specifies the type of boundary condition to be applied to a boundary
        /// </summary>
        /// 
        public enum BoundaryConditions
        {
            Dirichlet,
            Neumann,
            Mixed,
            Unknown
        }
        public enum BoundingBox
        {
            left,
            right,
            top,
            bottom
        }

        public enum Mode { Quiet, Verbose }
        public struct Boundary
        {
            public BoundingBox BoundaryLocation;
            public BoundaryConditions TypeBC;
            public BoundaryCondition_Del BoundaryFunction;
            public RVector PositionVaries;
            public double PositionFixed;
        }

        // Fields
        private string b_filename;
        private const string suffix = ".csv";

        private readonly double dt;
        private readonly double dx;
        private readonly double dy;
        private CompositionField2D CF_2D;
        private Boundary[] border;
        private Mode text_mode;
        private TridiagonalMatrix[] A;
        private TridiagonalMatrix[] B_row;
        private TridiagonalMatrix[] B_col;
        private RMatrix diffusivity;

        private readonly InitialCondition_Del I0;
        private readonly SourceTerm_MatrixD_Del gxt;

        // Properties
        public RMatrix C_Initial
        {
            get { return CF_2D.InitialComposition; }
            set
            {
                if (value.GetnCols > 0 && value.GetnRows > 0) { CF_2D.InitialComposition = value; }
            }
        }
        public RMatrix C_Final
        {
            get { return CF_2D.FinalComposition; }
            set { if (value.GetnCols > 0 && value.GetnRows > 0) { CF_2D.FinalComposition = value; } }
        }
        public RMatrix DiffusionCoefficientField
        {
            get { return CF_2D.DiffusionCoefficient; }
            set { if (value.GetnCols > 0 && value.GetnRows > 0) { CF_2D.DiffusionCoefficient = value; } }
        }
        public RMatrix X
        {
            get { return CF_2D.xposition_matrix; }
        }
        public RMatrix Y
        {
            get { return CF_2D.yposition_matrix; }
        }
        public Boundary[] BCs
        {
            get { return border; }
            set { border = value; }
        }
        public Mode Chat_mode
        {
            get { return text_mode; }
            set { text_mode = value; }
        }
        public string Base_filename
        {
            get { return b_filename;}
            set { b_filename = value; }
        }
        public RMatrix D
        {
            get { return diffusivity; }
            set { diffusivity = value; }
        }

        // Constructors
        public DiffusionSimulators_2D_MatrixD(RMatrix D, double dx, double dy, int nx, int ny, double dt, int nt,
            string[] Boundary_Conditions, BoundaryCondition_Del[] bc_s, InitialCondition_Del I0, SourceTerm_MatrixD_Del g)
        {
            this.dx = dx;
            this.dy = dy;
            this.dt = dt;
            this.diffusivity = D;
            CF_2D = new CompositionField2D(ny, nx);

            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    CF_2D.xposition_matrix[j, i] = i * dx;
                    CF_2D.yposition_matrix[j, i] = j * dy;
                    CF_2D.DiffusionCoefficient[j, i] = D[j, i];
                }
            }

            C_Initial = I0(CF_2D.xposition_matrix, CF_2D.yposition_matrix);
            this.I0 = I0;

            A = new TridiagonalMatrix[nx];
            B_row = new TridiagonalMatrix[nx];
            B_col = new TridiagonalMatrix[nx];

            double nu0 = dt / (2 * Math.Pow(dx, 2));

            RVector nu = new(nx);
            RVector off_d_val_l = new(nx - 1);
            RVector off_d_val_u = new(nx - 1);
            RVector main = new(nx);
            double D_minus, D_plus;
            for (int i = 0; i < nx; i++)
            {
                RVector D_row = D.GetRowVector(i);
                //nu = nu0 * D_row;
                //main = 1 - (2 * nu);
                // Working this for diffusion heterogeneous media
                //================================================
                // A Matrix first
                //================================================
                D_plus = 0.5 * (D_row[1] - D_row[0]);
                main[0] = main[0] = 1.0 - ((D_plus * nu0));

                D_minus = 0.5 * (D_row[nx - 2] - D_row[0]);
                main[nx - 1] = 1.0 - ((D_minus * nu0));

                for (int j = 1; j < nx - 2; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] - D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] - D_row[j]);
                    main[j] = 1.0 - ((D_minus * nu0) + (D_plus * nu0));
                }

                //off_d_val_l[0] = 0.5 * (D_row[0] - D_row[1]) * nu0; 
                //off_d_val_l[nx - 2] = nu[nx - 1];
                for (int j = 1; j < nx; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] - D_row[j]);
                    off_d_val_l[j - 1] = nu0 * D_minus;
                }

                //off_d_val_u[0] = nu[0]; off_d_val_u[nx - 2] = nu[nx - 1];
                for (int j = 0; j < nx - 1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] - D_row[j]);
                    off_d_val_u[j] = nu0 * D_plus;
                }

                A[i] = new TridiagonalMatrix(main, off_d_val_l, off_d_val_u);
                ////================================================
                //// The B-row Matrix next
                ////================================================
                //main = 1 + (2 * nu);

                D_plus = 0.5 * (D_row[1] - D_row[0]);
                main[0] = main[0] = 1.0 + ((D_plus * nu0));

                D_minus = 0.5 * (D_row[nx - 2] - D_row[0]);
                main[nx - 1] = 1.0 + ((D_minus * nu0));

                for (int j = 1; j < nx - 2; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] - D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] - D_row[j]);
                    main[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = 1; j < nx; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] - D_row[j]);
                    off_d_val_l[j - 1] = -nu0 * D_minus;
                }

                //off_d_val_u[0] = nu[0]; off_d_val_u[nx - 2] = nu[nx - 1];
                for (int j = 0; j < nx - 1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] - D_row[j]);
                    off_d_val_u[j] = -nu0 * D_plus;
                }

                B_row[i] = new TridiagonalMatrix(main, off_d_val_l, off_d_val_u);
                ////================================================
                //// The B-column Matrix last
                ////================================================
                RVector D_col = D.GetColVector(i);

                D_plus = 0.5 * (D_col[1] - D_col[0]);
                main[0] = main[0] = 1.0 + ((D_plus * nu0));

                D_minus = 0.5 * (D_col[nx - 2] - D_col[0]);
                main[nx - 1] = 1.0 + ((D_minus * nu0));

                for (int j = 1; j < nx - 2; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] - D_col[j]);
                    D_plus = 0.5 * (D_col[j + 1] - D_col[j]);
                    main[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = 1; j < nx; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] - D_col[j]);
                    off_d_val_l[j - 1] = -nu0 * D_minus;
                }

                //off_d_val_u[0] = nu[0]; off_d_val_u[nx - 2] = nu[nx - 1];
                for (int j = 0; j < nx - 1; j++)
                {
                    D_plus = 0.5 * (D_col[j + 1] - D_col[j]);
                    off_d_val_u[j] = -nu0 * D_plus;
                }

                ////for (int j = 0; j < nx - 1; j++) { off_d_val_l[j] = nu[j]; }
                ////for (int j = 1; j < nx; j++) { off_d_val_u[j - 1] = nu[j]; }

                B_col[i] = new TridiagonalMatrix(main, -off_d_val_l, -off_d_val_u);
            }

            C_Initial = I0(CF_2D.xposition_matrix, CF_2D.yposition_matrix);
            this.I0 = I0;

            int num_bounds = Boundary_Conditions.Length;
            border = new Boundary[num_bounds];
            RVector C0;
            for (int i = 0; i < num_bounds; i++)
            {
                border[i] = new Boundary
                {
                    BoundaryLocation = i switch
                    {
                        0 => BoundingBox.top,
                        1 => BoundingBox.right,
                        2 => BoundingBox.left,
                        3 => BoundingBox.bottom,
                        _ => BoundingBox.bottom,
                    },
                    TypeBC = ConvertStringToEnumBC(Boundary_Conditions[i]),
                    BoundaryFunction = bc_s[i]
                };
                switch (border[i].BoundaryLocation)
                {
                    case BoundingBox.top:
                        if (border[i].TypeBC == BoundaryConditions.Dirichlet)
                        {
                            border[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border[i].PositionFixed = Y[X.GetnRows - 1, 0];
                            C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(ny - 1) + C0;
                            C_Initial.ReplaceRow(C0, ny - 1);
                        }
                        break;
                    case BoundingBox.right:
                        if (border[i].TypeBC == BoundaryConditions.Dirichlet)
                        {
                            border[i].PositionVaries = Y.GetColVector(0);
                            border[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 1) + C0;
                            C_Initial.ReplaceCol(C0, nx - 1);
                        }
                        break;
                    case BoundingBox.left:
                        if (border[i].TypeBC == BoundaryConditions.Dirichlet)
                        {
                            border[i].PositionVaries = Y.GetColVector(0);
                            border[i].PositionFixed = X[0, 0];
                            C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(0) + C0;
                            C_Initial.ReplaceCol(C0, 0);
                        }
                        break;
                    case BoundingBox.bottom:
                        if (border[i].TypeBC == BoundaryConditions.Dirichlet)
                        {
                            border[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border[i].PositionFixed = Y[0, 0];
                            C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(0) + C0;
                            C_Initial.ReplaceRow(C0, 0);
                        }
                        break;
                }
            }

            gxt = g;
            Chat_mode = Mode.Quiet;
        }

        public DiffusionSimulators_2D_MatrixD(RMatrix D, double dx, double dy, int nx, int ny, double dt, int nt,
    string[] Boundary_Conditions, BoundaryCondition_Del[] bc_s, InitialCondition_Del I0, SourceTerm_MatrixD_Del g, Mode Tmode, string base_filename)
        {
            this.dx = dx;
            this.dy = dy;
            this.dt = dt;
            this.diffusivity = D;
            CF_2D = new CompositionField2D(ny, nx);

            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    CF_2D.xposition_matrix[j, i] = i * dx;
                    CF_2D.yposition_matrix[j, i] = j * dy;
                    CF_2D.DiffusionCoefficient[j, i] = D[j, i];
                }
            }

            C_Initial = I0(CF_2D.xposition_matrix, CF_2D.yposition_matrix);
            this.I0 = I0;

            A = new TridiagonalMatrix[nx];
            B_row = new TridiagonalMatrix[nx];
            B_col = new TridiagonalMatrix[nx];

            double nu0 = dt / (2 * Math.Pow(dx, 2));

            RVector nu = new(nx);
            RVector off_d_val_l = new(nx - 1);
            RVector off_d_val_u = new(nx - 1);
            RVector main = new(nx);
            double D_minus, D_plus;
            for (int i = 0; i < nx; i++)
            {
                RVector D_row = D.GetRowVector(i);
                //nu = nu0 * D_row;
                //main = 1 - (2 * nu);
                // Working this for diffusion heterogeneous media
                //================================================
                // A Matrix first
                //================================================
                D_plus = 0.5 * (D_row[1] - D_row[0]);
                main[0] = main[0] = 1.0 - ((D_plus * nu0));

                D_minus = 0.5 * (D_row[nx - 2] - D_row[0]);
                main[nx - 1] = 1.0 - ((D_minus * nu0));
                
                for (int j = 1; j < nx - 2; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] - D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] - D_row[j]);
                    main[j] = 1.0 - ((D_minus * nu0) + (D_plus * nu0));
                }

                //off_d_val_l[0] = 0.5 * (D_row[0] - D_row[1]) * nu0; 
                //off_d_val_l[nx - 2] = nu[nx - 1];
                for (int j = 1; j < nx; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] - D_row[j]);
                    off_d_val_l[j - 1] = nu0 * D_minus;
                }

                //off_d_val_u[0] = nu[0]; off_d_val_u[nx - 2] = nu[nx - 1];
                for (int j = 0; j < nx - 1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] - D_row[j]);
                    off_d_val_u[j] = nu0 * D_plus;
                }

                A[i] = new TridiagonalMatrix(main, off_d_val_l, off_d_val_u);
                ////================================================
                //// The B-row Matrix next
                ////================================================
                //main = 1 + (2 * nu);

                D_plus = 0.5 * (D_row[1] - D_row[0]);
                main[0] = main[0] = 1.0 + ((D_plus * nu0));

                D_minus = 0.5 * (D_row[nx - 2] - D_row[0]);
                main[nx - 1] = 1.0 + ((D_minus * nu0));

                for (int j = 1; j < nx - 2; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] - D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] - D_row[j]);
                    main[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = 1; j < nx; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] - D_row[j]);
                    off_d_val_l[j - 1] = -nu0 * D_minus;
                }

                //off_d_val_u[0] = nu[0]; off_d_val_u[nx - 2] = nu[nx - 1];
                for (int j = 0; j < nx - 1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] - D_row[j]);
                    off_d_val_u[j] = -nu0 * D_plus;
                }

                B_row[i] = new TridiagonalMatrix(main, off_d_val_l, off_d_val_u);
                ////================================================
                //// The B-column Matrix last
                ////================================================
                RVector D_col = D.GetColVector(i);
               
                D_plus = 0.5 * (D_col[1] - D_col[0]);
                main[0] = main[0] = 1.0 + ((D_plus * nu0));

                D_minus = 0.5 * (D_col[nx - 2] - D_col[0]);
                main[nx - 1] = 1.0 + ((D_minus * nu0));

                for (int j = 1; j < nx - 2; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] - D_col[j]);
                    D_plus = 0.5 * (D_col[j + 1] - D_col[j]);
                    main[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = 1; j < nx; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] - D_col[j]);
                    off_d_val_l[j - 1] = -nu0 * D_minus;
                }

                //off_d_val_u[0] = nu[0]; off_d_val_u[nx - 2] = nu[nx - 1];
                for (int j = 0; j < nx - 1; j++)
                {
                    D_plus = 0.5 * (D_col[j + 1] - D_col[j]);
                    off_d_val_u[j] = -nu0 * D_plus;
                }

                ////for (int j = 0; j < nx - 1; j++) { off_d_val_l[j] = nu[j]; }
                ////for (int j = 1; j < nx; j++) { off_d_val_u[j - 1] = nu[j]; }

                B_col[i] = new TridiagonalMatrix(main, -off_d_val_l, -off_d_val_u);
            }

            C_Initial = I0(CF_2D.xposition_matrix, CF_2D.yposition_matrix);
            this.I0 = I0;

            int num_bounds = Boundary_Conditions.Length;
            border = new Boundary[num_bounds];
            RVector C0;
            for (int i = 0; i < num_bounds; i++)
            {
                border[i] = new Boundary
                {
                    BoundaryLocation = i switch
                    {
                        0 => BoundingBox.top,
                        1 => BoundingBox.right,
                        2 => BoundingBox.left,
                        3 => BoundingBox.bottom,
                        _ => BoundingBox.bottom,
                    },
                    TypeBC = ConvertStringToEnumBC(Boundary_Conditions[i]),
                    BoundaryFunction = bc_s[i]
                };
                switch (border[i].BoundaryLocation)
                {
                    case BoundingBox.top:
                        if (border[i].TypeBC == BoundaryConditions.Dirichlet)
                        {
                            border[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border[i].PositionFixed = Y[X.GetnRows - 1, 0];
                            C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(ny - 1) + C0;
                            C_Initial.ReplaceRow(C0, ny - 1);
                        }
                        break;
                    case BoundingBox.right:
                        if (border[i].TypeBC == BoundaryConditions.Dirichlet)
                        {
                            border[i].PositionVaries = Y.GetColVector(0);
                            border[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 1) + C0;
                            C_Initial.ReplaceCol(C0, nx - 1);
                        }
                        break;
                    case BoundingBox.left:
                        if (border[i].TypeBC == BoundaryConditions.Dirichlet)
                        {
                            border[i].PositionVaries = Y.GetColVector(0);
                            border[i].PositionFixed = X[0, 0];
                            C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(0) + C0;
                            C_Initial.ReplaceCol(C0, 0);
                        }
                        break;
                    case BoundingBox.bottom:
                        if (border[i].TypeBC == BoundaryConditions.Dirichlet)
                        {
                            border[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border[i].PositionFixed = Y[0, 0];
                            C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(0) + C0;
                            C_Initial.ReplaceRow(C0, 0);
                        }
                        break;
                }
            }

            gxt = g;
            Chat_mode = Tmode;
            b_filename = base_filename;
        }

        // Solvers
        /// <summary>
        /// Method for solving the 2D diffusion equation using the Alternating Direction Implicit algorithm
        /// </summary>
        /// <param name="n_time_steps"></param>
        public void TwoD_ADI(int n_time_steps)
        {
            int nrows = CF_2D.InitialComposition.GetnRows;
            int ncols = CF_2D.InitialComposition.GetnCols;

            RMatrix C_Ex = new(nrows, ncols);
            RMatrix C_Im1 = new(nrows, ncols);
            RMatrix C_Im2 = new(nrows, ncols);

            RMatrix fn = new(nrows, ncols);
            RMatrix f0 = new(nrows, ncols);
            RMatrix f12 = new(nrows, ncols);

            RVector CL = new(nrows);
            RVector CR = new(nrows);
            RVector CT = new(ncols);
            RVector CB = new(ncols);

            double nu0 = dt / (2 * Math.Pow(dx, 2));

            //// Define the A matrix for the explicit steps
            //double off_d_val = nu;
            //double diag_val = 1 - (2 * nu);
            //TridiagonalMatrix A = new(ncols, diag_val, off_d_val, off_d_val);

            //// Define the B matrices for the implicit steps
            //off_d_val = -nu;
            //diag_val = 1 + (2 * nu);
            //TridiagonalMatrix B = new(nrows, diag_val, off_d_val, off_d_val);

            // Time evolution           
            for (int t = 0; t < n_time_steps; t++)
            {
                if (Chat_mode == Mode.Verbose)
                {
                    if (t % 10 == 0) { Console.WriteLine(t * dt); }
                }

                //if (t == n_time_steps - 1)
                //{
                //    Console.WriteLine(t);
                //}

                // ===================
                // Explicit time-step
                // ===================
                // 0 = top, 1 = right, 2 = left, 3 = bottom
                switch (BCs[1].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        CR = BCs[1].BoundaryFunction(t * dt, BCs[1].PositionVaries, BCs[1].PositionFixed);
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        double cUpper1, cUpper2, cLower1, cLower2;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetColVector(ncols - 1);
                            C2 = C_Initial.GetColVector(ncols - 2);

                            cLower1 = C_Initial[0, ncols - 1];
                            cLower2 = C_Initial[1, ncols - 2];
                            cUpper1 = C_Initial[nrows - 1, ncols - 1];
                            cUpper2 = C_Initial[nrows - 2, ncols - 2];
                        }
                        else
                        {
                            C1 = C_Im2.GetColVector(ncols - 1);
                            C2 = C_Im2.GetColVector(ncols - 2);

                            cLower1 = C_Im2[0, ncols - 1];
                            cLower2 = C_Im2[1, ncols - 2];
                            cUpper1 = C_Im2[nrows - 1, ncols - 1];
                            cUpper2 = C_Im2[nrows - 2, ncols - 2];
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        //CR = ((-2 * nu0) * C1) + ((2 * nu0) * C2); // + C3
                        //CR[0] = ((-2 * nu0) * cLower1) + ((2 * nu0) * cLower2);
                        //CR[ncols - 1] = ((-2 * nu0) * cUpper1) + ((2 * nu0) * cUpper2);

                        CR = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(ncols - 1), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(ncols - 2), C2); // + C3
                        CR[0] = ((-2 * nu0 * CF_2D.DiffusionCoefficient[0, ncols - 1]) * cLower1) + ((2 * nu0 * CF_2D.DiffusionCoefficient[1, ncols - 2]) * cLower2);
                        CR[ncols - 1] = ((-2 * nu0 * CF_2D.DiffusionCoefficient[nrows - 1, ncols - 1]) * cUpper1) + ((2 * nu0 * CF_2D.DiffusionCoefficient[nrows - 2, ncols - 2]) * cUpper2);
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                switch (BCs[2].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        CL = BCs[2].BoundaryFunction(t * dt, BCs[2].PositionVaries, BCs[2].PositionFixed);
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        double cUpper1, cUpper2, cLower1, cLower2;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetColVector(0);
                            C2 = C_Initial.GetColVector(1);

                            cLower1 = C_Initial[0, 0];
                            cLower2 = C_Initial[1, 1];
                            cUpper1 = C_Initial[nrows - 1, 0];
                            cUpper2 = C_Initial[nrows - 2, 1];
                        }
                        else
                        {
                            C1 = C_Im2.GetColVector(0);
                            C2 = C_Im2.GetColVector(1);

                            cLower1 = C_Im2[0, 0];
                            cLower2 = C_Im2[1, 1];
                            cUpper1 = C_Im2[nrows - 1, 0];
                            cUpper2 = C_Im2[nrows - 2, 1];
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);

                        CL = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(0),C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(1), C2); // + C3
                        CL[0] = ((-2 * nu0 * CF_2D.DiffusionCoefficient[0,0]) * cLower1) + ((2 * nu0 * CF_2D.DiffusionCoefficient[1, 1]) * cLower2);
                        CL[ncols - 1] = ((-2 * nu0 * CF_2D.DiffusionCoefficient[nrows - 1, 0]) * cUpper1) + ((2 * nu0 * CF_2D.DiffusionCoefficient[nrows - 2, 1]) * cUpper2);
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                for (int i = 1; i < nrows-1; i++)
                {
                    RVector xold;
                    if (t == 0) { xold = C_Initial.GetRowVector(i); }
                    else { xold = C_Im2.GetRowVector(i); } // C_Im2.GetRowVector(i); }

                    RVector v1 = A[i].Dot(xold);
                    C_Ex.ReplaceRow(v1, i);
                    C_Ex[i, 0] = CL[i]; //nu *
                    C_Ex[i, ncols - 1] = CR[i]; //nu *
                }

                // ===================
                // ===================
                // One-half implicit time-step
                // ===================
                // Source terms
                double t12 = (t + 0.5) * dt;
                if (t == 0) { fn = gxt(X, Y, t12, CF_2D.DiffusionCoefficient, C_Initial); }
                else { fn = gxt(X, Y, t12, CF_2D.DiffusionCoefficient, C_Im2); }

                f12 = (dt / 2.0) * fn;

                // BCs
                switch (BCs[0].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        RVector gn = BCs[0].BoundaryFunction((t + 1) * dt, BCs[0].PositionVaries, BCs[0].PositionFixed);
                        RVector g0 = BCs[0].BoundaryFunction(t * dt, BCs[0].PositionVaries, BCs[0].PositionFixed);
                        CT = ((0.5 * B_row[nrows-1].Dot(gn)) + (0.5 * A[nrows-1].Dot(g0))); ////((0.5 * Bm.Dot(gn)) + (0.5 * Bp.Dot(g0))); //nu * 
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetRowVector(nrows - 1);
                            C2 = C_Initial.GetRowVector(nrows - 2);
                        }
                        else
                        {
                            C1 = C_Im2.GetRowVector(nrows - 1);
                            C2 = C_Im2.GetRowVector(nrows - 2);
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        CT = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetRowVector(nrows-1), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetRowVector(nrows-2), C2); // + C3
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                switch (BCs[3].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        RVector gn = BCs[3].BoundaryFunction((t + 1) * dt, BCs[3].PositionVaries, BCs[3].PositionFixed);
                        RVector g0 = BCs[3].BoundaryFunction(t * dt, BCs[3].PositionVaries, BCs[3].PositionFixed);
                        CB = ((0.5 * B_row[0].Dot(gn)) + (0.5 * A[0].Dot(g0))); //((0.5 * Bm.Dot(gn)) + (0.5 * Bp.Dot(g0))); //nu * 
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetRowVector(0);
                            C2 = C_Initial.GetRowVector(1);
                        }
                        else
                        {
                            C1 = C_Im2.GetRowVector(0);
                            C2 = C_Im2.GetRowVector(1);
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        CB = RVector.Product(-2 * nu0*CF_2D.DiffusionCoefficient.GetRowVector(0), C1) + RVector.Product(2 * nu0*CF_2D.DiffusionCoefficient.GetRowVector(1), C2); // + C3
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                for (int j = 1; j < ncols-1; j++)
                {
                    RVector v1 = C_Ex.GetColVector(j);
                    v1[0] = CB[j]; //nu * 
                    v1[ncols - 1] = CT[j]; //nu * 

                    RVector f12s = f12.GetColVector(j); // 

                    RVector u12 = TridiagonalMatrix.Thomas_Algorithm(B_row[j], v1 + f12s);
                    C_Im1.ReplaceCol(u12, j);
                }

                // ===================
                // ===================
                // Full implicit time-step
                // ===================
                switch (BCs[1].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        CR = BCs[1].BoundaryFunction((t + 1) * dt, BCs[1].PositionVaries, BCs[1].PositionFixed);
                        CT = BCs[0].BoundaryFunction((t + 1) * dt, BCs[0].PositionVaries, BCs[0].PositionFixed);
                        RVector ctab = C_Initial.GetRowVector(0);
                        C_Im2.ReplaceRow(CT, nrows - 1);
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetColVector(ncols - 1);
                            C2 = C_Initial.GetColVector(ncols - 2);
                        }
                        else
                        {
                            C1 = C_Im2.GetColVector(ncols - 1);
                            C2 = C_Im2.GetColVector(ncols - 2);
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        CR = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(ncols-1), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(ncols-2), C2); // + C3
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                switch (BCs[2].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        CL = BCs[2].BoundaryFunction((t + 1) * dt, BCs[2].PositionVaries, BCs[2].PositionFixed);
                        CB = BCs[3].BoundaryFunction((t + 1) * dt, BCs[3].PositionVaries, BCs[3].PositionFixed);

                        RVector ctab = C_Initial.GetRowVector(nrows - 1);

                        C_Im2.ReplaceRow(CB, 0);
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetColVector(0);
                            C2 = C_Initial.GetColVector(1);
                        }
                        else
                        {
                            C1 = C_Im2.GetColVector(0);
                            C2 = C_Im2.GetColVector(1);
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        CL = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(0), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(1), C2); // + C3
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                for (int k = 1; k < nrows - 1; k++)
                {
                    RVector v1 = C_Ex.GetRowVector(k);
                    RVector u12 = C_Im1.GetRowVector(k);
                    RVector b = (2 * u12) - v1;
                    b[0] = CL[k]; //nu * 
                    b[ncols - 1] = CR[k]; //nu * 

                    RVector u1 = TridiagonalMatrix.Thomas_Algorithm(B_row[k], b);
                    u1[0] = CL[k];  //nu * 
                    u1[ncols - 1] = CR[k]; //nu * 
                    C_Im2.ReplaceRow(u1, k);
                }
                // ===================

                if (Chat_mode == Mode.Verbose) { if (t == n_time_steps - 1) { Console.WriteLine(t * dt); } }

            }

            // Setup the return composition field
            C_Final = C_Im2;
        }

        /// <summary>
        /// Method for solving the 2D diffusion equation using the Alternating Direction Implicit algorithm.  Will output time-steps elapsed.
        /// </summary>
        /// <param name="n_time_steps"></param>
        /// <param name="output_interval"></param>
        public void TwoD_ADI(int n_time_steps, int output_interval)
        {
            int nrows = CF_2D.InitialComposition.GetnRows;
            int ncols = CF_2D.InitialComposition.GetnCols;

            RMatrix C_Ex = new(nrows, ncols);
            RMatrix C_Im1 = new(nrows, ncols);
            RMatrix C_Im2 = new(nrows, ncols);

            RMatrix fn = new(nrows, ncols);
            RMatrix f0 = new(nrows, ncols);
            RMatrix f12 = new(nrows, ncols);

            RVector CL = new(nrows);
            RVector CR = new(nrows);
            RVector CT = new(ncols);
            RVector CB = new(ncols);

            double nu0 = dt / (2 * Math.Pow(dx, 2));
            string full_file_name;
            // Time evolution           
            for (int t = 0; t < n_time_steps; t++)
            {
                if (Chat_mode == Mode.Verbose && output_interval > 0 && Base_filename != null)
                {
                    if (t % output_interval == 0) 
                    {
                        decimal time = (decimal)(t * dt);
                        Console.WriteLine("{0}s have been simulated", time);
                        full_file_name = Base_filename + time.ToString() + suffix;
                        if (File.Exists(full_file_name)) { File.Delete(full_file_name); }
                        if (t == 0) { FileWriteData_CSV(full_file_name, X, Y, C_Initial); }
                        else { FileWriteData_CSV(full_file_name, X, Y, C_Im2); }                        
                    }
                }

                // ===================
                // Explicit time-step
                // ===================
                // 0 = top, 1 = right, 2 = left, 3 = bottom
                switch (BCs[1].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        CR = BCs[1].BoundaryFunction(t * dt, BCs[1].PositionVaries, BCs[1].PositionFixed);
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        double cUpper1, cUpper2, cLower1, cLower2;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetColVector(ncols - 1);
                            C2 = C_Initial.GetColVector(ncols - 2);

                            cLower1 = C_Initial[0, ncols - 1];
                            cLower2 = C_Initial[1, ncols - 2];
                            cUpper1 = C_Initial[nrows - 1, ncols - 1];
                            cUpper2 = C_Initial[nrows - 2, ncols - 2];
                        }
                        else
                        {
                            C1 = C_Im2.GetColVector(ncols - 1);
                            C2 = C_Im2.GetColVector(ncols - 2);

                            cLower1 = C_Im2[0, ncols - 1];
                            cLower2 = C_Im2[1, ncols - 2];
                            cUpper1 = C_Im2[nrows - 1, ncols - 1];
                            cUpper2 = C_Im2[nrows - 2, ncols - 2];
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        //CR = ((-2 * nu0) * C1) + ((2 * nu0) * C2); // + C3
                        //CR[0] = ((-2 * nu0) * cLower1) + ((2 * nu0) * cLower2);
                        //CR[ncols - 1] = ((-2 * nu0) * cUpper1) + ((2 * nu0) * cUpper2);

                        CR = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(ncols - 1), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(ncols - 2), C2); // + C3
                        CR[0] = ((-2 * nu0 * CF_2D.DiffusionCoefficient[0, ncols - 1]) * cLower1) + ((2 * nu0 * CF_2D.DiffusionCoefficient[1, ncols - 2]) * cLower2);
                        CR[ncols - 1] = ((-2 * nu0 * CF_2D.DiffusionCoefficient[nrows - 1, ncols - 1]) * cUpper1) + ((2 * nu0 * CF_2D.DiffusionCoefficient[nrows - 2, ncols - 2]) * cUpper2);
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                switch (BCs[2].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        CL = BCs[2].BoundaryFunction(t * dt, BCs[2].PositionVaries, BCs[2].PositionFixed);
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        double cUpper1, cUpper2, cLower1, cLower2;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetColVector(0);
                            C2 = C_Initial.GetColVector(1);

                            cLower1 = C_Initial[0, 0];
                            cLower2 = C_Initial[1, 1];
                            cUpper1 = C_Initial[nrows - 1, 0];
                            cUpper2 = C_Initial[nrows - 2, 1];
                        }
                        else
                        {
                            C1 = C_Im2.GetColVector(0);
                            C2 = C_Im2.GetColVector(1);

                            cLower1 = C_Im2[0, 0];
                            cLower2 = C_Im2[1, 1];
                            cUpper1 = C_Im2[nrows - 1, 0];
                            cUpper2 = C_Im2[nrows - 2, 1];
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);

                        CL = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(0), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(1), C2); // + C3
                        CL[0] = ((-2 * nu0 * CF_2D.DiffusionCoefficient[0, 0]) * cLower1) + ((2 * nu0 * CF_2D.DiffusionCoefficient[1, 1]) * cLower2);
                        CL[ncols - 1] = ((-2 * nu0 * CF_2D.DiffusionCoefficient[nrows - 1, 0]) * cUpper1) + ((2 * nu0 * CF_2D.DiffusionCoefficient[nrows - 2, 1]) * cUpper2);
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                for (int i = 1; i < nrows - 1; i++)
                {
                    RVector xold;
                    if (t == 0) { xold = C_Initial.GetRowVector(i); }
                    else { xold = C_Im2.GetRowVector(i); } // C_Im2.GetRowVector(i); }

                    RVector v1 = A[i].Dot(xold);
                    C_Ex.ReplaceRow(v1, i);
                    C_Ex[i, 0] = CL[i]; //nu *
                    C_Ex[i, ncols - 1] = CR[i]; //nu *
                }

                // ===================
                // ===================
                // One-half implicit time-step
                // ===================
                // Source terms
                double t12 = (t + 0.5) * dt;
                if (t == 0) { fn = gxt(X, Y, t12, CF_2D.DiffusionCoefficient, C_Initial); }
                else { fn = gxt(X, Y, t12, CF_2D.DiffusionCoefficient, C_Im2); }

                f12 = (dt / 2.0) * fn;

                // BCs
                switch (BCs[0].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        RVector gn = BCs[0].BoundaryFunction((t + 1) * dt, BCs[0].PositionVaries, BCs[0].PositionFixed);
                        RVector g0 = BCs[0].BoundaryFunction(t * dt, BCs[0].PositionVaries, BCs[0].PositionFixed);
                        CT = ((0.5 * B_row[nrows - 1].Dot(gn)) + (0.5 * A[nrows - 1].Dot(g0))); ////((0.5 * Bm.Dot(gn)) + (0.5 * Bp.Dot(g0))); //nu * 
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetRowVector(nrows - 1);
                            C2 = C_Initial.GetRowVector(nrows - 2);
                        }
                        else
                        {
                            C1 = C_Im2.GetRowVector(nrows - 1);
                            C2 = C_Im2.GetRowVector(nrows - 2);
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        CT = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetRowVector(nrows - 1), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetRowVector(nrows - 2), C2); // + C3
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                switch (BCs[3].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        RVector gn = BCs[3].BoundaryFunction((t + 1) * dt, BCs[3].PositionVaries, BCs[3].PositionFixed);
                        RVector g0 = BCs[3].BoundaryFunction(t * dt, BCs[3].PositionVaries, BCs[3].PositionFixed);
                        CB = ((0.5 * B_row[0].Dot(gn)) + (0.5 * A[0].Dot(g0))); //((0.5 * Bm.Dot(gn)) + (0.5 * Bp.Dot(g0))); //nu * 
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetRowVector(0);
                            C2 = C_Initial.GetRowVector(1);
                        }
                        else
                        {
                            C1 = C_Im2.GetRowVector(0);
                            C2 = C_Im2.GetRowVector(1);
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        CB = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetRowVector(0), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetRowVector(1), C2); // + C3
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                for (int j = 1; j < ncols - 1; j++)
                {
                    RVector v1 = C_Ex.GetColVector(j);
                    v1[0] = CB[j]; //nu * 
                    v1[ncols - 1] = CT[j]; //nu * 

                    RVector f12s = f12.GetColVector(j); // 

                    switch (BCs[0].TypeBC)
                    {
                        case BoundaryConditions.Dirichlet:
                            B_col[j][0, 0] = 1.0;
                            B_col[j][0, 1] = 0.0;
                            break;
                        default:
                            B_col[j][0, 0] = 1.0;
                            B_col[j][0, 1] = 0.0;
                            break;
                    }
                    switch (BCs[3].TypeBC)
                    {
                        case BoundaryConditions.Dirichlet:
                            B_col[j][ncols - 1, ncols - 1] = 1.0;
                            B_col[j][ncols - 1, ncols - 2] = 0.0;
                            break;
                        default:
                            B_col[j][ncols - 1, ncols - 1] = 1.0;
                            B_col[j][ncols - 1, ncols - 2] = 0.0;
                            break;
                    }
                    RVector u12 = TridiagonalMatrix.Thomas_Algorithm(B_col[j], v1 + f12s);
                    C_Im1.ReplaceCol(u12, j);
                }

                // ===================
                // ===================
                // Full implicit time-step
                // ===================
                switch (BCs[1].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        CR = BCs[1].BoundaryFunction((t + 1) * dt, BCs[1].PositionVaries, BCs[1].PositionFixed);
                        CT = BCs[0].BoundaryFunction((t + 1) * dt, BCs[0].PositionVaries, BCs[0].PositionFixed);
                        RVector ctab = C_Initial.GetRowVector(0);
                        C_Im2.ReplaceRow(CT, nrows - 1);
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetColVector(ncols - 1);
                            C2 = C_Initial.GetColVector(ncols - 2);
                        }
                        else
                        {
                            C1 = C_Im2.GetColVector(ncols - 1);
                            C2 = C_Im2.GetColVector(ncols - 2);
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        CR = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(ncols - 1), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(ncols - 2), C2); // + C3
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                switch (BCs[2].TypeBC)
                {
                    case BoundaryConditions.Dirichlet:
                        CL = BCs[2].BoundaryFunction((t + 1) * dt, BCs[2].PositionVaries, BCs[2].PositionFixed);
                        CB = BCs[3].BoundaryFunction((t + 1) * dt, BCs[3].PositionVaries, BCs[3].PositionFixed);

                        RVector ctab = C_Initial.GetRowVector(nrows - 1);

                        C_Im2.ReplaceRow(CB, 0);
                        break;
                    case BoundaryConditions.Neumann:
                        RVector C1, C2; //,C3;
                        if (t == 0)
                        {
                            C1 = C_Initial.GetColVector(0);
                            C2 = C_Initial.GetColVector(1);
                        }
                        else
                        {
                            C1 = C_Im2.GetColVector(0);
                            C2 = C_Im2.GetColVector(1);
                        }
                        //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
                        CL = RVector.Product(-2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(0), C1) + RVector.Product(2 * nu0 * CF_2D.DiffusionCoefficient.GetColVector(1), C2); // + C3
                        break;
                    case BoundaryConditions.Mixed:
                        break;
                    case BoundaryConditions.Unknown:
                        break;
                    default:
                        break;
                }
                for (int k = 1; k < nrows - 1; k++)
                {
                    RVector v1 = C_Ex.GetRowVector(k);
                    RVector u12 = C_Im1.GetRowVector(k);
                    RVector b = (2 * u12) - v1;
                    b[0] = CL[k]; //nu * 
                    b[ncols - 1] = CR[k]; //nu * 

                    switch (BCs[0].TypeBC)
                    {
                        case BoundaryConditions.Dirichlet:
                            B_row[k][0, 0] = 1.0;
                            B_row[k][0, 1] = 0.0;
                            break;
                        default:
                            B_row[k][0, 0] = 1.0;
                            B_row[k][0, 1] = 0.0;
                            break;
                    }
                    switch (BCs[3].TypeBC)
                    {
                        case BoundaryConditions.Dirichlet:
                            B_row[k][ncols - 1, ncols - 1] = 1.0;
                            B_row[k][ncols - 1, ncols - 2] = 0.0;
                            break;
                        default:
                            B_row[k][ncols - 1, ncols - 1] = 1.0;
                            B_row[k][ncols - 1, ncols - 2] = 0.0;
                            break;
                    }
                    RVector u1 = TridiagonalMatrix.Thomas_Algorithm(B_row[k], b);
                    u1[0] = CL[k];  //nu * 
                    u1[ncols - 1] = CR[k]; //nu * 
                    C_Im2.ReplaceRow(u1, k);
                }
                // ===================

            }

            // Setup the return composition field
            C_Final = C_Im2;
        }

        /// <summary>
        /// Method for outputting composition and position data to a csv file
        /// </summary>
        /// <param name="of">Filename</param>
        /// <param name="x">Position</param>
        /// <param name="c">Composition</param>
        public static void FileWriteData_CSV(string of, RMatrix x, RMatrix y, RMatrix c)
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

            string header = "x,y,c";
            sW.WriteLine(header);

            int nrows = x.GetnRows;
            int ncols = x.GetnCols;

            if (nrows > 0 && ncols > 0)
            {
                for (int i = 0; i < nrows; i++)
                {
                    for (int j = 0; j < ncols; j++)
                    {
                        string line = x[i, j].ToString() + "," + y[i, j].ToString() + "," + c[i, j].ToString();
                        sW.WriteLine(line);
                    }
                }

            }
            else
            {
                throw new Exception("No data available to write to the file!");
            }

            sW.Close();
            fS.Close();

        }

        // =======================================================================
        // Private methods
        // =======================================================================
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

        // =======================================================================

        //
        // =======================================================================
        // =======================================================================
        // =======================================================================


        //public DiffusionSimulators_2D_MatrixD(double dx, double dy, int nx, int ny, double dt, int nt,
        //    string[] Boundary_Conditions, BoundaryCondition_Del[] bc_s, InitialCondition_Del I0, SourceTerm_Del g, Mode chat)
        //{
        //    this.D = D;
        //    this.dx = dx;
        //    this.dy = dy;
        //    this.dt = dt;
        //    CF_2D = new CompositionField2D(ny, nx);

        //    Parallel.For(0, ny, i =>
        //    {
        //        for (int j = 0; j < nx; j++)
        //        {
        //            CF_2D.xposition_matrix[j, i] = i * dx;
        //            CF_2D.yposition_matrix[j, i] = j * dy;
        //        }
        //    });
        //    C_Initial = I0(CF_2D.xposition_matrix, CF_2D.yposition_matrix);

        //    int num_bounds = Boundary_Conditions.Length;
        //    border = new Boundary[num_bounds];
        //    RVector C0;
        //    for (int i = 0; i < num_bounds; i++)
        //    {
        //        border[i] = new Boundary
        //        {
        //            BoundaryLocation = i switch
        //            {
        //                0 => BoundingBox.top,
        //                1 => BoundingBox.right,
        //                2 => BoundingBox.left,
        //                3 => BoundingBox.bottom,
        //                _ => BoundingBox.bottom,
        //            },
        //            TypeBC = ConvertStringToEnumBC(Boundary_Conditions[i]),
        //            BoundaryFunction = bc_s[i]
        //        };
        //        switch (border[i].BoundaryLocation)
        //        {
        //            case BoundingBox.top:
        //                if (border[i].TypeBC == BoundaryConditions.Dirichlet)
        //                {
        //                    border[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
        //                    border[i].PositionFixed = Y[X.GetnRows - 1, 0];
        //                    C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

        //                    RVector Ctab = C_Initial.GetRowVector(ny - 1) + C0;
        //                    C_Initial.ReplaceRow(C0, ny - 1);
        //                }
        //                break;
        //            case BoundingBox.right:
        //                if (border[i].TypeBC == BoundaryConditions.Dirichlet)
        //                {
        //                    border[i].PositionVaries = Y.GetColVector(0);
        //                    border[i].PositionFixed = X[0, Y.GetnCols - 1];
        //                    C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

        //                    RVector Ctab = C_Initial.GetColVector(nx - 1) + C0;
        //                    C_Initial.ReplaceCol(C0, nx - 1);
        //                }
        //                break;
        //            case BoundingBox.left:
        //                if (border[i].TypeBC == BoundaryConditions.Dirichlet)
        //                {
        //                    border[i].PositionVaries = Y.GetColVector(0);
        //                    border[i].PositionFixed = X[0, 0];
        //                    C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

        //                    RVector Ctab = C_Initial.GetColVector(0) + C0;
        //                    C_Initial.ReplaceCol(C0, 0);
        //                }
        //                break;
        //            case BoundingBox.bottom:
        //                if (border[i].TypeBC == BoundaryConditions.Dirichlet)
        //                {
        //                    border[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
        //                    border[i].PositionFixed = Y[0, 0];
        //                    C0 = border[i].BoundaryFunction(0.0, border[i].PositionVaries, border[i].PositionFixed);

        //                    RVector Ctab = C_Initial.GetRowVector(0) + C0;
        //                    C_Initial.ReplaceRow(C0, 0);
        //                }
        //                break;
        //        }
        //    }
        //    this.I0 = I0;
        //    gxt = g;
        //    Chat_mode = chat;
        //}
        //public DiffusionSimulators_2D_MatrixD(double[] coeffs, int[] n, InitialCondition_Del I0)  //string Boundary_Conditions,
        //{
        //    if (coeffs.Length >= 3)
        //    {
        //        D = coeffs[0];
        //        dx = coeffs[1];
        //        dt = coeffs[2];
        //    }
        //    if (n.Length == 2)
        //    {
        //        CF_2D = new CompositionField2D(n[0], n[1]);
        //    }

        //    //this.Boundary_Conditions = Boundary_Conditions;
        //    C_Initial = I0(CF_2D.xposition_matrix, CF_2D.yposition_matrix);
        //}

    }
}
