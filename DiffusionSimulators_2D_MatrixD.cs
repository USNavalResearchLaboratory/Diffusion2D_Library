using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static Diffusion2D_Library.BoundaryCondition;
using static Diffusion2D_Library.DiffusionSimulators_1D;

namespace Diffusion2D_Library
{
    /// <summary>
    /// This class simulates diffusion in a 2D region with a non-constant diffusion coefficient 
    /// </summary>
   public class DiffusionSimulators_2D_MatrixD : DiffusionSimulator2D
    {
        // Fields
        private readonly TridiagonalMatrix[] A;
        private readonly TridiagonalMatrix[] B_row;
        private readonly TridiagonalMatrix[] B_col;
        private double[] NeumannBCs_L_A;
        private double[] NeumannBCs_R_A;
        

        // Properties
        public RMatrix D
        {
            get => cf_2D.DValues;
            set { if (value.GetnCols > 0 && value.GetnRows > 0) { cf_2D.DValues = value; } }
        }

        // ====================================================================
        // Constructors
        // ====================================================================
        /// <summary>
        /// Constructs a 2D solver for the diffusion equation with a matrix diffusivity
        /// </summary>
        /// <param name="D"></param>
        /// <param name="dx"></param>
        /// <param name="dy"></param>
        /// <param name="nx"></param>
        /// <param name="ny"></param>
        /// <param name="dt"></param>
        /// <param name="nt"></param>
        /// <param name="Boundary_Conditions"></param>
        /// <param name="bc_s"></param>
        /// <param name="I0"></param>
        /// <param name="g"></param>
        /// <param name="Tmode"></param>
        /// <param name="base_filename"></param>
        public DiffusionSimulators_2D_MatrixD(
            RMatrix D, 
            double dx, 
            double dy, 
            int nx, 
            int ny, 
            double dt, 
            int nt, 
            ABoundaryCondition[] Boundary_Conditions, 
            Del_BC_xy[] bc_s, 
            Del_IC_xy I0, 
            Del_Source_MatrixD g, 
            Mode Tmode, 
            string base_filename)
        {
            NX = nx;
            NY = ny;
            this.dx = dx;
            this.dy = dy;
            this.dt = dt;
            cf_2D = new CompositionField2D(ny, nx);

            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    cf_2D.XPositionValues[j, i] = i * dx;
                    cf_2D.YPositionValues[j, i] = j * dy;
                    cf_2D.DValues[j, i] = D[j, i];
                }
            }
            // ====================================
            // Create the solver method matrices
            // ====================================
            nx_less1 = NX - 1;
            // Total points minus the left/right or up/down boundaries
            nx_less2 = NX - 2;
            // nx_less1 - one point because the total entries are one fewer because these are on the off-diagonals
            nx_less3 = nx_less2 - 1;

            end_idx1 = nx_less3;
            end_idx2 = nx_less3 - 1;
            start_idx1 = 0;
            start_idx2 = 1;

            A = new TridiagonalMatrix[nx];
            B_row = new TridiagonalMatrix[nx];
            B_col = new TridiagonalMatrix[nx];
            NeumannBCs_L_A = new double[nx];
            NeumannBCs_R_A = new double[nx];

            double nu0 = dt / (2 * Math.Pow(dx, 2));

            RVector nu = new(nx_less2);
            RVector off_d_val_l0 = new(nx_less3);
            RVector off_d_val_u0 = new(nx_less3);
            RVector main0 = new(nx_less2);
            RVector off_d_val_l1 = new(nx_less3);
            RVector off_d_val_u1 = new(nx_less3);
            RVector main1 = new(nx_less2);
            RVector off_d_val_l2 = new(nx_less3);
            RVector off_d_val_u2 = new(nx_less3);
            RVector main2 = new(nx_less2);

            double D_minus, D_plus;
            for (int i = start_idx1; i < NX; i++)
            {
                RVector D_row = D.GetRowVector(i);
                double D0 = D_row[start_idx1];
                double D1 = D_row[start_idx2];
                double Dend = D_row[nx_less1];
                double Dend1 = D_row[nx_less2];

                // Working this for diffusion heterogeneous media
                //================================================
                // A Matrix first
                //================================================
                #region A
                D_plus = 0.5 * (D1 + D0);
                main0[start_idx1] = 1.0 - 2.0 * D_plus * nu0;
                NeumannBCs_L_A[i] = 1.0 - D_plus * nu0;
                
                D_minus = 0.5 * (Dend1 + Dend);
                main0[end_idx1] = 1.0 - 2.0 * D_minus * nu0;
                NeumannBCs_R_A[i] = 1.0 - D_minus * nu0;

                for (int j = start_idx2; j < end_idx1; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    main0[j] = 1.0 - ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = start_idx2; j < end_idx1+1; j++)
                {
                    D_minus = 0.5 * ( D_row[j] + D_row[j - 1]);
                    off_d_val_l0[j - 1] = nu0 * D_minus;
                }

                for (int j = start_idx1; j < end_idx1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    off_d_val_u0[j] = nu0 * D_plus;
                }

                A[i] = new TridiagonalMatrix(main0, off_d_val_l0, off_d_val_u0);
                #endregion
                ////================================================
                //// The B-row Matrix next
                ////================================================
                #region B_row
                D_plus = 0.5 * (D1 + D0);
                main1[start_idx1] = 1.0 + 2.0 * D_plus * nu0;

                D_minus = 0.5 * (Dend1 + Dend);
                main1[end_idx1] = 1.0 + 2.0 * D_minus * nu0;

                for (int j = start_idx2; j < end_idx1; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    main1[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = start_idx2; j < end_idx1+1; j++)
                {
                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
                    off_d_val_l1[j - 1] = -nu0 * D_minus;
                }

                for (int j = start_idx1; j < end_idx1; j++)
                {
                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
                    off_d_val_u1[j] = -nu0 * D_plus;
                }

                B_row[i] = new TridiagonalMatrix(main1, off_d_val_l1, off_d_val_u1);
                #endregion
                ////================================================
                //// The B-column Matrix last
                ////================================================
                #region B_col
                RVector D_col = D.GetColVector(i);
                D0 = D_col[0];
                D1 = D_col[1];
                Dend = D_col[nx - 1];
                Dend1 = D_col[nx - 2];

                D_plus = 0.5 * (D1 + D0);
                main2[start_idx1] = 1.0 + 2.0 * D_plus * nu0;

                D_minus = 0.5 * (Dend1 + Dend);
                main2[end_idx1] = 1.0 + 2.0 * D_minus * nu0;

                for (int j = start_idx2; j < end_idx1; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
                    D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
                    main2[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
                }

                for (int j = start_idx2; j < end_idx1+1; j++)
                {
                    D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
                    off_d_val_l2[j - 1] = -nu0 * D_minus;
                }

                for (int j = start_idx1; j < end_idx1; j++)
                {
                    D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
                    off_d_val_u2[j] = -nu0 * D_plus;
                }

                B_col[i] = new TridiagonalMatrix(main2, off_d_val_l2, off_d_val_u2);
                #endregion
            }

            // ====================================
            // Set-up the initial condition function and start establishing the initial composition field
            // ===================================
            this.I0 = I0;
            C_Initial = I0(cf_2D.XPositionValues, cf_2D.YPositionValues);            
            // ====================================
            // Set the boundary conditions and apply them to the initial composition field, if needed
            // ====================================
            int num_bounds = Boundary_Conditions.Length;
            border_with_function = new BoundaryWithFunction[num_bounds];
            RVector C0;
            for (int i = 0; i < num_bounds; i++)
            {
                border_with_function[i] = new BoundaryWithFunction
                {
                    BoundaryLocation = i switch
                    {
                        0 => BoundingBox.top,
                        1 => BoundingBox.right,
                        2 => BoundingBox.left,
                        3 => BoundingBox.bottom,
                        _ => BoundingBox.bottom,
                    },
                    TypeBC = Boundary_Conditions[i],
                    BoundaryFunction = bc_s[i]
                };
                switch (border_with_function[i].BoundaryLocation)
                {
                    case BoundingBox.top:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[X.GetnRows - 1, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(ny - 2) + C0;
                            C_Initial.ReplaceRow(Ctab, ny - 1);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {   
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[X.GetnRows - 1, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(ny - 2) + C0;
                            C_Initial.ReplaceRow(Ctab, ny - 1);
                        }
                        break;
                    case BoundingBox.right:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 2) + C0;
                            C_Initial.ReplaceCol(Ctab, nx - 1);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, Y.GetnCols - 1];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(nx - 2) + C0;
                            C_Initial.ReplaceCol(Ctab, nx - 1);
                        }
                        break;
                    case BoundingBox.left:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(1) + C0;
                            C_Initial.ReplaceCol(Ctab, 0);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = Y.GetColVector(0);
                            border_with_function[i].PositionFixed = X[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetColVector(1) + C0;
                            C_Initial.ReplaceCol(C0, 0);
                        }
                        break;
                    case BoundingBox.bottom:
                        if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(1) + C0;
                            C_Initial.ReplaceRow(Ctab, 0);
                        }
                        else if (border_with_function[i].TypeBC == ABoundaryCondition.neumann)
                        {
                            border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
                            border_with_function[i].PositionFixed = Y[0, 0];
                            C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

                            RVector Ctab = C_Initial.GetRowVector(1) + C0;
                            C_Initial.ReplaceRow(Ctab, 0);
                        }
                        break;
                }
            }

            gxt_function = g;
            Chat_mode = Tmode;
            b_filename = base_filename;
            Errors = new string[nt];
        }

        // ====================================================================
        // Solvers
        // ====================================================================
        /// <summary>
        /// Method for solving the 2D diffusion equation using the 2D Alternating Direction Implicit algorithm.  Will output time-steps elapsed.
        /// </summary>
        /// <param name="n_time_steps"></param>
        /// <param name="output_interval"></param>
        public void Solve(int n_time_steps, int output_interval)
        {
            int nrows = cf_2D.InitialCompositionValues.GetnRows;
            int ncols = cf_2D.InitialCompositionValues.GetnCols;

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
            RVector CT0 = new(ncols);
            RVector CT1 = new(ncols);
            RVector CB0 = new(ncols);
            RVector CB1 = new(ncols);

            double nu0 = dt / (2 * Math.Pow(dx, 2));
            string full_file_name;
            // Time evolution           
            for (int t = 0; t < n_time_steps; t++)
            {
                try
                {
                    if (Chat_mode == Mode.verbose && output_interval > 0 && Base_filename != null)
                    {
                        if (t % output_interval == 0)
                        {
                            decimal time = (decimal)(t * dt);
                            Console.ForegroundColor = ConsoleColor.Green;
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
                    CR = BCs_Functions[1].BoundaryFunction(t * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
                    CL = BCs_Functions[2].BoundaryFunction(t * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
                    double holdA0 = 0.0, holdA1 = 0.0;
                    for (int i = start_idx2; i < nx_less1; i++)
                    {
                        RVector xold;
                        if (t == 0) { xold = C_Initial.GetRowVector(i, start_idx2, nx_less1); }
                        else { xold = C_Im2.GetRowVector(i, start_idx2, nx_less1); }
                        if (BCs_Functions[2].TypeBC == ABoundaryCondition.neumann) { holdA0 = A[i][start_idx1, start_idx1];  A[i][start_idx1, start_idx1] = NeumannBCs_L_A[i]; }
                        if (BCs_Functions[1].TypeBC == ABoundaryCondition.neumann) { holdA1 = A[i][end_idx1, end_idx1];  A[i][end_idx1, end_idx1] = NeumannBCs_R_A[i]; }
                       
                        RVector v1 = A[i].Dot(xold);
                        C_Ex.ReplaceRow(v1, i, start_idx2, end_idx1);
                                                                        
                        switch (BCs_Functions[2].TypeBC)
                        {
                            case ABoundaryCondition.dirichlet:
                                C_Ex[i, start_idx1] = CL[i];
                                C_Ex[i, start_idx2] = C_Ex[i, start_idx2] + (-(NeumannBCs_L_A[i]-1)) * CL[i];
                                break;
                            case ABoundaryCondition.neumann:                              
                                C_Ex[i, start_idx2] = C_Ex[i, start_idx2] + (-(NeumannBCs_L_A[i] - 1)) * dx * CL[i];
                                C_Ex[i, start_idx1] = C_Ex[i, start_idx2];
                                break;
                        }
                        switch (BCs_Functions[1].TypeBC)
                        {
                            case ABoundaryCondition.dirichlet:
                                C_Ex[i, nx_less1] = CR[i];
                                C_Ex[i, nx_less2] = C_Ex[i, nx_less2] + (-(NeumannBCs_R_A[i] - 1)) * CR[i];
                                break;
                            case ABoundaryCondition.neumann:
                                C_Ex[i, nx_less2] = C_Ex[i, nx_less2] + (-(NeumannBCs_R_A[i] - 1)) * dx * CR[i];
                                C_Ex[i, nx_less1] = C_Ex[i, nx_less2];
                                break;
                        }

                        if (BCs_Functions[2].TypeBC == ABoundaryCondition.neumann) { A[i][start_idx1, start_idx1] = holdA0; }
                        if (BCs_Functions[1].TypeBC == ABoundaryCondition.neumann) { A[i][end_idx1, end_idx1] = holdA1; }
                    }
                    // ===================
                    // ===================
                    // One-half implicit time-step
                    // ===================
                    // Source terms
                    double t12 = (t + 0.5) * dt;
                    if (t == 0) { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Initial); }
                    else { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Im2); }

                    f12 = dt / 2.0 * fn;

                    // BCs
                    RVector gn, g0;
                    gn = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                    g0 = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                    CT = (0.5 * B_row[nx_less1].Dot(gn.Section(start_idx2, nx_less1))) + (0.5 * A[nx_less1].Dot(g0.Section(start_idx2, nx_less1)));
                    
                    gn = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                    g0 = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                    CB = (0.5 * B_row[start_idx1].Dot(gn.Section(start_idx2, nx_less1))) + (0.5 * A[start_idx1].Dot(g0.Section(start_idx2, nx_less1)));

                    TridiagonalMatrix B_col2;
                    for (int j = start_idx2; j < nx_less1; j++)
                    {
                        RVector v1 = C_Ex.GetColVector(j);
                        v1[0] = CB[j]; //nu * 
                        v1[ncols - 1] = CT[j]; //nu * 

                        RVector f12s = f12.GetColVector(j);
                        RVector v2, u12;
                        double s = -B_col[j][1, 2];

                        // Bottom dirichlet, Top dirichlet
                        if (BCs_Functions[3].TypeBC == ABoundaryCondition.dirichlet && BCs_Functions[0].TypeBC == ABoundaryCondition.dirichlet) 
                        {
                            v2 = (v1 + f12s).Section(start_idx2, nx_less1);
                            v2[0] = s * (v1[0] + f12s[0]);
                            v2[v2.GetRVectorSize - 1] = s * (v1[nx_less1] + f12s[nx_less1]);
                            B_col2 = B_col[j];
                            
                            u12 = TridiagonalMatrix.Thomas_Algorithm(B_col2, v2);
                            
                            int ns = (nrows - 1) - u12.GetRVectorSize;
                            C_Im1.ReplaceCol(u12, j, ns, u12.GetRVectorSize + 1);
                            C_Im1[0, j] = CB[j];
                            C_Im1[nrows - 1, j] = CT[j];
                        }
                        // Bottom neumann, Top dirichlet
                        else if (BCs_Functions[3].TypeBC == ABoundaryCondition.neumann && BCs_Functions[0].TypeBC == ABoundaryCondition.dirichlet) 
                        { 
                            v2 = (v1 + f12s).Section(start_idx1, nx_less1);
                            v2[0] = 2 * s * dx * (v1[0] + f12s[0]);
                            v2[v2.GetRVectorSize - 1] = s * (v1[nx_less1] + f12s[nx_less1]);
                            
                            B_col2 = new(v2.GetRVectorSize, v2.GetRVectorSize);
                            B_col2.FitMainUpperLower(start_idx2, nx_less1, start_idx2, nx_less1-1, start_idx2, nx_less1-1, B_col[j].GetMain(), B_col[j].GetUpper(), B_col[j].GetLower());
                            B_col2[0, 1] = -2 * s;
                            B_col2[0, 0] = B_col2[1, 1];
                            B_col2[1, 0] = B_col2[1, 2];

                            u12 = TridiagonalMatrix.Thomas_Algorithm(B_col2, v2);
                            
                            int ns = start_idx1;
                            C_Im1.ReplaceCol(u12, j, ns, u12.GetRVectorSize);
                            C_Im1[nrows - 1, j] = CT[j];
                        }
                        // Bottom dirichlet, Top neumann
                        else if (BCs_Functions[3].TypeBC == ABoundaryCondition.dirichlet && BCs_Functions[0].TypeBC == ABoundaryCondition.neumann) 
                        { 
                            v2 = (v1 + f12s).Section(start_idx2, NX);
                            v2[0] = s * (v1[0] + f12s[0]);
                            v2[v2.GetRVectorSize - 1] = 2 * s * dx * (v1[nx_less1] + f12s[nx_less1]);

                            B_col2 = new(v2.GetRVectorSize, v2.GetRVectorSize);
                            B_col2.FitMainUpperLower(start_idx1, nx_less2, start_idx1, nx_less2 - 1, start_idx1, nx_less2 - 1, B_col[j].GetMain(), B_col[j].GetUpper(), B_col[j].GetLower());
                            B_col2[nx_less2, nx_less3] = -2 * s;
                            B_col2[nx_less2, nx_less2] = B_col2[nx_less3, nx_less3];
                            B_col2[nx_less3, nx_less2] = B_col2[nx_less3, nx_less3 - 1];

                            u12 = TridiagonalMatrix.Thomas_Algorithm(B_col2, v2);

                            int ns = start_idx2;
                            C_Im1.ReplaceCol(u12, j, ns, u12.GetRVectorSize);
                            C_Im1[0, j] = CB[j];
                        }
                        // Bottom neumann, top neumann
                        else 
                        { 
                            v2 = v1 + f12s;
                            v2[0] = 2 * s * dx * (v1[0] + f12s[0]);
                            v2[v2.GetRVectorSize - 1] = 2 * s * dx * (v1[nx_less1] + f12s[nx_less1]);

                            B_col2 = new(v2.GetRVectorSize, v2.GetRVectorSize);
                            B_col2.FitMainUpperLower(start_idx2, nx_less1, start_idx2, nx_less1 - 1, start_idx2, nx_less1 - 1, B_col[j].GetMain(), B_col[j].GetUpper(), B_col[j].GetLower());
                            
                            B_col2[0, 1] = -2 * s;
                            B_col2[0, 0] = B_col2[1, 1];
                            B_col2[1, 0] = B_col2[1, 2];

                            B_col2[nx_less1, nx_less2] = -2 * s;
                            B_col2[nx_less1, nx_less1] = B_col2[nx_less3, nx_less3];
                            B_col2[nx_less2, nx_less1] = B_col2[nx_less3, nx_less3 - 1];

                            u12 = TridiagonalMatrix.Thomas_Algorithm(B_col2, v2);

                            int ns = start_idx1;
                            C_Im1.ReplaceCol(u12, j, ns, u12.GetRVectorSize);
                        }                                               
                    }
                    // ===================
                    // ===================
                    // Full implicit time-step
                    // ===================
                    switch (BCs_Functions[1].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            CR = BCs_Functions[1].BoundaryFunction((t + 1) * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
                            //CT = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                            //RVector ctab = C_Initial.GetRowVector(0);
                            //C_Im2.ReplaceRow(CT, nrows - 1);
                            break;
                        case ABoundaryCondition.neumann:
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
                            CR = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(ncols - 1), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(ncols - 2), C2); // + C3
                            //CT = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                            //C_Im2.ReplaceRow(CT, nrows - 1);
                            break;
                        default:
                            break;
                    }
                    switch (BCs_Functions[2].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            CL = BCs_Functions[2].BoundaryFunction((t + 1) * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
                            //CB = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);

                            //RVector ctab = C_Initial.GetRowVector(nrows - 1);

                            //C_Im2.ReplaceRow(CB, 0);
                            break;
                        case ABoundaryCondition.neumann:
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
                            CL = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(0), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(1), C2); // + C3
                            //CB = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                            //C_Im2.ReplaceRow(CB, 0);
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

                        switch (BCs_Functions[0].TypeBC)
                        {
                            case ABoundaryCondition.dirichlet:
                                B_row[k][0, 0] = 1.0;
                                B_row[k][0, 1] = 0.0;
                                break;
                            case ABoundaryCondition.neumann:
                                B_row[k][0, 0] = nu0;
                                //B_row[k][0, 1] = 0.0;
                                break;
                            default:
                                B_row[k][0, 0] = 1.0;
                                B_row[k][0, 1] = 0.0;
                                break;
                        }
                        switch (BCs_Functions[3].TypeBC)
                        {
                            case ABoundaryCondition.dirichlet:
                                B_row[k][ncols - 1, ncols - 1] = 1.0;
                                B_row[k][ncols - 1, ncols - 2] = 0.0;
                                break;
                            case ABoundaryCondition.neumann:
                                B_row[k][ncols - 1, ncols - 1] = nu0;
                                //B_row[k][0, 1] = 0.0;
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
                    // Set boundary conditions on the C_Im2 matrix
                    // ===================
                    switch (BCs_Functions[0].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            CT = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
                            C_Im2.ReplaceRow(CT, nrows - 1);
                            break;
                        case ABoundaryCondition.neumann:
                            break;
                    }
                    switch (BCs_Functions[3].TypeBC)
                    {
                        case ABoundaryCondition.dirichlet:
                            CB = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
                            C_Im2.ReplaceRow(CB, 0);
                            break;
                        case ABoundaryCondition.neumann:
                            break;
                    }
                    // ===================
                }
                catch (Exception e)
                {
                    Errors[t] = e.Message;
                    error_flag = true;
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine("Fault!");
                }
            }
            // Setup the final composition field
            if (error_flag == false) { C_Final = C_Im2; }
            else { C_Final = C_Initial; }
        }
        //public void Solve(int n_time_steps, int output_interval)
        //{            
        //    int nrows = cf_2D.InitialCompositionValues.GetnRows;
        //    int ncols = cf_2D.InitialCompositionValues.GetnCols;

        //    RMatrix C_Ex = new(nrows, ncols);
        //    RMatrix C_Im1 = new(nrows, ncols);
        //    RMatrix C_Im2 = new(nrows, ncols);

        //    RMatrix fn = new(nrows, ncols);
        //    RMatrix f0 = new(nrows, ncols);
        //    RMatrix f12 = new(nrows, ncols);

        //    RVector CL = new(nrows);
        //    RVector CR = new(nrows);
        //    RVector CT = new(ncols);
        //    RVector CB = new(ncols);
        //    RVector CT0 = new(ncols);
        //    RVector CT1 = new(ncols);
        //    RVector CB0 = new(ncols);
        //    RVector CB1 = new(ncols);

        //    double nu0 = dt / (2 * Math.Pow(dx, 2));
        //    string full_file_name;
        //    // Time evolution           
        //    for (int t = 0; t < n_time_steps; t++)
        //    {
        //        try
        //        {
        //            if (Chat_mode == Mode.verbose && output_interval > 0 && Base_filename != null)
        //            {
        //                if (t % output_interval == 0)
        //                {
        //                    decimal time = (decimal)(t * dt);
        //                    Console.ForegroundColor = ConsoleColor.Green;
        //                    Console.WriteLine("{0}s have been simulated", time);
        //                    full_file_name = Base_filename + time.ToString() + suffix;
        //                    if (File.Exists(full_file_name)) { File.Delete(full_file_name); }
        //                    if (t == 0) { FileWriteData_CSV(full_file_name, X, Y, C_Initial); }
        //                    else { FileWriteData_CSV(full_file_name, X, Y, C_Im2); }
        //                }
        //            }

        //            // ===================
        //            // Explicit time-step
        //            // ===================
        //            // 0 = top, 1 = right, 2 = left, 3 = bottom
        //            switch (BCs_Functions[1].TypeBC)
        //            {
        //                case ABoundaryCondition.dirichlet:
        //                    CR = BCs_Functions[1].BoundaryFunction(t * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
        //                    break;
        //                case ABoundaryCondition.neumann:
        //                    RVector C1, C2; //,C3;
        //                    double cUpper1, cUpper2, cLower1, cLower2;
        //                    if (t == 0)
        //                    {
        //                        C1 = C_Initial.GetColVector(ncols - 1);
        //                        C2 = C_Initial.GetColVector(ncols - 2);

        //                        cLower1 = C_Initial[0, ncols - 1];
        //                        cLower2 = C_Initial[1, ncols - 2];
        //                        cUpper1 = C_Initial[nrows - 1, ncols - 1];
        //                        cUpper2 = C_Initial[nrows - 2, ncols - 2];
        //                    }
        //                    else
        //                    {
        //                        C1 = C_Im2.GetColVector(ncols - 1);
        //                        C2 = C_Im2.GetColVector(ncols - 2);

        //                        cLower1 = C_Im2[0, ncols - 1];
        //                        cLower2 = C_Im2[1, ncols - 2];
        //                        cUpper1 = C_Im2[nrows - 1, ncols - 1];
        //                        cUpper2 = C_Im2[nrows - 2, ncols - 2];
        //                    }
        //                    CR = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(ncols - 1), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(ncols - 2), C2); // + C3                            
        //                    CR[0] = (-2 * nu0 * cf_2D.DValues[0, ncols - 1] * cLower1) + (2 * nu0 * cf_2D.DValues[1, ncols - 2] * cLower2);
        //                    CR[ncols - 1] = (-2 * nu0 * cf_2D.DValues[nrows - 1, ncols - 1] * cUpper1) + (2 * nu0 * cf_2D.DValues[nrows - 2, ncols - 2] * cUpper2);
        //                    break;
        //                default:
        //                    break;
        //            }
        //            switch (BCs_Functions[2].TypeBC)
        //            {
        //                case ABoundaryCondition.dirichlet:
        //                    CL = BCs_Functions[2].BoundaryFunction(t * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
        //                    break;
        //                case ABoundaryCondition.neumann:
        //                    RVector C1, C2; //,C3;
        //                    double cUpper1, cUpper2, cLower1, cLower2;
        //                    if (t == 0)
        //                    {
        //                        C1 = C_Initial.GetColVector(0);
        //                        C2 = C_Initial.GetColVector(1);

        //                        cLower1 = C_Initial[0, 0];
        //                        cLower2 = C_Initial[1, 1];
        //                        cUpper1 = C_Initial[nrows - 1, 0];
        //                        cUpper2 = C_Initial[nrows - 2, 1];
        //                    }
        //                    else
        //                    {
        //                        C1 = C_Im2.GetColVector(0);
        //                        C2 = C_Im2.GetColVector(1);

        //                        cLower1 = C_Im2[0, 0];
        //                        cLower2 = C_Im2[1, 1];
        //                        cUpper1 = C_Im2[nrows - 1, 0];
        //                        cUpper2 = C_Im2[nrows - 2, 1];
        //                    }
        //                    CL = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(0), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(1), C2); // + C3
        //                    CL[0] = (-2 * nu0 * cf_2D.DValues[0, 0] * cLower1) + (2 * nu0 * cf_2D.DValues[1, 1] * cLower2);
        //                    CL[ncols - 1] = (-2 * nu0 * cf_2D.DValues[nrows - 1, 0] * cUpper1) + (2 * nu0 * cf_2D.DValues[nrows - 2, 1] * cUpper2);
        //                    break;

        //                default:
        //                    break;
        //            }

        //            for (int i = 1; i < nrows - 1; i++)
        //            {
        //                RVector xold;
        //                if (t == 0) { xold = C_Initial.GetRowVector(i); }
        //                else { xold = C_Im2.GetRowVector(i); } // C_Im2.GetRowVector(i); }
                        
        //                switch (BCs_Functions[1].TypeBC)
        //                {
        //                    case ABoundaryCondition.dirichlet:
        //                        A[i][nrows - 1, ncols - 1] = 1.0;
        //                        A[i][nrows - 1, ncols - 2] = 0.0;
        //                        break;
        //                    case ABoundaryCondition.neumann:
        //                        break;
        //                }
        //                switch (BCs_Functions[1].TypeBC)
        //                {
        //                    case ABoundaryCondition.dirichlet:
        //                        A[i][0, 0] = 1.0;
        //                        A[i][0, 1] = 0.0;
        //                        break;
        //                    case ABoundaryCondition.neumann:
        //                        break;
        //                }

        //                RVector v1 = A[i].Dot(xold);
        //                C_Ex.ReplaceRow(v1, i);
        //                C_Ex[i, 0] = CL[i]; //nu *
        //                C_Ex[i, ncols - 1] = CR[i]; //nu *                        
        //            }
        //            // ===================
        //            // ===================
        //            // One-half implicit time-step
        //            // ===================
        //            // Source terms
        //            double t12 = (t + 0.5) * dt;
        //            if (t == 0) { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Initial); }
        //            else { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Im2); }

        //            f12 = dt / 2.0 * fn;

        //            // BCs
        //            switch (BCs_Functions[0].TypeBC)
        //            {
        //                case ABoundaryCondition.dirichlet:
        //                    RVector gn = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
        //                    RVector g0 = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
        //                    CT = (0.5 * B_row[nrows - 1].Dot(gn)) + (0.5 * A[nrows - 1].Dot(g0));  
        //                    break;
        //                case ABoundaryCondition.neumann:
        //                    RVector C1, C2; //,C3;
        //                    double cUpper1, cUpper2, cLower1, cLower2;
        //                    if (t == 0)
        //                    {
        //                        C1 = C_Initial.GetRowVector(nrows - 1);
        //                        C2 = C_Initial.GetRowVector(nrows - 2);

        //                        cLower1 = C_Initial[0, ncols - 1];
        //                        cLower2 = C_Initial[1, ncols - 2];
        //                        cUpper1 = C_Initial[nrows - 1, ncols - 1];
        //                        cUpper2 = C_Initial[nrows - 2, ncols - 2];
        //                    }
        //                    else
        //                    {
        //                        C1 = C_Im2.GetRowVector(nrows - 1);
        //                        C2 = C_Im2.GetRowVector(nrows - 2);
                                
        //                        cLower1 = C_Im2[0, ncols - 1];
        //                        cLower2 = C_Im2[1, ncols - 2];
        //                        cUpper1 = C_Im2[nrows - 1, ncols - 1];
        //                        cUpper2 = C_Im2[nrows - 2, ncols - 2];
        //                    }
        //                    CT0 = RVector.Product(-2 * nu0 * cf_2D.DValues.GetRowVector(nrows - 1), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetRowVector(nrows - 2), C2); // + C3
        //                    CT1 = A[nrows - 1].Dot(RVector.Product(-2 * nu0 * cf_2D.DValues.GetRowVector(nrows - 1), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetRowVector(nrows - 2), C2)); // + C3
        //                    CT = (0.5 * CT0) + (0.5 * CT1);
        //                    break;
        //                default:
        //                    break;
        //            }
        //            switch (BCs_Functions[3].TypeBC)
        //            {
        //                case ABoundaryCondition.dirichlet:
        //                    RVector gn = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
        //                    RVector g0 = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
        //                    CB = (0.5 * B_row[0].Dot(gn)) + (0.5 * A[0].Dot(g0)); //((0.5 * Bm.Dot(gn)) + (0.5 * Bp.Dot(g0))); //nu * 
        //                    break;
        //                case ABoundaryCondition.neumann:
        //                    RVector C1, C2; //,C3;
        //                    if (t == 0)
        //                    {
        //                        C1 = C_Initial.GetRowVector(0);
        //                        C2 = C_Initial.GetRowVector(1);
        //                    }
        //                    else
        //                    {
        //                        C1 = C_Im2.GetRowVector(0);
        //                        C2 = C_Im2.GetRowVector(1);
        //                    }
        //                    CB0 = RVector.Product(-2 * nu0 * cf_2D.DValues.GetRowVector(0), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetRowVector(1), C2); // + C3
        //                    CB1 = A[nrows - 1].Dot(RVector.Product(-2 * nu0 * cf_2D.DValues.GetRowVector(0), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetRowVector(1), C2)); // + C3
        //                    CB = (0.5 * CB0) + (0.5 * CB1);
        //                    //Console.WriteLine(CB.ToString());
        //                    //Console.ReadKey();
        //                    break;
        //                default:
        //                    break;
        //            }
        //            for (int j = 1; j < ncols - 1; j++)
        //            {
        //                RVector v1 = C_Ex.GetColVector(j);
        //                v1[0] = CB[j]; //nu * 
        //                v1[ncols - 1] = CT[j]; //nu * 

        //                RVector f12s = f12.GetColVector(j); // 

        //                switch (BCs_Functions[0].TypeBC)
        //                {
        //                    case ABoundaryCondition.dirichlet:
        //                        B_col[j][0, 0] = 1.0;
        //                        B_col[j][0, 1] = 0.0;
        //                        break;
        //                    case ABoundaryCondition.neumann:
        //                        B_col[j][0, 0] = nu0;
        //                        //B_col[j][0, 1] = -nu0;
        //                        break;
        //                    default:
        //                        B_col[j][0, 0] = 1.0;
        //                        B_col[j][0, 1] = 0.0;
        //                        break;
        //                }
        //                switch (BCs_Functions[3].TypeBC)
        //                {
        //                    case ABoundaryCondition.dirichlet:
        //                        B_col[j][ncols - 1, ncols - 1] = 1.0;
        //                        B_col[j][ncols - 1, ncols - 2] = 0.0;
        //                        break;
        //                    case ABoundaryCondition.neumann:
        //                        B_col[j][ncols - 1, ncols - 1] = nu0;
        //                        //B_col[j][0, 1] = -nu0;
        //                        break;
        //                    default:
        //                        B_col[j][ncols - 1, ncols - 1] = 1.0;
        //                        B_col[j][ncols - 1, ncols - 2] = 0.0;
        //                        break;
        //                }
        //                RVector u12 = TridiagonalMatrix.Thomas_Algorithm(B_col[j], v1 + f12s);
        //                C_Im1.ReplaceCol(u12, j);                        
        //            }

        //            // ===================
        //            // ===================
        //            // Full implicit time-step
        //            // ===================
        //            switch (BCs_Functions[1].TypeBC)
        //            {
        //                case ABoundaryCondition.dirichlet:
        //                    CR = BCs_Functions[1].BoundaryFunction((t + 1) * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
        //                    //CT = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
        //                    //RVector ctab = C_Initial.GetRowVector(0);
        //                    //C_Im2.ReplaceRow(CT, nrows - 1);
        //                    break;
        //                case ABoundaryCondition.neumann:
        //                    RVector C1, C2; //,C3;
        //                    if (t == 0)
        //                    {
        //                        C1 = C_Initial.GetColVector(ncols - 1);
        //                        C2 = C_Initial.GetColVector(ncols - 2);
        //                    }
        //                    else
        //                    {
        //                        C1 = C_Im2.GetColVector(ncols - 1);
        //                        C2 = C_Im2.GetColVector(ncols - 2);
        //                    }
        //                    //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
        //                    CR = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(ncols - 1), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(ncols - 2), C2); // + C3
        //                    //CT = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
        //                    //C_Im2.ReplaceRow(CT, nrows - 1);
        //                    break;
        //                default:
        //                    break;
        //            }
        //            switch (BCs_Functions[2].TypeBC)
        //            {
        //                case ABoundaryCondition.dirichlet:
        //                    CL = BCs_Functions[2].BoundaryFunction((t + 1) * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
        //                    //CB = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);

        //                    //RVector ctab = C_Initial.GetRowVector(nrows - 1);

        //                    //C_Im2.ReplaceRow(CB, 0);
        //                    break;
        //                case ABoundaryCondition.neumann:
        //                    RVector C1, C2; //,C3;
        //                    if (t == 0)
        //                    {
        //                        C1 = C_Initial.GetColVector(0);
        //                        C2 = C_Initial.GetColVector(1);
        //                    }
        //                    else
        //                    {
        //                        C1 = C_Im2.GetColVector(0);
        //                        C2 = C_Im2.GetColVector(1);
        //                    }
        //                    //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
        //                    CL = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(0), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(1), C2); // + C3
        //                    //CB = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
        //                    //C_Im2.ReplaceRow(CB, 0);
        //                    break;
        //                default:
        //                    break;
        //            }
        //            for (int k = 1; k < nrows - 1; k++)
        //            {
        //                RVector v1 = C_Ex.GetRowVector(k);
        //                RVector u12 = C_Im1.GetRowVector(k);
        //                RVector b = (2 * u12) - v1;
        //                b[0] = CL[k]; //nu * 
        //                b[ncols - 1] = CR[k]; //nu * 

        //                switch (BCs_Functions[0].TypeBC)
        //                {
        //                    case ABoundaryCondition.dirichlet:
        //                        B_row[k][0, 0] = 1.0;
        //                        B_row[k][0, 1] = 0.0;
        //                        break;
        //                    case ABoundaryCondition.neumann:
        //                        B_row[k][0, 0] = nu0;
        //                        //B_row[k][0, 1] = 0.0;
        //                        break;
        //                    default:
        //                        B_row[k][0, 0] = 1.0;
        //                        B_row[k][0, 1] = 0.0;
        //                        break;
        //                }
        //                switch (BCs_Functions[3].TypeBC)
        //                {
        //                    case ABoundaryCondition.dirichlet:
        //                        B_row[k][ncols - 1, ncols - 1] = 1.0;
        //                        B_row[k][ncols - 1, ncols - 2] = 0.0;
        //                        break;
        //                    case ABoundaryCondition.neumann:
        //                        B_row[k][ncols - 1, ncols - 1] = nu0;
        //                        //B_row[k][0, 1] = 0.0;
        //                        break;
        //                    default:
        //                        B_row[k][ncols - 1, ncols - 1] = 1.0;
        //                        B_row[k][ncols - 1, ncols - 2] = 0.0;
        //                        break;
        //                }
        //                RVector u1 = TridiagonalMatrix.Thomas_Algorithm(B_row[k], b);
        //                u1[0] = CL[k];  //nu * 
        //                u1[ncols - 1] = CR[k]; //nu * 
        //                C_Im2.ReplaceRow(u1, k);
        //            }
        //            // ===================
        //            // Set boundary conditions on the C_Im2 matrix
        //            // ===================
        //            switch (BCs_Functions[0].TypeBC)
        //            {
        //                case ABoundaryCondition.dirichlet:
        //                    CT = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
        //                    C_Im2.ReplaceRow(CT, nrows - 1);
        //                    break;
        //                case ABoundaryCondition.neumann:
        //                    break;
        //            }
        //            switch (BCs_Functions[3].TypeBC)
        //            {
        //                case ABoundaryCondition.dirichlet:
        //                    CB = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
        //                    C_Im2.ReplaceRow(CB, 0);
        //                    break;
        //                case ABoundaryCondition.neumann:
        //                    break;
        //            }
        //            // ===================
        //        }
        //        catch (Exception e)
        //        {
        //            Errors[t] = e.Message;
        //            error_flag = true;
        //            Console.ForegroundColor = ConsoleColor.Red;
        //            Console.WriteLine("Fault!");
        //        }
        //    }
        //    // Setup the final composition field
        //    if (error_flag == false) { C_Final = C_Im2; }
        //    else { C_Final = C_Initial; }            
        //}




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

        //        public DiffusionSimulators_2D_MatrixD(RMatrix D, double dx, double dy, int nx, int ny, double dt, int nt,
        //ABoundaryCondition[] Boundary_Conditions, RVector[] bc_s, InitialCondition_Del I0, SourceTerm_MatrixD_Del g, Mode Tmode, string base_filename)
        //        {
        //            this.dx = dx;
        //            this.dy = dy;
        //            this.dt = dt;
        //            this.diffusivity = D;
        //            CF_2D = new CompositionField2D(ny, nx);

        //            for (int i = 0; i < ny; i++)
        //            {
        //                for (int j = 0; j < nx; j++)
        //                {
        //                    CF_2D.xposition_matrix[j, i] = i * dx;
        //                    CF_2D.yposition_matrix[j, i] = j * dy;
        //                    CF_2D.DiffusionCoefficient[j, i] = D[j, i];
        //                }
        //            }

        //            //C_Initial = I0;
        //            C_Initial = I0(CF_2D.xposition_matrix, CF_2D.yposition_matrix);
        //            this.I0 = I0;

        //            A = new TridiagonalMatrix[nx];
        //            B_row = new TridiagonalMatrix[nx];
        //            B_col = new TridiagonalMatrix[nx];

        //            double nu0 = dt / (2 * Math.Pow(dx, 2));

        //            RVector nu = new(nx);
        //            RVector off_d_val_l0 = new(nx - 1);
        //            RVector off_d_val_u0 = new(nx - 1);
        //            RVector main0 = new(nx);
        //            RVector off_d_val_l1 = new(nx - 1);
        //            RVector off_d_val_u1 = new(nx - 1);
        //            RVector main1 = new(nx);
        //            RVector off_d_val_l2 = new(nx - 1);
        //            RVector off_d_val_u2 = new(nx - 1);
        //            RVector main2 = new(nx);

        //            double D_minus, D_plus;
        //            for (int i = 0; i < nx; i++)
        //            {
        //                RVector D_row = D.GetRowVector(i);
        //                //nu = nu0 * D_row;
        //                //main = 1 - (2 * nu);
        //                // Working this for diffusion heterogeneous media
        //                //================================================
        //                // A Matrix first
        //                //================================================
        //                D_plus = 0.5 * (D_row[1] + D_row[0]);
        //                main0[0] = 1.0 - ((D_plus * nu0));

        //                D_minus = 0.5 * (D_row[nx - 2] + D_row[0]);
        //                main0[nx - 1] = 1.0 - ((D_minus * nu0));

        //                for (int j = 1; j < nx - 1; j++)
        //                {
        //                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
        //                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
        //                    main0[j] = 1.0 - ((D_minus * nu0) + (D_plus * nu0));
        //                }

        //                for (int j = 1; j < nx; j++)
        //                {
        //                    D_minus = 0.5 * (D_row[j] + D_row[j - 1]);
        //                    off_d_val_l0[j - 1] = nu0 * D_minus;
        //                }

        //                for (int j = 0; j < nx - 1; j++)
        //                {
        //                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
        //                    off_d_val_u0[j] = nu0 * D_plus;
        //                }

        //                A[i] = new TridiagonalMatrix(main0, off_d_val_l0, off_d_val_u0);
        //                ////================================================
        //                //// The B-row Matrix next
        //                ////================================================
        //                //main = 1 + (2 * nu);

        //                D_plus = 0.5 * (D_row[1] + D_row[0]);
        //                main1[0] = 1.0 + ((D_plus * nu0));

        //                D_minus = 0.5 * (D_row[nx - 2] + D_row[0]);
        //                main1[nx - 1] = 1.0 + ((D_minus * nu0));

        //                for (int j = 1; j < nx - 1; j++)
        //                {
        //                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
        //                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
        //                    main1[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
        //                }

        //                for (int j = 1; j < nx; j++)
        //                {
        //                    D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
        //                    off_d_val_l1[j - 1] = -nu0 * D_minus;
        //                }

        //                for (int j = 0; j < nx - 1; j++)
        //                {
        //                    D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
        //                    off_d_val_u1[j] = -nu0 * D_plus;
        //                }

        //                B_row[i] = new TridiagonalMatrix(main1, off_d_val_l1, off_d_val_u1);
        //                ////================================================
        //                //// The B-column Matrix last
        //                ////================================================
        //                RVector D_col = D.GetColVector(i);

        //                D_plus = 0.5 * (D_col[1] + D_col[0]);
        //                main2[0] = 1.0 + ((D_plus * nu0));

        //                D_minus = 0.5 * (D_col[nx - 2] + D_col[0]);
        //                main2[nx - 1] = 1.0 + ((D_minus * nu0));

        //                for (int j = 1; j < nx - 1; j++)
        //                {
        //                    D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
        //                    D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
        //                    main2[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
        //                }

        //                for (int j = 1; j < nx; j++)
        //                {
        //                    D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
        //                    off_d_val_l2[j - 1] = -nu0 * D_minus;
        //                }

        //                for (int j = 0; j < nx - 1; j++)
        //                {
        //                    D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
        //                    off_d_val_u2[j] = -nu0 * D_plus;
        //                }

        //                ////for (int j = 0; j < nx - 1; j++) { off_d_val_l[j] = nu[j]; }
        //                ////for (int j = 1; j < nx; j++) { off_d_val_u[j - 1] = nu[j]; }

        //                B_col[i] = new TridiagonalMatrix(main2, off_d_val_l2, off_d_val_u2);
        //            }

        //            int num_bounds = Boundary_Conditions.Length;
        //            border_with_vector = new BoundaryWithVector[num_bounds];
        //            RVector C0;
        //            for (int i = 0; i < num_bounds; i++)
        //            {
        //                border_with_vector[i] = new BoundaryWithVector
        //                {
        //                    BoundaryLocation = i switch
        //                    {
        //                        0 => BoundingBox.top,
        //                        1 => BoundingBox.right,
        //                        2 => BoundingBox.left,
        //                        3 => BoundingBox.bottom,
        //                        _ => BoundingBox.bottom,
        //                    },
        //                    TypeBC = Boundary_Conditions[i],
        //                    FunctionValues = bc_s[i]
        //                };
        //                switch (border_with_vector[i].BoundaryLocation)
        //                {
        //                    case BoundingBox.top:
        //                        if (border_with_vector[i].TypeBC == ABoundaryCondition.dirichlet)
        //                        {
        //                            C0 = border_with_vector[i].FunctionValues;
        //                            RVector Ctab = C_Initial.GetRowVector(ny - 1) + C0;
        //                            C_Initial.ReplaceRow(C0, ny - 1);
        //                        }
        //                        break;
        //                    case BoundingBox.right:
        //                        if (border_with_vector[i].TypeBC == ABoundaryCondition.dirichlet)
        //                        {
        //                            C0 = border_with_vector[i].FunctionValues;
        //                            RVector Ctab = C_Initial.GetColVector(nx - 1) + C0;
        //                            C_Initial.ReplaceCol(C0, nx - 1);
        //                        }
        //                        break;
        //                    case BoundingBox.left:
        //                        if (border_with_vector[i].TypeBC == ABoundaryCondition.dirichlet)
        //                        {
        //                            C0 = border_with_vector[i].FunctionValues;
        //                            RVector Ctab = C_Initial.GetColVector(0) + C0;
        //                            C_Initial.ReplaceCol(C0, 0);
        //                        }
        //                        break;
        //                    case BoundingBox.bottom:
        //                        if (border_with_vector[i].TypeBC == ABoundaryCondition.dirichlet)
        //                        {
        //                            C0 = border_with_vector[i].FunctionValues;
        //                            RVector Ctab = C_Initial.GetRowVector(0) + C0;
        //                            C_Initial.ReplaceRow(C0, 0);
        //                        }
        //                        break;
        //                }
        //            }

        //            gxt_function = g;
        //            Chat_mode = Tmode;
        //            b_filename = base_filename;
        //        }
        /// <summary>
        /// Method for solving the 2D diffusion equation using the 2D Alternating Direction Implicit algorithm
        /// </summary>
        /// <param name="n_time_steps"></param>
        //public void Solve(int n_time_steps)
        //{
        //    int nrows = cf_2D.InitialCompositionValues.GetnRows;
        //    int ncols = cf_2D.InitialCompositionValues.GetnCols;

        //    RMatrix C_Ex = new(nrows, ncols);
        //    RMatrix C_Im1 = new(nrows, ncols);
        //    RMatrix C_Im2 = new(nrows, ncols);

        //    RMatrix fn = new(nrows, ncols);
        //    RMatrix f0 = new(nrows, ncols);
        //    RMatrix f12 = new(nrows, ncols);

        //    RVector CL = new(nrows);
        //    RVector CR = new(nrows);
        //    RVector CT = new(ncols);
        //    RVector CB = new(ncols);

        //    double nu0 = dt / (2 * Math.Pow(dx, 2));

        //    //// Define the A matrix for the explicit steps
        //    //double off_d_val = nu;
        //    //double diag_val = 1 - (2 * nu);
        //    //TridiagonalMatrix A = new(ncols, diag_val, off_d_val, off_d_val);

        //    //// Define the B matrices for the implicit steps
        //    //off_d_val = -nu;
        //    //diag_val = 1 + (2 * nu);
        //    //TridiagonalMatrix B = new(nrows, diag_val, off_d_val, off_d_val);

        //    // Time evolution           
        //    for (int t = 0; t < n_time_steps; t++)
        //    {
        //        if (Chat_mode == Mode.verbose)
        //        {
        //            if (t % 10 == 0) { Console.WriteLine(t * dt); }
        //        }

        //        //if (t == n_time_steps - 1)
        //        //{
        //        //    Console.WriteLine(t);
        //        //}

        //        // ===================
        //        // Explicit time-step
        //        // ===================
        //        // 0 = top, 1 = right, 2 = left, 3 = bottom
        //        switch (BCs_Functions[1].TypeBC)
        //        {
        //            case ABoundaryCondition.dirichlet:
        //                CR = BCs_Functions[1].BoundaryFunction(t * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
        //                break;
        //            case ABoundaryCondition.neumann:
        //                RVector C1, C2; //,C3;
        //                double cUpper1, cUpper2, cLower1, cLower2;
        //                if (t == 0)
        //                {
        //                    C1 = C_Initial.GetColVector(ncols - 1);
        //                    C2 = C_Initial.GetColVector(ncols - 2);

        //                    cLower1 = C_Initial[0, ncols - 1];
        //                    cLower2 = C_Initial[1, ncols - 2];
        //                    cUpper1 = C_Initial[nrows - 1, ncols - 1];
        //                    cUpper2 = C_Initial[nrows - 2, ncols - 2];
        //                }
        //                else
        //                {
        //                    C1 = C_Im2.GetColVector(ncols - 1);
        //                    C2 = C_Im2.GetColVector(ncols - 2);

        //                    cLower1 = C_Im2[0, ncols - 1];
        //                    cLower2 = C_Im2[1, ncols - 2];
        //                    cUpper1 = C_Im2[nrows - 1, ncols - 1];
        //                    cUpper2 = C_Im2[nrows - 2, ncols - 2];
        //                }
        //                //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
        //                //CR = ((-2 * nu0) * C1) + ((2 * nu0) * C2); // + C3
        //                //CR[0] = ((-2 * nu0) * cLower1) + ((2 * nu0) * cLower2);
        //                //CR[ncols - 1] = ((-2 * nu0) * cUpper1) + ((2 * nu0) * cUpper2);

        //                CR = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(ncols - 1), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(ncols - 2), C2); // + C3
        //                CR[0] = ((-2 * nu0 * cf_2D.DValues[0, ncols - 1]) * cLower1) + ((2 * nu0 * cf_2D.DValues[1, ncols - 2]) * cLower2);
        //                CR[ncols - 1] = ((-2 * nu0 * cf_2D.DValues[nrows - 1, ncols - 1]) * cUpper1) + ((2 * nu0 * cf_2D.DValues[nrows - 2, ncols - 2]) * cUpper2);
        //                break;
        //            default:
        //                break;
        //        }
        //        switch (BCs_Functions[2].TypeBC)
        //        {
        //            case ABoundaryCondition.dirichlet:
        //                CL = BCs_Functions[2].BoundaryFunction(t * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
        //                break;
        //            case ABoundaryCondition.neumann:
        //                RVector C1, C2; //,C3;
        //                double cUpper1, cUpper2, cLower1, cLower2;
        //                if (t == 0)
        //                {
        //                    C1 = C_Initial.GetColVector(0);
        //                    C2 = C_Initial.GetColVector(1);

        //                    cLower1 = C_Initial[0, 0];
        //                    cLower2 = C_Initial[1, 1];
        //                    cUpper1 = C_Initial[nrows - 1, 0];
        //                    cUpper2 = C_Initial[nrows - 2, 1];
        //                }
        //                else
        //                {
        //                    C1 = C_Im2.GetColVector(0);
        //                    C2 = C_Im2.GetColVector(1);

        //                    cLower1 = C_Im2[0, 0];
        //                    cLower2 = C_Im2[1, 1];
        //                    cUpper1 = C_Im2[nrows - 1, 0];
        //                    cUpper2 = C_Im2[nrows - 2, 1];
        //                }
        //                //C3 = BCs[1].BoundaryFunction(t * dt, ncols);

        //                CL = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(0), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(1), C2); // + C3
        //                CL[0] = ((-2 * nu0 * cf_2D.DValues[0, 0]) * cLower1) + ((2 * nu0 * cf_2D.DValues[1, 1]) * cLower2);
        //                CL[ncols - 1] = ((-2 * nu0 * cf_2D.DValues[nrows - 1, 0]) * cUpper1) + ((2 * nu0 * cf_2D.DValues[nrows - 2, 1]) * cUpper2);
        //                break;
        //            default:
        //                break;
        //        }
        //        for (int i = 1; i < nrows - 1; i++)
        //        {
        //            RVector xold;
        //            if (t == 0) { xold = C_Initial.GetRowVector(i); }
        //            else { xold = C_Im2.GetRowVector(i); } // C_Im2.GetRowVector(i); }

        //            RVector v1 = A[i].Dot(xold);
        //            C_Ex.ReplaceRow(v1, i);
        //            C_Ex[i, 0] = CL[i]; //nu *
        //            C_Ex[i, ncols - 1] = CR[i]; //nu *
        //        }

        //        // ===================
        //        // ===================
        //        // One-half implicit time-step
        //        // ===================
        //        // Source terms
        //        double t12 = (t + 0.5) * dt;
        //        if (t == 0) { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Initial); }
        //        else { fn = gxt_function(X, Y, t12, cf_2D.DValues, C_Im2); }

        //        f12 = (dt / 2.0) * fn;

        //        // BCs
        //        switch (BCs_Functions[0].TypeBC)
        //        {
        //            case ABoundaryCondition.dirichlet:
        //                RVector gn = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
        //                RVector g0 = BCs_Functions[0].BoundaryFunction(t * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
        //                CT = ((0.5 * B_row[nrows - 1].Dot(gn)) + (0.5 * A[nrows - 1].Dot(g0))); ////((0.5 * Bm.Dot(gn)) + (0.5 * Bp.Dot(g0))); //nu * 
        //                break;
        //            case ABoundaryCondition.neumann:
        //                RVector C1, C2; //,C3;
        //                if (t == 0)
        //                {
        //                    C1 = C_Initial.GetRowVector(nrows - 1);
        //                    C2 = C_Initial.GetRowVector(nrows - 2);
        //                }
        //                else
        //                {
        //                    C1 = C_Im2.GetRowVector(nrows - 1);
        //                    C2 = C_Im2.GetRowVector(nrows - 2);
        //                }
        //                //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
        //                CT = RVector.Product(-2 * nu0 * cf_2D.DValues.GetRowVector(nrows - 1), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetRowVector(nrows - 2), C2); // + C3
        //                break;
        //            default:
        //                break;
        //        }
        //        switch (BCs_Functions[3].TypeBC)
        //        {
        //            case ABoundaryCondition.dirichlet:
        //                RVector gn = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
        //                RVector g0 = BCs_Functions[3].BoundaryFunction(t * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);
        //                CB = ((0.5 * B_row[0].Dot(gn)) + (0.5 * A[0].Dot(g0))); //((0.5 * Bm.Dot(gn)) + (0.5 * Bp.Dot(g0))); //nu * 
        //                break;
        //            case ABoundaryCondition.neumann:
        //                RVector C1, C2; //,C3;
        //                if (t == 0)
        //                {
        //                    C1 = C_Initial.GetRowVector(0);
        //                    C2 = C_Initial.GetRowVector(1);
        //                }
        //                else
        //                {
        //                    C1 = C_Im2.GetRowVector(0);
        //                    C2 = C_Im2.GetRowVector(1);
        //                }
        //                //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
        //                CB = RVector.Product(-2 * nu0 * cf_2D.DValues.GetRowVector(0), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetRowVector(1), C2); // + C3
        //                break;
        //            default:
        //                break;
        //        }
        //        for (int j = 1; j < ncols - 1; j++)
        //        {
        //            RVector v1 = C_Ex.GetColVector(j);
        //            v1[0] = CB[j]; //nu * 
        //            v1[ncols - 1] = CT[j]; //nu * 

        //            RVector f12s = f12.GetColVector(j); // 

        //            RVector u12 = TridiagonalMatrix.Thomas_Algorithm(B_row[j], v1 + f12s);
        //            C_Im1.ReplaceCol(u12, j);
        //        }

        //        // ===================
        //        // ===================
        //        // Full implicit time-step
        //        // ===================
        //        switch (BCs_Functions[1].TypeBC)
        //        {
        //            case ABoundaryCondition.dirichlet:
        //                CR = BCs_Functions[1].BoundaryFunction((t + 1) * dt, BCs_Functions[1].PositionVaries, BCs_Functions[1].PositionFixed);
        //                CT = BCs_Functions[0].BoundaryFunction((t + 1) * dt, BCs_Functions[0].PositionVaries, BCs_Functions[0].PositionFixed);
        //                RVector ctab = C_Initial.GetRowVector(0);
        //                C_Im2.ReplaceRow(CT, nrows - 1);
        //                break;
        //            case ABoundaryCondition.neumann:
        //                RVector C1, C2; //,C3;
        //                if (t == 0)
        //                {
        //                    C1 = C_Initial.GetColVector(ncols - 1);
        //                    C2 = C_Initial.GetColVector(ncols - 2);
        //                }
        //                else
        //                {
        //                    C1 = C_Im2.GetColVector(ncols - 1);
        //                    C2 = C_Im2.GetColVector(ncols - 2);
        //                }
        //                //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
        //                CR = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(ncols - 1), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(ncols - 2), C2); // + C3
        //                break;
        //            default:
        //                break;
        //        }
        //        switch (BCs_Functions[2].TypeBC)
        //        {
        //            case ABoundaryCondition.dirichlet:
        //                CL = BCs_Functions[2].BoundaryFunction((t + 1) * dt, BCs_Functions[2].PositionVaries, BCs_Functions[2].PositionFixed);
        //                CB = BCs_Functions[3].BoundaryFunction((t + 1) * dt, BCs_Functions[3].PositionVaries, BCs_Functions[3].PositionFixed);

        //                RVector ctab = C_Initial.GetRowVector(nrows - 1);

        //                C_Im2.ReplaceRow(CB, 0);
        //                break;
        //            case ABoundaryCondition.neumann:
        //                RVector C1, C2; //,C3;
        //                if (t == 0)
        //                {
        //                    C1 = C_Initial.GetColVector(0);
        //                    C2 = C_Initial.GetColVector(1);
        //                }
        //                else
        //                {
        //                    C1 = C_Im2.GetColVector(0);
        //                    C2 = C_Im2.GetColVector(1);
        //                }
        //                //C3 = BCs[1].BoundaryFunction(t * dt, ncols);
        //                CL = RVector.Product(-2 * nu0 * cf_2D.DValues.GetColVector(0), C1) + RVector.Product(2 * nu0 * cf_2D.DValues.GetColVector(1), C2); // + C3
        //                break;
        //            default:
        //                break;
        //        }
        //        for (int k = 1; k < nrows - 1; k++)
        //        {
        //            RVector v1 = C_Ex.GetRowVector(k);
        //            RVector u12 = C_Im1.GetRowVector(k);
        //            RVector b = (2 * u12) - v1;
        //            b[0] = CL[k]; //nu * 
        //            b[ncols - 1] = CR[k]; //nu * 

        //            RVector u1 = TridiagonalMatrix.Thomas_Algorithm(B_row[k], b);
        //            u1[0] = CL[k];  //nu * 
        //            u1[ncols - 1] = CR[k]; //nu * 
        //            C_Im2.ReplaceRow(u1, k);
        //        }
        //        // ===================

        //        if (Chat_mode == Mode.verbose) { if (t == n_time_steps - 1) { Console.WriteLine(t * dt); } }

        //    }

        //    // Setup the return composition field
        //    C_Final = C_Im2;
        //}

        //public DiffusionSimulators_2D_MatrixD(RMatrix D, double dx, double dy, int nx, int ny, double dt, int nt, ABoundaryCondition[] Boundary_Conditions, Del_BC_xy[] bc_s, Del_IC_xy I0, Del_Source_MatrixD g)
        //{
        //    this.dx = dx;
        //    this.dy = dy;
        //    this.dt = dt;
        //    cf_2D = new CompositionField2D(ny, nx);

        //    for (int i = 0; i < ny; i++)
        //    {
        //        for (int j = 0; j < nx; j++)
        //        {
        //            cf_2D.XPositionValues[j, i] = i * dx;
        //            cf_2D.YPositionValues[j, i] = j * dy;
        //            cf_2D.DValues[j, i] = D[j, i];
        //        }
        //    }

        //    C_Initial = I0(cf_2D.XPositionValues, cf_2D.YPositionValues);
        //    this.I0 = I0;

        //    A = new TridiagonalMatrix[nx];
        //    B_row = new TridiagonalMatrix[nx];
        //    B_col = new TridiagonalMatrix[nx];

        //    double nu0 = dt / (2 * Math.Pow(dx, 2));

        //    RVector nu = new(nx);
        //    RVector off_d_val_l0 = new(nx - 1);
        //    RVector off_d_val_u0 = new(nx - 1);
        //    RVector main0 = new(nx);
        //    RVector off_d_val_l1 = new(nx - 1);
        //    RVector off_d_val_u1 = new(nx - 1);
        //    RVector main1 = new(nx);
        //    RVector off_d_val_l2 = new(nx - 1);
        //    RVector off_d_val_u2 = new(nx - 1);
        //    RVector main2 = new(nx);

        //    double D_minus, D_plus;
        //    for (int i = 0; i < nx; i++)
        //    {
        //        RVector D_row = D.GetRowVector(i);
        //        //nu = nu0 * D_row;
        //        //main = 1 - (2 * nu);
        //        // Working this for diffusion heterogeneous media
        //        //================================================
        //        // A Matrix first
        //        //================================================
        //        D_plus = 0.5 * (D_row[1] + D_row[0]);
        //        main0[0] = 1.0 - ((D_plus * nu0));

        //        D_minus = 0.5 * (D_row[nx - 2] + D_row[0]);
        //        main0[nx - 1] = 1.0 - ((D_minus * nu0));

        //        for (int j = 1; j < nx - 1; j++)
        //        {
        //            D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
        //            D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
        //            main0[j] = 1.0 - ((D_minus * nu0) + (D_plus * nu0));
        //        }

        //        for (int j = 1; j < nx; j++)
        //        {
        //            D_minus = 0.5 * (D_row[j] + D_row[j - 1]);
        //            off_d_val_l0[j - 1] = nu0 * D_minus;
        //        }

        //        for (int j = 0; j < nx - 1; j++)
        //        {
        //            D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
        //            off_d_val_u0[j] = nu0 * D_plus;
        //        }

        //        A[i] = new TridiagonalMatrix(main0, off_d_val_l0, off_d_val_u0);
        //        ////================================================
        //        //// The B-row Matrix next
        //        ////================================================
        //        //main = 1 + (2 * nu);

        //        D_plus = 0.5 * (D_row[1] + D_row[0]);
        //        main1[0] = 1.0 + ((D_plus * nu0));

        //        D_minus = 0.5 * (D_row[nx - 2] + D_row[0]);
        //        main1[nx - 1] = 1.0 + ((D_minus * nu0));

        //        for (int j = 1; j < nx - 1; j++)
        //        {
        //            D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
        //            D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
        //            main1[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
        //        }

        //        for (int j = 1; j < nx; j++)
        //        {
        //            D_minus = 0.5 * (D_row[j - 1] + D_row[j]);
        //            off_d_val_l1[j - 1] = -nu0 * D_minus;
        //        }

        //        for (int j = 0; j < nx - 1; j++)
        //        {
        //            D_plus = 0.5 * (D_row[j + 1] + D_row[j]);
        //            off_d_val_u1[j] = -nu0 * D_plus;
        //        }

        //        B_row[i] = new TridiagonalMatrix(main1, off_d_val_l1, off_d_val_u1);
        //        ////================================================
        //        //// The B-column Matrix last
        //        ////================================================
        //        RVector D_col = D.GetColVector(i);

        //        D_plus = 0.5 * (D_col[1] + D_col[0]);
        //        main2[0] = 1.0 + ((D_plus * nu0));

        //        D_minus = 0.5 * (D_col[nx - 2] + D_col[0]);
        //        main2[nx - 1] = 1.0 + ((D_minus * nu0));

        //        for (int j = 1; j < nx - 1; j++)
        //        {
        //            D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
        //            D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
        //            main2[j] = 1.0 + ((D_minus * nu0) + (D_plus * nu0));
        //        }

        //        for (int j = 1; j < nx; j++)
        //        {
        //            D_minus = 0.5 * (D_col[j - 1] + D_col[j]);
        //            off_d_val_l2[j - 1] = -nu0 * D_minus;
        //        }

        //        for (int j = 0; j < nx - 1; j++)
        //        {
        //            D_plus = 0.5 * (D_col[j + 1] + D_col[j]);
        //            off_d_val_u2[j] = -nu0 * D_plus;
        //        }

        //        ////for (int j = 0; j < nx - 1; j++) { off_d_val_l[j] = nu[j]; }
        //        ////for (int j = 1; j < nx; j++) { off_d_val_u[j - 1] = nu[j]; }

        //        B_col[i] = new TridiagonalMatrix(main2, off_d_val_l2, off_d_val_u2);
        //    }

        //    C_Initial = I0(cf_2D.XPositionValues, cf_2D.YPositionValues);
        //    this.I0 = I0;

        //    int num_bounds = Boundary_Conditions.Length;
        //    border_with_function = new BoundaryWithFunction[num_bounds];
        //    RVector C0;
        //    for (int i = 0; i < num_bounds; i++)
        //    {
        //        border_with_function[i] = new BoundaryWithFunction
        //        {
        //            BoundaryLocation = i switch
        //            {
        //                0 => BoundingBox.top,
        //                1 => BoundingBox.right,
        //                2 => BoundingBox.left,
        //                3 => BoundingBox.bottom,
        //                _ => BoundingBox.bottom,
        //            },
        //            TypeBC = Boundary_Conditions[i],
        //            BoundaryFunction = bc_s[i]
        //        };
        //        switch (border_with_function[i].BoundaryLocation)
        //        {
        //            case BoundingBox.top:
        //                if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
        //                {
        //                    border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
        //                    border_with_function[i].PositionFixed = Y[X.GetnRows - 1, 0];
        //                    C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

        //                    RVector Ctab = C_Initial.GetRowVector(ny - 1) + C0;
        //                    C_Initial.ReplaceRow(C0, ny - 1);
        //                }
        //                break;
        //            case BoundingBox.right:
        //                if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
        //                {
        //                    border_with_function[i].PositionVaries = Y.GetColVector(0);
        //                    border_with_function[i].PositionFixed = X[0, Y.GetnCols - 1];
        //                    C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

        //                    RVector Ctab = C_Initial.GetColVector(nx - 1) + C0;
        //                    C_Initial.ReplaceCol(C0, nx - 1);
        //                }
        //                break;
        //            case BoundingBox.left:
        //                if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
        //                {
        //                    border_with_function[i].PositionVaries = Y.GetColVector(0);
        //                    border_with_function[i].PositionFixed = X[0, 0];
        //                    C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

        //                    RVector Ctab = C_Initial.GetColVector(0) + C0;
        //                    C_Initial.ReplaceCol(C0, 0);
        //                }
        //                break;
        //            case BoundingBox.bottom:
        //                if (border_with_function[i].TypeBC == ABoundaryCondition.dirichlet)
        //                {
        //                    border_with_function[i].PositionVaries = X.GetRowVector(X.GetnRows - 1);
        //                    border_with_function[i].PositionFixed = Y[0, 0];
        //                    C0 = border_with_function[i].BoundaryFunction(0.0, border_with_function[i].PositionVaries, border_with_function[i].PositionFixed);

        //                    RVector Ctab = C_Initial.GetRowVector(0) + C0;
        //                    C_Initial.ReplaceRow(C0, 0);
        //                }
        //                break;
        //        }
        //    }

        //    gxt_function = g;
        //    Chat_mode = Mode.quiet;
        //    Errors = new string[nt];
        //}

    }
}
