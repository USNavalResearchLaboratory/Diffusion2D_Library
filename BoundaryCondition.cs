using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Diffusion2D_Library
{
    public class BoundaryCondition
    {
        public enum ABoundaryCondition { dirichlet, neumann };

        // =======================================================================
        // Private methods
        // =======================================================================
        /// <summary>
        /// Converts the enumerated boundary condition to a string
        /// </summary>
        /// <param name="bc">variable of type BoundaryCondition</param>
        /// <returns></returns>
        public static string ConvertEnumBCToString(ABoundaryCondition bc)
        {
            string result = bc switch
            {
                ABoundaryCondition.dirichlet => "Dirichlet",
                ABoundaryCondition.neumann => "Neumann",
                _ => "Dirichlet",
            };
            return result;
        }
        /// <summary>
        /// Method to convert a string value to one of the enumerated boundary conditions
        /// </summary>
        /// <param name="value">string variable to be converted to the boundary condition enum</param>
        /// <returns></returns>
        public static ABoundaryCondition ConvertStringToEnumBC(string value)
        {
            ABoundaryCondition bc = value.ToLower() switch
            {
                "dirichlet" => ABoundaryCondition.dirichlet,
                "neumann" => ABoundaryCondition.neumann,
                _ => ABoundaryCondition.dirichlet,
            };
            return bc;
        }

        // =======================================================================
    }
}
