// Copyright (C) 2011-2017 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for more
// details.
//
// You should have received a copy of the European Union Public Licence (EUPL) v1.2
// along with HiFlow3.  If not, see <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

/// \author Teresa Beck, Philipp Gerstner

#ifndef _CYL_CELL_VISUALIZATION_H_
#    define _CYL_CELL_VISUALIZATION_H_

#    include <map>
#    include <vector>
#    include <string>
#    include <cmath>
#    include <limits>
#    include <sstream>

#    include "space/cell_visualization.h"
#    include "common/log.h"
#    include "fem/cell_transformation.h"
#    include "mesh/entity.h"
#    include "mesh/iterator.h"
#    include "space/vector_space.h"

#    include "../met_flow_main.h"
#    include "eval_cyl_fe.h"

#    include "../../../../contrib/tinyxml/tinyxml.h"

namespace hiflow
{

    /// \brief Description of a square 2d grid or cubic 3d grid.
    ///
    /// \todo Use of variable grid size in each direction.
    /// \todo Visualization for simplex grids.

    class RegularGrid
    {
      public:
        RegularGrid ( int gdim, int num_intervals, double origin, double side_length );
        int num_points ( ) const;
        int num_cells ( ) const;
        const std::vector<int>& vertices_of_cell ( int i ) const;
        const std::vector<double>& coords ( ) const;
        int gdim ( ) const;
      private:
        void init2d ( double origin, double dx );
        void init3d ( double origin, double dx );
        int index ( int i ) const;
        int index ( int i, int j ) const;
        int index ( int i, int j, int k ) const;
        int gdim_;
        int num_intervals_;
        std::vector< std::vector<int> > cells_;
        std::vector<double> coords_;
    };
}

namespace hiflow
{
    namespace mesh
    {

        /// \brief Maps parameter space to geometrical space
        /// \param[inout] vector with coordinates

        inline void parameter_space_to_geometrical_space ( std::vector<double>& coords,
                                                           double parameter_r,
                                                           double parameter_R )
        {
            // parameter space, e.g. (0..2pi)
            const double pi = M_PI; //acos(-1.);
            double phi_min = 0.; // e1-min in parameter space
            double phi_max = 2 * pi; // e1-max in parameter space

            // physical space, e.g. sphere
            double phi = coords.at ( 0 ) / ( phi_max - phi_min )*( 2. * pi );
            double r = coords.at ( 1 ); ///(  r_max- r_min)*(parameter_R-parameter_r);

            double x = r * cos ( phi );
            double y = r * sin ( phi );

            coords.at ( 0 ) = x;
            coords.at ( 1 ) = y;

            //      assert(coords.at(2) == 0);
        }

        class CylCellVisualization : public ParallelCellVisualization<DATATYPE>
        {
            /// \brief CellVisualization to be used for periodic domains
            /// \see CellVisualization
          private:
            const VectorSpace<DATATYPE>& space_;
            RegularGrid grid_;
            std::map< std::string, std::vector<double> > functions_;
            std::map< std::string, int > num_components_;
          public:

            CylCellVisualization ( const VectorSpace<DATATYPE>& space, int num_intervals, const MPI_Comm& comm, int master_rank ) : ParallelCellVisualization<DATATYPE> ( space, num_intervals, comm, master_rank ), space_ ( space ), grid_ ( space.mesh ( ).gdim ( ), num_intervals, 0., 1. )
            {
            };
            void write_sequential ( const std::string& filename ) const;
            void write ( const std::string& filename ) const;
            //#ifdef USE_HYPRE
            //      void visualize(const EvalFeFunctionHypre<LAD>& fun, const std::string& name, int num_components);
            //#else
            void visualize ( const EvalFeFunction<LAD>& fun, const std::string& name, int num_components );
            //#endif
            void visualize ( const EvalCylFeFunction& fun, const std::string& name, int num_components );

            void visualize_cell_data ( const std::vector<double>& cell_data, const std::string& name )
            {
                return CellVisualization<DATATYPE>::visualize_cell_data ( cell_data, name );

            };
        };
    };
};

#endif
