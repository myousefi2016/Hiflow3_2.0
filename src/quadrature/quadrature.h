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

#ifndef __QUADRATURE_QUADRATURE_H_
#    define __QUADRATURE_QUADRATURE_H_

#    include <cassert>
#    include <iostream>
#    include <string>
#    include <vector>

#    include "quadrature/quadraturetype.h"

namespace hiflow
{

    ///
    /// \class Quadrature quadrature.h
    /// \brief Holds all necessary information about the desired quadrature rule
    ///        and is templated by a DataType (e.g. double)
    /// \author Michael Schick
    ///

    template<class DataType>
    class Quadrature
    {
      public:

        /// Default constructor (setting all pointers to NULL)
        Quadrature ( );
        /// Default destructor (clearing all pointers)
        ~Quadrature ( );
        /// Copy constructor (deep copying object)
        Quadrature ( const Quadrature<DataType>& q );
        /// Assignment operator (deep copying object)
        const Quadrature<DataType>& operator= ( const Quadrature<DataType>& q );

        /// Returns the size (number of quadrature points) of the quadrature

        int size ( ) const
        {
            return size_;
        }
        /// Returns the name of the used quadrature

        const std::string& name ( ) const
        {
            return name_;
        }

        /// Setting up which quadrature rule should be used, e.g. GaussHexahedron
        /// and setting up the desired size
        void set_quadrature ( const std::string& name, int size );

        /// \brief Sets the quadrature rule by the desired order of accuracy.
        void set_quadrature_by_order ( const std::string& name, int order );

        /// Set up facet quadrature rule by mapping a lower-dimensional
        /// quadrature rule to a higher-dimensional cell.
        void set_facet_quadrature ( const Quadrature<DataType>& base_quadrature,
                                    int cell_type,
                                    int facet_number );

        void set_custom_quadrature ( const std::vector<DataType>& x_coords,
                                     const std::vector<DataType>& y_coords,
                                     const std::vector<DataType>& z_coords,
                                     const std::vector<DataType>& weights );
        /// Set cell type
        inline void set_cell_type ( int cell_type );

        // Get cell type

        int get_cell_type ( ) const
        {
            return cell_type_;
        }

        /// Get the x value of the quadrature point with number qpt_index
        inline DataType x ( int qpt_index ) const;
        /// Get the y value of the quadrature point with number qpt_index
        inline DataType y ( int qpt_index ) const;
        /// Get the z value of the quadrature point with number qpt_index
        inline DataType z ( int qpt_index ) const;
        /// Get the weight value of the quadrature point with number qpt_index
        inline DataType w ( int qpt_index ) const;

        /// Print some status information about the used quadrature
        void print_status ( ) const;

      protected:

        /// Reseting all quadrature pointers to NULL
        void clear ( );

        /// Name of quadrature
        std::string name_;
        /// Size of quadrature, i.e. number of quadrature points
        int size_;

        /// Cell Type of the quadrature
        int cell_type_;

        /// Holds the information about the used quadrature
        const QuadratureType<DataType>* quad_;

        /// Pointer to the x - values of the quadrature points
        const std::vector<DataType>* qpts_x_;
        /// Pointer to the y - values of the quadrature points
        const std::vector<DataType>* qpts_y_;
        /// Pointer to the z - values of the quadrature points
        const std::vector<DataType>* qpts_z_;
        /// Pointer to the weight values of the quadrature points
        const std::vector<DataType>* weights_;
    };

    // INLINE FUNCTIONS FOR QUADRATURE

    template<class DataType>
    void Quadrature<DataType>::set_cell_type ( int cell_type )
    {
        cell_type_ = cell_type;
    }

    template<class DataType>
    DataType Quadrature<DataType>::x ( int qpt_index ) const
    {
        return (*qpts_x_ )[qpt_index];
    }

    template<class DataType>
    DataType Quadrature<DataType>::y ( int qpt_index ) const
    {
        return (*qpts_y_ )[qpt_index];
    }

    template<class DataType>
    DataType Quadrature<DataType>::z ( int qpt_index ) const
    {
        return (*qpts_z_ )[qpt_index];
    }

    template<class DataType>
    DataType Quadrature<DataType>::w ( int qpt_index ) const
    {
        return (*weights_ )[qpt_index];
    }

} // namespace hiflow

#endif
