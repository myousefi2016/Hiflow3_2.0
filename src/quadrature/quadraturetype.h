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

#ifndef __QUADRATURE_QUADRATURETYPE_H_
#    define __QUADRATURE_QUADRATURETYPE_H_

#    include <vector>
#    include <map>
#    include <iostream>

namespace hiflow
{

    ///
    /// \class QuadratureType quadraturetype.h
    /// \brief Ancestor class of all implemented quadrature types
    ///        which is templated by a DataType (e.g. double)
    /// \author Michael Schick
    ///

    template<class DataType>
    class QuadratureType
    {
      public:

        /// Default constructor
        QuadratureType ( );

        virtual ~QuadratureType ( )
        {
        }

        /// Polymorphic cloning function, returning a deep copy of itself.
        virtual QuadratureType<DataType>* clone ( ) const = 0;

        /// Returns const pointer to the stored x values of the quadrature
        const std::vector<DataType>* x ( int size ) const;
        /// Returns const pointer to the stored y values of the quadrature
        const std::vector<DataType>* y ( int size ) const;
        /// Returns const pointer to the stored z values of the quadrature
        const std::vector<DataType>* z ( int size ) const;
        /// Returns const pointer to the stored weight values of the quadrature
        const std::vector<DataType>* w ( int size ) const;

        /// Returns the size for the quadrature of the requested order.
        int size_for_order ( int order ) const;

      protected:

        /// Find the corresponding vector index for the desired size
        int find_quadrature ( int size ) const;

        /// Stored x values of the quadrature
        std::vector<std::vector<DataType> > x_;
        /// Stored y values of the quadrature
        std::vector<std::vector<DataType> > y_;
        /// Stored z values of the quadrature
        std::vector<std::vector<DataType> > z_;
        /// Stored weight values of the quadrature
        std::vector<std::vector<DataType> > w_;

        /// Mapping for finding the correct vector index for a desired quadrature size
        std::map<int, int> index_field_;

        /// Map from order to size
        std::vector<int> order_size_map_;
    };

} // namespace hiflow

#endif
