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

#include "quadrature.h"

#include "custom_quadrature_type.h"
#include "qgausstetrahedron.h"
#include "qgausshexahedron.h"
#include "qgausstriangle.h"
#include "qgaussquadrilateral.h"
#include "qgaussline.h"
#include "qeconomicalgaussquadrilateral.h"
#include "qeconomicalgausshexahedron.h"
#include "mapped_quadrature_type.h"
#include "qgausspyramid.h"

#include <sstream>

namespace hiflow
{

    /// Factory function for quadrature types.

    template<class DataType>
    QuadratureType<DataType>* create_quadrature_type ( const std::string& name )
    {
        if ( name == "GaussTetrahedron" )
        {
            return new QuadratureGaussTetrahedron<DataType>;
        }
        else if ( name == "GaussHexahedron" )
        {
            return new QuadratureGaussHexahedron<DataType>;
        }
        else if ( name == "EconomicalGaussHexahedron" )
        {
            return new QuadratureEconomicalGaussHexahedron<DataType>;
        }
        else if ( name == "GaussTriangle" )
        {
            return new QuadratureGaussTriangle<DataType>;
        }
        else if ( name == "GaussQuadrilateral" )
        {
            return new QuadratureGaussQuadrilateral<DataType>;
        }
        else if ( name == "EconomicalGaussQuadrilateral" )
        {
            return new QuadratureEconomicalGaussQuadrilateral<DataType>;
        }
        else if ( name == "GaussLine" )
        {
            return new QuadratureGaussLine<DataType>;
        }
        else if ( name == "GaussPyramid" )
        {
            return new QuadratureGaussPyramid<DataType>;
        }
        else
        {
            assert ( 0 );
        }
        return 0;
    }

    template<class DataType>
    Quadrature<DataType>::Quadrature ( )
    : name_ ( "" ), size_ ( 0 ), cell_type_ ( -1 ), quad_ ( 0 ), qpts_x_ ( 0 ), qpts_y_ ( 0 ), qpts_z_ ( 0 ), weights_ ( 0 )
    {
    }

    template<class DataType>
    Quadrature<DataType>::~Quadrature ( )
    {
        clear ( );
    }

    template<class DataType>
    Quadrature<DataType>::Quadrature ( const Quadrature<DataType>& q )
    : name_ ( q.name ( ) ), size_ ( q.size ( ) ), quad_ ( 0 )
    {
        if ( q.quad_ )
        {
            quad_ = q.quad_->clone ( );

            qpts_x_ = quad_->x ( size_ );
            qpts_y_ = quad_->y ( size_ );
            qpts_z_ = quad_->z ( size_ );
            weights_ = quad_->w ( size_ );
        }
    }

    template<class DataType>
    const Quadrature<DataType>& Quadrature<DataType>::operator= ( const Quadrature<DataType>& q )
    {
        if ( &q != this )
        {
            clear ( );

            size_ = q.size ( );
            name_ = q.name ( );

            quad_ = q.quad_->clone ( );

            qpts_x_ = quad_->x ( size_ );
            qpts_y_ = quad_->y ( size_ );
            qpts_z_ = quad_->z ( size_ );
            weights_ = quad_->w ( size_ );
        }
        return *this;
    }

    template<class DataType>
    void Quadrature<DataType>::clear ( )
    {
        if ( quad_ )
        {
            delete quad_;
        }
        quad_ = 0;
        qpts_x_ = 0;
        qpts_y_ = 0;
        qpts_z_ = 0;
        weights_ = 0;
    }

    template<class DataType>
    void Quadrature<DataType>::set_quadrature ( const std::string& name, int size )
    {
        clear ( );
        quad_ = create_quadrature_type<DataType>( name );

        assert ( quad_ != 0 );

        qpts_x_ = quad_->x ( size );
        qpts_y_ = quad_->y ( size );
        qpts_z_ = quad_->z ( size );
        weights_ = quad_->w ( size );

        assert ( qpts_x_ );
        assert ( qpts_y_ );
        assert ( qpts_z_ );
        assert ( weights_ );

        size_ = size;
        name_ = name;
    }

    template<class DataType>
    void Quadrature<DataType>::set_quadrature_by_order ( const std::string& name, int order )
    {
        clear ( );

        quad_ = create_quadrature_type<DataType>( name );

        assert ( quad_ != 0 );

        const int size = quad_->size_for_order ( order );

        qpts_x_ = quad_->x ( size );
        qpts_y_ = quad_->y ( size );
        qpts_z_ = quad_->z ( size );
        weights_ = quad_->w ( size );

        assert ( qpts_x_ );
        assert ( qpts_y_ );
        assert ( qpts_z_ );
        assert ( weights_ );

        size_ = size;
        name_ = name;
    }

    template<class DataType>
    void Quadrature<DataType>::set_facet_quadrature ( const Quadrature<DataType>& base_quadrature,
                                                      int cell_type,
                                                      int facet_number )
    {
        clear ( );

        const std::string base_name = base_quadrature.name ( );
        const int base_size = base_quadrature.size ( );

        QuadratureMapping<DataType>* mapping;

        switch ( cell_type )
        {
            case 2: // Triangle
                mapping = new TriangleQuadratureMapping<DataType>( facet_number );
                break;
            case 3: // Quadrilateral
                mapping = new QuadrilateralQuadratureMapping<DataType>( facet_number );
                break;
            case 4: // Tetrahedron
                mapping = new TetrahedralQuadratureMapping<DataType>( facet_number );
                break;
            case 5: // Hexahedron
                mapping = new HexahedralQuadratureMapping<DataType>( facet_number );
                break;
            case 6: // Pyramid
                mapping = new PyramidQuadratureMapping<DataType>( facet_number );
                break;
            default:
                // TODO: not yet implemented
                std::cerr << "No mapping implemented for cell type " << cell_type << ".\n";
                assert ( false );
        };

        assert ( mapping != 0 );
        quad_ = new MappedQuadratureType<DataType>( base_quadrature, *mapping );
        size_ = base_size;
        qpts_x_ = quad_->x ( size_ );
        qpts_y_ = quad_->y ( size_ );
        qpts_z_ = quad_->z ( size_ );
        weights_ = quad_->w ( size_ );
        assert ( qpts_x_ != 0 );
        assert ( qpts_y_ != 0 );
        assert ( qpts_z_ != 0 );
        assert ( weights_ != 0 );
        // create unique name for quadrature
        std::stringstream name_stream;
        name_stream << base_name << "_" << base_size
                << "_Mapped_" << cell_type
                << "_Facet_" << facet_number;
        name_ = name_stream.str ( );
        if ( mapping )
        {
            delete mapping;
        }
    }

    template<class DataType>
    void Quadrature<DataType>::set_custom_quadrature ( const std::vector<DataType>& x_coords,
                                                       const std::vector<DataType>& y_coords,
                                                       const std::vector<DataType>& z_coords,
                                                       const std::vector<DataType>& weights )
    {
        clear ( );

        size_ = weights.size ( );

        assert ( x_coords.size ( ) == size_ );
        assert ( y_coords.size ( ) == size_ );
        assert ( z_coords.size ( ) == size_ );

        quad_ = new CustomQuadratureType<DataType>( x_coords, y_coords, z_coords, weights );

        qpts_x_ = quad_->x ( size_ );
        qpts_y_ = quad_->y ( size_ );
        qpts_z_ = quad_->z ( size_ );
        weights_ = quad_->w ( size_ );

        assert ( qpts_x_ != 0 );
        assert ( qpts_y_ != 0 );
        assert ( qpts_z_ != 0 );
        assert ( weights_ != 0 );

        name_ = "Custom";
    }

    template<class DataType>
    void Quadrature<DataType>::print_status ( ) const
    {
        std::cout << "Quadrature points for class Quadrature" << std::endl;

        std::cout << "Name of quadrature:\t" << name_ << std::endl;
        std::cout << "Size of quadrature:\t" << size_ << std::endl;

        if ( qpts_x_ != 0 )
        {
            std::cout << "x values:" << std::endl;
            for ( int j = 0; j < size_; ++j )
                std::cout << qpts_x_->at ( j ) << "\t";
            std::cout << std::endl;
            std::cout << "-------------------------------------" << std::endl;
        }
        else
            std::cout << "No x values defined for this quadrature!" << std::endl;

        if ( qpts_y_ != 0 )
        {
            std::cout << "y values:" << std::endl;
            for ( int j = 0; j < size_; ++j )
                std::cout << qpts_y_->at ( j ) << "\t";
            std::cout << std::endl;
            std::cout << "-------------------------------------" << std::endl;
        }
        else
            std::cout << "No y values defined for this quadrature!" << std::endl;

        if ( qpts_z_ != 0 )
        {
            std::cout << "z values:" << std::endl;
            for ( int j = 0; j < size_; ++j )
                std::cout << qpts_z_->at ( j ) << "\t";
            std::cout << std::endl;
            std::cout << "-------------------------------------" << std::endl;
        }
        else
            std::cout << "No z values defined for this quadrature!" << std::endl;

        if ( weights_ != 0 )
        {
            std::cout << "weight values:" << std::endl;
            for ( int j = 0; j < size_; ++j )
                std::cout << weights_->at ( j ) << "\t";
            std::cout << std::endl;
            std::cout << "-------------------------------------" << std::endl;
        }
        else
            std::cout << "WARNING: No weights defined for this quadrature!" << std::endl;
    }

    // template instanciation
    template class Quadrature<double>;
    template class Quadrature<float>;

} // namespace hiflow
