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

#include "fetype.h"

#include <cmath>

#include "common/macros.h"
#include "common/vector_algebra.h"

#include "../mesh/refinement.h"
#include "common/log.h"

/// \author Michael Schick<br>Martin Baumann

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        FEType<DataType>::FEType ( )
        : init_status_ ( false ), ref_cell_ ( 0 ), tdim_ ( 0 ), fe_deg_ ( -1 )
        {
        }

        template<class DataType>
        FEType<DataType>::~FEType ( )
        {
        }

        template<class DataType>
        int FEType<DataType>::get_nb_dof_on_cell ( ) const
        {
            return coord_.size ( );
        }

        template<class DataType>
        int FEType<DataType>::get_nb_subentity ( int tdim ) const
        {
            assert ( tdim >= 0 && tdim < ref_cell_->tdim ( ) );
            return coord_on_subentity_[tdim].size ( );
        }

        template<class DataType>
        int FEType<DataType>::get_nb_dof_on_subentity ( int tdim, int index ) const
        {
            assert ( coord_vertices_.size ( ) > 0 );
            assert ( tdim >= 0 && tdim < ref_cell_->tdim ( ) );
            assert ( index >= 0 && index < coord_on_subentity_[tdim].size ( ) );
            return coord_on_subentity_[tdim][index].size ( );
        }

        template<class DataType>
        std::vector<std::vector<DataType> > const& FEType<DataType>::get_coord_on_subentity ( int tdim, int index ) const
        {
            assert ( coord_vertices_.size ( ) > 0 );
            assert ( tdim >= 0 && tdim < ref_cell_->tdim ( ) );
            assert ( index >= 0 && index < coord_on_subentity_[tdim].size ( ) );
            return coord_on_subentity_[tdim][index];
        }

        template<class DataType>
        std::vector<DofID> const& FEType<DataType>::get_dof_on_subentity ( int tdim, int index ) const
        {
            assert ( coord_vertices_.size ( ) > 0 );
            assert ( tdim >= 0 && tdim < ref_cell_->tdim ( ) );
            assert ( index >= 0 && index < dof_on_subentity_[tdim].size ( ) );
            return dof_on_subentity_[tdim][index];
        }

        template<class DataType>
        std::vector<std::vector<DataType> > const& FEType<DataType>::get_coord_on_cell ( ) const
        {
            return coord_;
        }

        template<class DataType>
        void FEType<DataType>::init ( )
        {
            // calculate the coordinates of the DoF points
            this->init_coord ( );

            // tdim_ must be set by init_coord()
            assert ( tdim_ != 0 );

            // filter the DoF points for the cell's faces, edges or vertices
            if ( coord_vertices_.size ( ) > 0 )
                this->init_dofs_on_subentities ( );

            init_status_ = true;
        }

        template<class DataType>
        bool FEType<DataType>::operator== ( const FEType<DataType>& fe_slave ) const
        {
            if ( this->my_id_ == 0 || fe_slave.my_id_ == 0 )
                assert ( 0 );

            if ( this->my_id_ == fe_slave.my_id_ &&
                 this->fe_deg_ == fe_slave.fe_deg_ )
                return true;
            else
                return false;
        }

        template<class DataType>
        bool FEType<DataType>::operator< ( const FEType<DataType>& fe_slave ) const
        {
            if ( this->my_id_ == 0 || fe_slave.my_id_ == 0 )
                assert ( 0 );

            if ( this->my_id_ < fe_slave.my_id_ )
                return true;
            else if ( this->my_id_ == fe_slave.my_id_ )
            {
                if ( this->fe_deg_ < fe_slave.fe_deg_ )
                    return true;
                else
                    return false;
            }
            else
                return false;
        }

        template<class DataType>
        void FEType<DataType>::set_ref_vtx_coords ( std::vector<Coord> ref_coords )
        {
            for ( size_t i = 0, e_i = ref_coords.size ( ); i != e_i; ++i )
                coord_vertices_.push_back ( ref_coords[i] );

        }

        /// \details This function should only be used once. This is made sure by
        ///          checking (assert) that the fe_deg_ is not set (equals to -1)
        ///          before executing the rest of the function. It is not the intention
        ///          of the authors that this function is used for p refinement, because
        ///          of the structure in FEInstance \see FEInstance

        template<class DataType>
        void FEType<DataType>::set_fe_deg ( int degree )
        {
            assert ( fe_deg_ == -1 );

            fe_deg_ = degree;
        }

        /// This function filters all DoFs and fills the data structures
        /// coord_on_subentity_ and dof_on_subentity such that it is possible
        /// to find out, f.e. which DoF IDs are located on the cell's second
        /// edge. The numbering of the subentities is the one given by the
        /// mesh::CellType. Both data structures use the following indices:
        /// info_on_subentity[tdim][index_subentity][index_point], where
        /// tdim =  0..dim-1 (0->point, 1->edge, 2->face) <br>
        /// Example: dof_on_subentity_[1][3][2] will return the point 2 of
        ///          the cells subentity 3, which is an edge (-> 1).

        template<class DataType>
        void FEType<DataType>::init_dofs_on_subentities ( )
        {

            int verbatim_mode = 0; // 0 -> no output
            // 1 -> some output

            // clear data structures

            coord_on_subentity_.clear ( );
            coord_on_subentity_.resize ( ref_cell_->tdim ( ) );

            dof_on_subentity_.clear ( );
            dof_on_subentity_.resize ( ref_cell_->tdim ( ) );

            assert ( coord_on_subentity_.size ( ) == dof_on_subentity_.size ( ) );
            assert ( coord_on_subentity_.size ( ) == ref_cell_->tdim ( ) );

            for ( int tdim = 0; tdim < ref_cell_->tdim ( ); ++tdim )
            {
                coord_on_subentity_[tdim].resize ( ref_cell_->num_entities ( tdim ) );
                dof_on_subentity_[tdim].resize ( ref_cell_->num_entities ( tdim ) );
            }

            // get point coordinates for this entity, including refined entities

            int gdim = coord_vertices_[0].size ( );
            std::vector<double> cv_; // coord_vertices_ in other numeration
            for ( size_t p = 0, e_p = coord_vertices_.size ( ); p != e_p; ++p )
                for ( size_t c = 0, e_c = gdim; c != e_c; ++c )
                    cv_.push_back ( coord_vertices_[p][c] );

            std::vector<double> coord_vertices_refined;
            mesh::compute_refined_vertices ( *ref_cell_, gdim, cv_, coord_vertices_refined );

            if ( verbatim_mode > 0 )
            {
                for ( size_t i = 0, e_i = coord_vertices_refined.size ( ); i != e_i; ++i )
                    std::cout << "\t" << coord_vertices_refined[i] << std::endl;
            }

            // insert filtered DoFs

            /// \todo also tdim = ref_cell_->tdim() should be available
            /// \todo coord_ could also be stored by coord_on_subentity_

            for ( size_t tdim = 0, e_tdim = ref_cell_->tdim ( ); tdim != e_tdim; ++tdim )
            {
                // for output purposes

                std::string entity_type;
                if ( tdim == 0 ) entity_type = "point";
                if ( tdim == 1 ) entity_type = "line";
                if ( tdim == 2 ) entity_type = "face";
                if ( tdim == 3 ) entity_type = "volume";

                for ( size_t idx_entity = 0, e_idx_entity = ref_cell_->num_entities ( tdim ); idx_entity != e_idx_entity; ++idx_entity )
                {

                    if ( verbatim_mode > 0 )
                        std::cout << "DoF points on " << entity_type << " " << idx_entity << ":" << std::endl;

                    // get point coordinates for this entity

                    std::vector<Coord> points;
                    std::vector<int> vertex = ref_cell_->local_vertices_of_entity ( tdim, idx_entity );

                    // points can't be handled using local_vertices_of_entity()

                    if ( tdim == 0 )
                    {
                        vertex.clear ( );
                        vertex.push_back ( idx_entity );
                    }

                    // fill list of points that define the entity
                    for ( size_t point = 0, e_point = vertex.size ( ); point != e_point; ++point )
                    {
                        Coord temp;
                        for ( int d = 0; d < gdim; ++d )
                            temp.push_back ( coord_vertices_refined[vertex[point] * gdim + d] );
                        points.push_back ( temp );
                    }

                    if ( verbatim_mode > 0 )
                    {
                        std::cout << "defining points: " << std::endl;
                        for ( int p = 0; p < points.size ( ); ++p )
                        {
                            for ( int d = 0; d < points[0].size ( ); ++d )
                                std::cout << "\t" << points[p][d];
                            std::cout << "\t --- ";
                        }
                        std::cout << std::endl;
                    }

                    // filter points that are on subentity
                    if ( verbatim_mode > 0 )
                    {
                        std::cout << "Filtering points on subentity (dim = " << tdim << ", idx = " << idx_entity << ")\n";
                    }

                    for ( size_t dof_point = 0, e_dof_point = coord_.size ( ); dof_point != e_dof_point; ++dof_point )
                    {

                        Coord dof_coord = coord_[dof_point];

                        // filter DoF
                        if ( is_point_on_subentity ( dof_coord, points ) == true )
                        {
                            coord_on_subentity_[tdim][idx_entity].push_back ( dof_coord );
                            dof_on_subentity_[tdim][idx_entity].push_back ( dof_point );

                            if ( verbatim_mode > 0 )
                            {
                                std::cout << dof_point << "\t-> ";
                                for ( int d = 0; d < dof_coord.size ( ); ++d )
                                    std::cout << dof_coord[d] << "\t";
                                std::cout << std::endl;
                            }
                        } // if (is_point_on_subentity(dof_coord, points) == true)
                        else
                        {
                            if ( verbatim_mode > 0 )
                            {
                                std::cout << dof_point << " is not on subentity.\n";
                            }
                        }
                        if ( verbatim_mode > 0 )
                        {
                            std::cout << "\n\n";
                        }
                    } // for (int dof_point=0; dof_point<coord_.size(); ++dof_point)
                } // for (int idx_entity=0; idx_entity<ref_cell_->num_entities(tdim); ++idx_entity)

            } // for (int tdim=0; tdim<ref_cell_->tdim(); ++tdim)
            if ( verbatim_mode > 0 )
            {
                std::cout << "\n\n";
            }

            // cell type depending tests

            if ( ( ref_cell_->tag ( ) == mesh::CellType::LINE ) ||
                 ( ref_cell_->tag ( ) == mesh::CellType::TRIANGLE ) ||
                 ( ref_cell_->tag ( ) == mesh::CellType::QUADRILATERAL ) ||
                 ( ref_cell_->tag ( ) == mesh::CellType::TETRAHEDRON ) ||
                 ( ref_cell_->tag ( ) == mesh::CellType::HEXAHEDRON ) ||
                 ( ref_cell_->tag ( ) == mesh::CellType::PYRAMID ) )
            {

                // 1. Test: For the listed cell types, a test is that the number of DoFs
                //          lying on each facet should be the same

                bool ok = true;

                // pyramid element's first facet has one more dof than the others
                int number_dofs;
                if ( ref_cell_->tag ( ) == mesh::CellType::PYRAMID )
                {
                    number_dofs = dof_on_subentity_[ref_cell_->tdim ( ) - 1][1].size ( );
                }
                else
                {
                    number_dofs = dof_on_subentity_[ref_cell_->tdim ( ) - 1][0].size ( );
                }

                for ( int i = 1, e_i = ref_cell_->num_regular_entities ( ref_cell_->tdim ( ) - 1 ); i < e_i; ++i )
                {
                    if ( number_dofs != dof_on_subentity_[ref_cell_->tdim ( ) - 1][i].size ( ) )
                    {

                        ok = false;
                        std::cerr << "Only " << dof_on_subentity_[ref_cell_->tdim ( ) - 1][i].size ( ) << " dofs on facet " << i
                                << ", expected " << number_dofs << "\n";
                        quit_program ( );
                    }
                }

                if ( ok == false )
                {
                    std::cout << "Error: number of subentities wrong!" << std::endl
                            << "       'eps' parameter in is_point_on_subentity() might be "
                            << "to small" << std::endl;
                    quit_program ( );
                }

                // 2. Test: The number of DoFs lying on subfacets should be checked as
                //          well.

                /// \todo Write test.
            }
            else
            {
                // here could stand some other criterion ...
                quit_program ( );
            }

            // print summary
            if ( verbatim_mode > 0 )
            {
                std::cout << "Summary for fe deg " << fe_deg_ << ": " << std::endl;
                std::cout << "========" << std::endl;
                for ( int tdim = 0; tdim < ref_cell_->tdim ( ); ++tdim )
                {
                    for ( int idx_entity = 0; idx_entity < dof_on_subentity_[tdim].size ( ); ++idx_entity )
                    {
                        std::cout << "TDim " << tdim << ", " << idx_entity << " ";
                        if ( idx_entity >= ref_cell_->num_regular_entities ( tdim ) )
                            std::cout << "(refined entity) ";
                        else
                        {
                            std::cout << "(regular entity) ";
                        }
                        std::cout << dof_on_subentity_[tdim][idx_entity].size ( )
                                << " DoFs:" << std::endl;
                        std::cout << string_from_range ( dof_on_subentity_[tdim][idx_entity].begin ( ), dof_on_subentity_[tdim][idx_entity].end ( ) ) << "\n\n";
                    }
                    std::cout << "========" << std::endl;
                }
                std::cout << std::endl;
            }
        }

        template<class DataType>
        bool FEType<DataType>::is_point_on_subentity ( Coord point, const std::vector<Coord>& points ) const
        {

            int verbatim_mode = 0; // 0 -> no output
            // 1 -> some output

            if ( verbatim_mode > 0 )
            {
                std::cout << "Is point (";
                for ( int k = 0; k < point.size ( ); ++k )
                {
                    std::cout << point[k] << " ";
                }
                std::cout << ") on subentity?\n";
            }

            assert ( points.size ( ) > 0 );

            // geometrical dimension should be equal
            assert ( points[0].size ( ) == point.size ( ) );

            //  if (point.size() == 2 && points.at(0).size() == 3)
            //    point.push_back(0.);
            //  assert(point.size() == points.at(0).size());

            DataType eps = 1.0e-6;
            DataType delta = 1.0e-6;
            if ( typeid (DataType ) == typeid (double ) )
                delta = 1.0e-13;

            DataType dist = 1.0; // value to be calculated

            if ( points.size ( ) == 1 )
            {
                // single point

                DataType temp = 0.0;
                for ( size_t dim = 0, e_dim = point.size ( ); dim != e_dim; ++dim )
                    temp += ( point[dim] - points[0][dim] )
                    * ( point[dim] - points[0][dim] );

                assert ( temp >= 0.0 );
                dist = std::sqrt ( temp );
            }
            else if ( points.size ( ) == 2 )
            {
                // line

                // distance between point P and even (defined by point A and B)

                Coord AP;
                AP.push_back ( point[0] - points[0][0] );
                AP.push_back ( point[1] - points[0][1] );
                if ( point.size ( ) == 3 )
                    AP.push_back ( point[2] - points[0][2] );
                else
                    AP.push_back ( 0.0 );

                DataType tmp = AP[0] * AP[0]
                        + AP[1] * AP[1]
                        + AP[2] * AP[2];

                assert ( tmp >= 0.0 );
                DataType norm_AP = std::sqrt ( tmp );

                Coord AB;
                AB.push_back ( points[1][0] - points[0][0] );
                AB.push_back ( points[1][1] - points[0][1] );
                if ( point.size ( ) == 3 )
                    AB.push_back ( points[1][2] - points[0][2] );
                else
                    AB.push_back ( 0.0 );

                tmp = AB[0] * AB[0]
                        + AB[1] * AB[1]
                        + AB[2] * AB[2];

                assert ( tmp >= 0.0 );
                DataType norm_AB = std::sqrt ( tmp );

                DataType AP_times_AB = AP[0] * AB[0]
                        + AP[1] * AB[1]
                        + AP[2] * AB[2];

                tmp = norm_AP * norm_AP - ( AP_times_AB * AP_times_AB ) / ( norm_AB * norm_AB );

                tmp = std::abs ( tmp );
                assert ( tmp >= 0.0 || tmp <= eps );
                dist = std::sqrt ( tmp );

                if ( verbatim_mode > 0 )
                {
                    std::cout << "dist to line = " << dist << "\n";
                }

                if ( dist < delta )
                {

                    // further test: point between the two defining points of the line?

                    // set vectors
                    Vec<3, DataType> A ( points[0] );
                    Vec<3, DataType> B ( points[1] );
                    Vec<3, DataType> P ( point );

                    // Compute alpha such that P = A + alpha * AB . Then
                    // alpha == 0 => P = A
                    // alpha == 1 => P = B
                    // 0 < alpha < 1 => P lies between A and B
                    // otherwise P lies outside A and B .

                    const Vec<3, DataType> AB = B - A;
                    const Vec<3, DataType> AP = P - A;
                    const DataType alpha = dot ( AB, AP ) / dot ( AB, AB );

                    if ( alpha < -eps || alpha > 1 + eps )
                    {
                        return false;
                    }

#if 0   // old version, does not work on diagonal of triangle
                    // normed vector pointing from A to B
                    Vec<3, DataType> AB_normed ( B - A );
                    DataType norm = std::abs ( hiflow::norm ( AB_normed ) );
                    assert ( std::abs ( norm ) > 1.0e-10 );
                    AB_normed = AB_normed * ( 1. / norm );

                    if ( verbatim_mode > 0 )
                    {
                        std::cout << "A:" << A
                                << "\tB:" << B
                                << "\tP:" << P
                                << "\tAB_u:" << AB_normed << "\n";
                    }

                    // if point B is between A and B then there exists an alpha between 0 and 1
                    // such that P = A + alpha*AB_normed
                    DataType alpha = 0.;
                    bool first_time = true;
                    for ( int d = 0; d < 3; ++d )
                    {
                        if ( std::abs ( AB_normed[d] ) > 1.0e-12 )
                        {
                            //          double temp_alpha = (P[d]-A[d])/AB_normed[d];
                            DataType temp_alpha = ( P[d] - A[d] );
                            if ( std::abs ( temp_alpha ) > 1.0e-14 * AB_normed[d] )
                            {
                                if ( first_time == true )
                                {
                                    alpha = temp_alpha;
                                    first_time = false;
                                }
                                else
                                {
                                    assert ( std::abs ( alpha - temp_alpha ) < 1.0e-12 );
                                }
                            }
                        }
                    }

                    if ( verbatim_mode > 0 )
                    {
                        std::cout << "alpha = " << alpha << "\n";
                    }

                    if ( alpha >= 0. || alpha <= 1. )
                        dist = 0.;
                    if ( alpha < 0 )
                        dist = -alpha;
                    if ( alpha > 1 )
                        dist = alpha - 1.;

                    if ( verbatim_mode > 0 )
                    {
                        std::cout << "dist from A = " << dist << "\n";
                    }
#endif
                }

            }
            else
            {
                // face

                assert ( point.size ( ) == 3 );

                // 1. calculate normed normal vector of face (p0, p1, p2)

                Coord v; // v := points.at(1) - points.at(0)
                v.push_back ( points[1][0] - points[0][0] );
                v.push_back ( points[1][1] - points[0][1] );
                v.push_back ( points[1][2] - points[0][2] );

                Coord w; // w := points.at(2) - points.at(0)
                w.push_back ( points[2][0] - points[0][0] );
                w.push_back ( points[2][1] - points[0][1] );
                w.push_back ( points[2][2] - points[0][2] );

                Coord normal; // cross product
                normal.push_back ( v[1] * w[2] - v[2] * w[1] );
                normal.push_back ( v[2] * w[0] - v[0] * w[2] );
                normal.push_back ( v[0] * w[1] - v[1] * w[0] );

                // normalize normal vector
                DataType norm = sqrt ( normal[0] * normal[0]
                                       + normal[1] * normal[1]
                                       + normal[2] * normal[2] );
                normal[0] = normal[0] / norm;
                normal[1] = normal[1] / norm;
                normal[2] = normal[2] / norm;

                // 2. calculate parameter d

                DataType d = points[0][0] * normal[0]
                        + points[0][1] * normal[1]
                        + points[0][2] * normal[2];

                // 3. calculate value for considered point

                DataType d_point = point[0] * normal[0]
                        + point[1] * normal[1]
                        + point[2] * normal[2];

                // 4. calculate distance point to face

                dist = d - d_point;
                dist = std::abs ( dist );

                // 5. in case of more than 3 points

                for ( int i = 3, e_i = points.size ( ); i < e_i; ++i )
                {
                    DataType d_temp = points[i][0] * normal[0]
                            + points[i][1] * normal[1]
                            + points[i][2] * normal[2];
                    assert ( std::abs ( d - d_temp ) < eps );
                }

                // 6. for refined cells

                // until now it has been checked that the DoF point is on the plane
                // defined by the given points, now it is further checked whether
                // the DoF point lies within the convex envelope

                /// Idea of point on face test:
                /// Calculate the normal vectors by sucessive cross products
                /// \f \vec{n}_i(p_{i}-p_{i-1})\times (p-p_{i-1})\quad i=1\cdots n_points \f
                /// and check whether these normal vectors are parallel or zero.

                /// TODO (staffan): an alternative, but similar algorithm that one might try out would be as follows:
                ///
                /// For each edge compute edge vector e_i = p_i - p_{i-1} and
                /// point vector v_i = p - p_{i-1} (like above). Then compute a
                /// vector u_i = cross(e_i, v_i) and the normal to a plane
                /// containing e_i: n_i = cross(e_i, u_i) . Finally check the sign
                /// of dot(v_i, e_i) to find out what side of the plane that p
                /// lies on. This should remain the same for all edges.
                ///
                /// A special case is when p lies on the edge in question. As in
                /// the current algorithm, this would have to be tested for
                /// separately. Perhaps both algorithms actually work the same?

                if ( dist < eps )
                {
                    // normal vector
                    Vec<3, DataType> vec_normal;
                    bool first_time = true;

                    for ( int i = 1, e_i = points.size ( ); i < e_i; ++i )
                    {
                        // set vectors

                        Vec<3, DataType> vec_edge; // p_i - p_{i-1}
                        for ( int d = 0; d < 3; ++d )
                            vec_edge[d] = points[i][d] - points[i - 1][d];

                        Vec<3, DataType> vec_to_point; // point - p_{i-1}
                        for ( int d = 0; d < 3; ++d )
                            vec_to_point[d] = point[d] - points[i - 1][d];

                        // cross product

                        Vec<3, DataType> vec_normal_temp;
                        vec_normal_temp = cross ( vec_edge, vec_to_point );
                        DataType norm = std::abs ( hiflow::norm ( vec_normal_temp ) );
                        if ( verbatim_mode > 0 )
                        {
                            std::cout << "vec_edge = " << vec_edge << ", " << "\nvec_to_point = " << vec_to_point << ",\nvec_normal_temp = " << vec_normal_temp << "\n";
                            std::cout << "norm = " << norm << "\n";
                        }
                        if ( norm > delta )
                            vec_normal_temp = vec_normal_temp * ( 1. / hiflow::norm ( vec_normal_temp ) );
                        if ( norm > delta )
                        {
                            if ( first_time )
                            {
                                vec_normal = vec_normal_temp;
                                first_time = false;
                            }
                            else
                            {
                                DataType diff = 0;
                                for ( int d = 0; d < 3; ++d )
                                    diff += std::abs ( vec_normal[d] - vec_normal_temp[d] );

                                if ( diff > delta )
                                    dist = 1.0; // not within the convex envelope
                            }
                        }
                    }
                }
                else
                {
                    if ( verbatim_mode > 0 )
                    {
                        std::cout << "dist > eps: dist = " << dist << ", eps = " << eps << "\n";
                    }
                }
            }

            if ( dist < eps )
                return true;

            return false;
        }

        template class FEType<double>;
        template class FEType<float>;

    } // namespace doffem
} // namespace hiflow
