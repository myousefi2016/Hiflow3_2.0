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

#ifndef HIFLOW_GENERIC_ASSEMBLY_ALGORITHM
#    define HIFLOW_GENERIC_ASSEMBLY_ALGORITHM

#    include "fem/femanager.h"
#    include "fem/fetype.h"
#    include "mesh/cell_type.h"
#    include "mesh/entity.h"
#    include "space/element.h"
#    include "assembly.h"

/// \file generic_assembly_algorithm.h
///
/// \author Staffan Ronnas, Simon Gawlok
///
/// \details This header contains helper templates for implementing global
/// assembly algorithms in a structured way. It should only be
/// included by source files in which classes derived from
/// GlobalAssembler are derived.
///
/// The heavy use of templates in this file, as well as the different
/// assembly implementations is motivated by the wish to avoid code
/// duplication for the assembly algorithms. Without the present
/// construction, the code for each assembly procedure would look
/// quite similar. The difference between vector and matrix assembly
/// is basically only what is done with the locally assembled object
/// for each element. The scalar assembly also follows the same
/// general pattern, with some differences concerning which elements
/// are actually visited.
///
/// Assembly over the boundary and the interior of the domain follow
/// the same general pattern, but differ sufficiently to warrant two
/// separate functions. The treatment of the locally assembled objects
/// is however identical between both algorithms.
///
/// In order to avoid having to code 6 = 3 types of global objects
/// (scalar, vector, matrix) * 2 assembly algorithms (interior and
/// boundary), we split the different parts into separate classes and
/// get away with only 5 = 3 + 2. This might seem like a small gain,
/// but for future expansion it is significant. Adding a new assembly
/// algorithm (i.e. integral over interior faces) now only requires
/// writing one new function instead of three. Additionally,
/// alternative implementations (new subclasses of GlobalAssembler)
/// can in many cases use the algorithms here to shorten its code.
///
/// How is the code organized? The overall assembly algorithms are
/// implemented as template classes taking a policy class as its
/// arguments. The algorithm template defines a function assemble()
/// which take a LocalAssembler functor and a QuadraturePolicy functor
/// as template arguments. This corresponds to the two user-definable
/// entities in the GlobalAssembler interface class. The assemble()
/// function performs the global assembly, but delegates the central
/// tasks to the AssemblyPolicy class.
///
/// The AssemblyPolicy class must fulfil a range of requirements. The
/// first, and perhaps most surprising thing is that it should derive
/// from algorithm template instantiated with the AssemblyPolicy class
/// itself, so that you will see the following:
///
/// class Policy : public Algorithm<Policy> {...}
///
/// This is a standard technique in C++ template programming, known as
/// the "Curiously Recurring Template Pattern"
/// (http://en.wikipedia.org/wiki/Curiously_Recurring_Template_Pattern).
/// It makes it possible to cast the "this" pointer to a pointer to
/// the policy class, and use this for calling the required functions.
/// It would not be absolutely necessary here, but it avoids
/// having to add another parameter to the assemble() function, or a
/// member to the algorithm class.
///
/// Furthermore, the AssemblyPolicy must define two types:
///
/// - LocalObjectType is the type of the assembly target for the local
/// assembly (such as double, LocalVector or LocalMatrix). This is
/// coupled to the type of LocalAssembler that is used.
///
/// - QuadratureType is the type of the quadrature formula. This is
/// coupled to the type of QuadraturePolicy that is used.
///
/// AssemblyPolicy must also define three functions:
///
/// - has_next() returns true as long as there are more elements to
/// traverse
///
/// - next() returns a const reference to the next element to visit
/// for local assembly
///
/// - reset(LocalObjectType&) is called before the local assembly, and
/// should reset the local assembly object to zero.
///
/// - add(const Element&, const LocalObjectType&) is called after the
/// local assembly and should add the value(s) in the local object to
/// the global object.
///
/// The first three functions have been implemented in the
/// AssemblyAlgorithmBase class, which can be derived from to quickly
/// implement new AssemblyPolicy classes.
///
/// \see AssemblyAlgorithmBase

namespace hiflow
{

    /// \brief Generic assembly algorithm for global assembly on the
    /// interior of the mesh.

    template <class AssemblyPolicy, class DataType>

    class InteriorAssemblyAlgorithm
    {
      public:

        template<class LocalAssembler, class QuadraturePolicy>
        void assemble ( LocalAssembler& local_asm,
                        QuadraturePolicy& quadrature_policy )
        {
            AssemblyPolicy& ref = static_cast < AssemblyPolicy& > ( *this );

            typename AssemblyPolicy::LocalObjectType local_obj;
            typename AssemblyPolicy::QuadratureType quadrature;

            while ( ref.has_next ( ) )
            {
                const Element<DataType>& elem = ref.next ( );

                // Choose quadrature via policy.
                quadrature_policy ( elem, quadrature );

                // Assemble locally.
                ref.reset ( local_obj );
                local_asm ( elem, quadrature, local_obj );

                // Assemble globally.
                ref.add ( elem, local_obj );
            }
        }
    };

    /// \brief Generic assembly algorithm for global assembly on the
    /// interior of the mesh.

    template <class AssemblyPolicy, class DataType>
    struct BoundaryAssemblyAlgorithm
    {

        template<class LocalAssembler, class QuadraturePolicy>
        void assemble ( LocalAssembler& local_asm,
                        QuadraturePolicy& quadrature_policy )
        {

            AssemblyPolicy& ref = static_cast < AssemblyPolicy& > ( *this );

            typename AssemblyPolicy::LocalObjectType local_obj;

            while ( ref.has_next ( ) )
            {
                const Element<DataType>& elem = ref.next ( );

                if ( elem.is_boundary ( ) )
                {
                    // Get facets of cell that lie on the boundary.
                    const std::vector<int> bdy_facet_numbers
                            = elem.boundary_facet_numbers ( );

                    // Reset local object
                    ref.reset ( local_obj );

                    for ( std::vector<int>::const_iterator it = bdy_facet_numbers.begin ( ),
                          end = bdy_facet_numbers.end ( ); it != end; ++it )
                    {

                        // Create a mapped quadrature based on the facet quadrature.
                        typename AssemblyPolicy::QuadratureType mapped_quadrature;
                        quadrature_policy ( elem, mapped_quadrature, *it );

                        // Assemble locally, without resetting.
                        local_asm ( elem, *it, mapped_quadrature, local_obj );
                    }

                    // Assemble globally.
                    ref.add ( elem, local_obj );
                }
            }
        }
    };

    //using doffem::FEManager;
    // using doffem::FEType;
    using namespace doffem;

    /// \TODO do not base ordering only on first variable.

    template<class DataType>
    struct ElementCmp
    {

        ElementCmp ( const FEManager<DataType>& man ) : fe_man_ ( man )
        {
        }

        bool operator() ( int i, int j ) const
        {
            // compare based on FEType of first variable
            const FEType<DataType>* fe_i = fe_man_.get_fe_on_cell ( i, 0 );
            const FEType<DataType>* fe_j = fe_man_.get_fe_on_cell ( j, 0 );
            return fe_i->get_global_id ( ) < fe_j->get_global_id ( );
        }
        const FEManager<DataType>& fe_man_;
    };

    /// \brief Initialze traversal

    template<class DataType>
    inline void init_traversal ( const VectorSpace<DataType>& space, std::vector<int>& iteration_order )
    {
        const mesh::Mesh& mesh = space.mesh ( );
        const int num_elements = mesh.num_entities ( mesh.tdim ( ) );

        iteration_order.resize ( num_elements );
        for ( size_t i = 0; i != num_elements; ++i )
        {
            iteration_order[i] = i;
        }
    }

    /// \brief Compute order of assembly, based on FEType and degree.

    template<class DataType>
    inline void sort_elements ( const VectorSpace<DataType>& space, std::vector<int>& iteration_order )
    {

        std::sort ( iteration_order.begin ( ), iteration_order.end ( ), ElementCmp<DataType>( space.fe_manager ( ) ) );
    }

    /// \brief Basic services for the implementation of an assembly
    /// algorithm.
    ///
    /// \details Provides some common functionality, such as sorted
    /// traversal, and obtaining the dofs of the current element. This
    /// class can be used together with both interior assembly and
    /// boundary assembly. The main interface to subclasses is through
    /// protected variables.
    ///
    /// In order to support all variations of algorithm type and
    /// AssemblyPolicy classes, this base class takes both of them as
    /// template arguments (the algorithm type is a _template_
    /// template argument), and derives from AlgorithmType<AssemblyPolicy> .

    template<template<class, class> class AlgorithmType, class AssemblyPolicy, class DataType>
    class AssemblyAlgorithmBase : public AlgorithmType<AssemblyPolicy, DataType>
    {
      public:

        AssemblyAlgorithmBase ( const VectorSpace<DataType>& space )
        : space_ ( space ), curr_ ( 0 ), elem_ ( space, 0 )
        {

            //initialize traversal_
            init_traversal ( space, this->traversal_ );
        }

        bool has_next ( )
        {
            return this->curr_ < this->traversal_.size ( );
        }

        void reset ( DataType& local_val )
        {
            local_val = 0.;
        }

        void reset ( typename GlobalAssembler<DataType>::LocalVector& local_vec )
        {
            local_vec.clear ( );
            local_vec.resize ( this->dof_.size ( ), 0. );
        }

        void reset ( typename GlobalAssembler<DataType>::LocalMatrix& local_mat )
        {
            local_mat.Resize ( this->dof_.size ( ), this->dof_.size ( ) );
            local_mat.Zeros ( );
        }

        const Element<DataType>& next ( )
        {
            assert ( this->has_next ( ) );

            this->elem_ = Element<DataType>( this->space_, this->traversal_[this->curr_] );
            this->elem_.get_dof_indices ( this->dof_ );

            ++( this->curr_ );
            return this->elem_;
        }

      protected:
        const VectorSpace<DataType>& space_;
        std::vector<int> traversal_;
        std::vector<int> dof_;
        int curr_;
        Element<DataType> elem_;

    };
}

#endif
