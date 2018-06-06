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

#ifndef HIFLOW_ADAPTIVITY_TIME_MESH
#    define HIFLOW_ADAPTIVITY_TIME_MESH

/// \author Philipp Gerstner

#    include <string>
#    include <vector>
#    include "common/sorted_array.h"

namespace hiflow
{
    /// \class TimeMesh time_mesh.h
    /// \brief Sequential, 1D mesh used for time discretization
    ///

    template<class DataType>
    class TimeMesh
    {
      public:
        TimeMesh ( );

        ~TimeMesh ( )
        {
            this->clear ( );
        }

        /// \brief reset internal variables to default values
        virtual void clear ( );

        /// \brief add a set of new points in time to mesh
        /// @param[in] times set of time points to add
        /// @param[in] history_index history index of mesh to which points should be added
        virtual void add_times ( std::vector<DataType>& times, int history_index );

        /// \brief add a set of new points in time to mesh
        /// @param[in] time time point to add
        /// @param[in] history_index history index of mesh to which points should be added
        virtual void add_time ( DataType time, int history_index );

        /// \brief change time point 
        /// @param[in] time_step index of time point to be changed
        /// @param[in] add_time new time point
        /// @param[in] history_index history index of mesh to which points should be added
        virtual void set_time ( int time_step, DataType add_time, int history_index );

        /// \brief get first time point of mesh
        /// @param[in] history_index index of considered mesh
        /// @return start point
        virtual DataType start ( int history_index ) const;

        /// \brief get last time point of mesh
        /// @param[in] history_index index of considered mesh
        /// @return end point
        virtual DataType end ( int history_index ) const;

        /// \brief get first time point of most recent mesh
        /// @return start point
        virtual DataType start ( ) const;

        /// \brief get last time point of most recent mesh
        /// @return end point        
        virtual DataType end ( ) const;

        /// \brief get number of intervals in considered mesh
        /// @param[in] history_index index of considered mesh
        /// @return number of intervals 
        virtual int num_intervals ( int history_index ) const;

        /// \brief get number of intervals in most recent mesh
        /// @return number of intervals         
        virtual int num_intervals ( ) const;

        /// \brief get number meshes in hierarchy
        /// @return number of meshes
        virtual int num_levels ( ) const;

        /// \brief get time step size of considered mesh
        /// @param[in] time_step index of subinterval
        /// @param[in] history_index index of considered mesh
        /// @return length of sub interval
        virtual DataType delta_t ( int time_step, int history_index ) const;

        /// \brief get specific time point
        /// @param[in] time_step index of time point
        /// @param[in] history_index index of considered mesh
        /// @return time point
        virtual DataType time ( int time_step, int history_index ) const;

        /// \brief get specific time point in most recent mesh
        /// @param[in] time_step index of time point
        /// @return time point
        virtual DataType time ( int time_step ) const;

        /// \brief get all time points in specific mesh
        /// @param[in] history_index index of considered mesh
        /// @return time points
        virtual std::vector<DataType> get_all_times ( int history_index ) const;

        /// \brief adapt specific mesh
        /// @param[in] markers vector of refinement markers for each sub interval in specified mesh \br
        /// 1: refine, 0: keep, <0: coarsen
        /// @param[in] history_index index of considered mesh
        virtual void refine ( const std::vector<int>& markers, int history_index );

        /// \brief load complete hierarchy of time meshes
        /// @param[in] filename path to saved meshed
        /// \param[in] num_mesh number of meshes to be loaded
        virtual void load_all ( const std::string& filename, int num_mesh );

        /// \brief load specific mesh
        /// @param[in] filename path to saved meshed
        /// @param[in] history_index index of considered mesh
        virtual void load ( const std::string& filename, int history_index );

        /// \brief save complete hierarchy of meshes
        /// @param[in] filename path to saved meshed
        virtual void save_all ( const std::string& filename ) const;

        /// \brief save specific mesh
        /// @param[in] filename path to saved meshed
        /// @param[in] history_index index of considered mesh
        virtual void save ( const std::string& filename, int history_index ) const;

        /// \brief ensure to have an even number of sub intervals in mesh
        /// @param[in] history_index index of considered mesh
        virtual void make_regular ( int history_index );

      protected:

        /// data structure that stores time points
        std::vector< SortedArray<DataType> > times_;

    };
}

#endif
