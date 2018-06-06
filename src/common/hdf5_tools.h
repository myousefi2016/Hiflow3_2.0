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

#ifndef HIFLOW_COMMON_HDF5_TOOLS_H
#    define HIFLOW_COMMON_HDF5_TOOLS_H

#    include "config.h"

#    ifdef WITH_HDF5

#        include "hdf5.h"

#        include <algorithm>
#        include <cstddef>
#        include <string>
#        include <tr1/memory>
#        include <map>
#        include <sstream>

/// \brief Classes and functions which help reading and writing data from/to HDF5 files.
/// \author Teresa Beck, Simon Gawlok

namespace hiflow
{

    /// \brief Maximum length of std::string objects that
    /// can be read/written from/to HDF5 files.
    const unsigned int MaxStrLength = 512;

    ///
    /// \brief Structure for handling strings.
    ///
    /// HDF5 does not support writing strings. To do so, they must be
    /// converted into character arrays, called StringContainer, first.
    /// Maximum string size is 512 Byte.

    struct StringContainer
    {
        char string[MaxStrLength];
    };

    /// \brief Function to get length of a string.
    ///
    /// \param input String
    ///
    /// \returns Size of string

    inline int get_stringlength ( const std::string &input )
    {
        return input.size ( );
    }

    /// \brief Function to get maximum length of a vector of strings.
    ///
    /// \param input Vector of strings
    ///
    /// \returns Maximum size of vector of strings.

    inline int get_max_stringlength ( const std::vector<std::string> &vi )
    {
        int maxsize = 0;
        for ( size_t i = 0; i != vi.size ( ); ++i )
        {
            if ( get_stringlength ( vi[i] ) > maxsize )
            {
                maxsize = get_stringlength ( vi[i] );
            }
        }
        return maxsize;
    }

    ///
    /// \brief Class for HDF5 property lists.
    ///
    /// Property lists are a mechanism to modify the default
    /// behavior when creating or accessing objects such as files,
    /// groups or datasets.
    /// The default property list is specified by H5P_DEFAULT and handels
    /// the most common needs. Thus in some cases, such as accessing a file,
    /// transferring data or creating an attribute, more specific
    /// property lists are required which can be created with H5Propertylist.

    class H5Propertylist
    {
      public:
        /// \brief Constructor

        H5Propertylist ( ) : id_ ( -1 )
        {
        }

        /// \brief Constructor
        ///
        /// \param purpose HDF5 Purpose of property list
        ///
        /// Exemplary purposes are H5P_DATASET_CREATE (create a dataset),
        /// H5P_DATASET_ACCESS (access data), H5P_DATASET_XFER (transfer data),...

        H5Propertylist ( hid_t purpose )
        {
            id_ = H5Pcreate ( purpose );
        }

        H5Propertylist ( const H5Propertylist& type )
        {
            id_ = H5Pcopy ( type.id_ );
        }

        H5Propertylist& operator= ( const H5Propertylist& type )
        {
            if ( this != &type )
            {
                if ( id_ > 0 )
                {
                    H5Pclose ( id_ );
                }
                id_ = H5Pcopy ( type.id_ );
            }
            return *this;
        }

        /// \brief Destructor

        ~H5Propertylist ( )
        {
            if ( id_ > 0 )
            {
                H5Pclose ( id_ );
            }
        }

        /// \brief Access the identifier of a property list.
        ///
        /// \return The identifier

        const hid_t& id ( ) const
        {
            return id_;
        }

      private:
        hid_t id_;
    };

    ///
    /// \brief Class for HDF5 files
    ///
    /// An HDF5 file is a binary file that contains data and supporting metadata.
    /// To create or access a HDF5 file, one must specify a file name and the
    /// mode (read or write).

    class H5File
    {
      public:
        /// \brief Constructor
        ///
        /// \param filename Name of file
        /// \param mode Read or write
        ///
        /// If the file already exists, it will be opened.
        /// If it does not and in case of using the write mode, a
        /// new file will be created. Trying to read from an inexistent
        /// file will return an error.

        H5File ( const std::string& filename, const std::string& mode, const MPI_Comm& comm )
        {
            H5Propertylist faplist ( H5P_FILE_ACCESS );
            H5Pset_fapl_mpio ( faplist.id ( ), comm, MPI_INFO_NULL );

            FILE *file = fopen ( filename.c_str ( ), "r" );
            if ( file == NULL )
            {
                if ( mode == "w" ) id_ = H5Fcreate ( filename.c_str ( ),
                                                     H5F_ACC_TRUNC, H5P_DEFAULT, faplist.id ( ) );
                else std::cerr << "HDF5 ERROR:  H5File does not exist! \n";
            }
            else
            {
                id_ = H5Fopen ( filename.c_str ( ), H5F_ACC_RDWR, faplist.id ( ) );
                fclose ( file );
            }
        }

        /// \brief Destructor

        ~H5File ( )
        {
            H5Fclose ( id_ );
        }
        /// \brief Access the identifier of a file
        ///
        /// \return The identifier

        const hid_t& id ( ) const
        {
            return id_;
        }

        hid_t id_;

    };

    typedef std::tr1::shared_ptr<H5File> H5FilePtr;

    ///
    /// \brief Class for HDF5 groups
    ///
    /// An HDF5 group is a structure containing zero or more HDF5 objects.
    /// An object can either be a dataset or a group.
    /// To create or access a HDF5 group, one must obtain the location
    /// identifier where to create the group. This can be the identifier of
    /// a file or the identifier of a group. Furthermore, the name of group
    /// and the mode (read or write) must be specified.

    class H5Group
    {
      public:
        /// \brief Constructor
        ///
        /// \param loc_id Identifier of place where to create the group (can be a file identifier or a group identifier)
        /// \param groupname Name of group
        /// \param mode Read or write
        ///
        /// In write mode a group will be opened if it already exists
        /// and created if it does not. Accessing an inexistent
        /// group in read mode will return an error.

        H5Group ( H5FilePtr file_ptr,
                  const std::string& groupname,
                  const std::string& mode )
        {
            H5Propertylist lapl_id ( H5P_LINK_ACCESS );
            std::string newpath, checkpath;
            size_t first, last;
            hid_t temp_id;
            htri_t link = 0;
            // Init pathname and location identifier with groupname and file identifier
            file_ptr_ = file_ptr;
            newpath = groupname.c_str ( );
            temp_id = file_ptr->id ( );

            // Open root group
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
            id_ = H5Gopen ( temp_id, "/" );
#        else
            id_ = H5Gopen1 ( temp_id, "/" );
#        endif
            temp_id = id_;

            // assimilate strings
            if ( newpath.find_first_of ( "/" ) == 0 ) newpath.erase ( 0, 1 );
            if ( newpath.find_last_of ( "/" ) != newpath.size ( ) - 1 ) newpath.append ( "/" );
            if ( newpath == "" ) newpath = "/";

            // Check if further links to root group exist recursively
            if ( newpath != "/" )
            {
                do
                {
                    first = newpath.find_first_of ( "/" );
                    last = newpath.find_last_of ( "/" );
                    checkpath = newpath.substr ( 0, first );

                    link = H5Lexists ( temp_id, checkpath.c_str ( ), lapl_id.id ( ) );
                    if ( link <= 0 )
                    {
                        if ( mode == "w" )
                        {
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
                            id_ = H5Gcreate ( temp_id, checkpath.c_str ( ), 0 );
#        else
                            id_ = H5Gcreate1 ( temp_id, checkpath.c_str ( ), 0 );
#        endif
                        }
                        else
                        {
                            std::cerr << "HDF5 ERROR:  H5Group does not exist! \n";
                        }
                    }
                    else
                    {
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
                        id_ = H5Gopen ( temp_id, checkpath.c_str ( ) );
#        else
                        id_ = H5Gopen1 ( temp_id, checkpath.c_str ( ) );
#        endif
                    }

                    H5Gclose ( temp_id );

                    temp_id = id_;
                    newpath = newpath.substr ( first + 1, groupname.size ( ) );
                }
                while ( first != last );
            }
        }

        /// \brief Destructor

        ~H5Group ( )
        {
            H5Gclose ( id_ );
        }

        /// \brief Access the identifier of a group.
        ///
        /// \return The identifier

        const hid_t& id ( ) const
        {
            return id_;
        }

        hid_t id_;
        H5FilePtr file_ptr_;

    };

    typedef std::tr1::shared_ptr<H5Group> H5GroupPtr;

    ///
    /// \brief Class for HDF5 Dataspace
    ///
    /// A dataspace describes the dimensionality of a dataset.
    /// Currently, it is a regular one-dimensional array of N data points.

    class H5Dataspace
    {
      public:
        /// \brief Constructor

        H5Dataspace ( ) : id_ ( -1 )
        {
        }

        /// \brief Constructor
        ///
        /// \param dims Dimension of one-dimensional dataspace

        H5Dataspace ( hsize_t dims )
        {
            hsize_t maxdims[1] = { H5S_UNLIMITED };

            id_ = H5Screate_simple ( 1, &dims, maxdims );
        }

        H5Dataspace ( hid_t dataset_id )
        {
            id_ = H5Dget_space ( dataset_id );
        }

        H5Dataspace ( const H5Dataspace& type )
        {
            id_ = H5Scopy ( type.id_ );
        }

        H5Dataspace& operator= ( const H5Dataspace& type )
        {
            if ( this != &type )
            {
                if ( id_ > 0 )
                {
                    H5Sclose ( id_ );
                }
                id_ = H5Scopy ( type.id_ );
            }
            return *this;
        }

        /// \brief Destructor

        ~H5Dataspace ( )
        {
            if ( id_ > 0 )
            {
                H5Sclose ( id_ );
            }
        }

        /// \brief Get number of elements in a dataspace.
        ///
        /// \return Number of elements in a dataspace

        hsize_t number_of_elements ( )
        {
            return H5Sget_simple_extent_npoints ( id_ );
        }

        /// \brief Access the identifier of a file.
        ///
        /// \return The identifier

        const hid_t& id ( ) const
        {
            return id_;
        }

      private:
        hid_t id_;
    };

    ///
    /// \brief Class for HDF5 datatypes
    ///

    class H5Datatype
    {
      public:
        /// \brief Constructor

        H5Datatype ( ) : id_ ( -1 )
        {
        }

        H5Datatype ( hid_t id )
        {
            id_ = H5Tcopy ( id );
        }

        H5Datatype ( const H5Datatype& type )
        {
            id_ = H5Tcopy ( type.id_ );
        }

        H5Datatype& operator= ( const H5Datatype& type )
        {
            if ( this != &type )
            {
                if ( id_ > 0 )
                {
                    H5Tclose ( id_ );
                }
                id_ = H5Tcopy ( type.id_ );
            }
            return *this;
        }

        /// Destructor

        ~H5Datatype ( )
        {
            if ( id_ > 0 )
            {
                H5Tclose ( id_ );
            }
        }

        /// \brief Access the identifier of a datatype.
        ///
        /// \return The identifier

        const hid_t& id ( ) const
        {
            return id_;
        }

      private:
        hid_t id_;
    };

    /// \brief Get HDF5 datatype
    ///
    /// \return HDF5 datatype
    ///
    /// Supported datatypes are integer, long integer,
    /// float, double, character and std::string.

    template<class T>
    inline H5Datatype get_hdf5_data_type ( const T& )
    {
        assert ( false );
        throw "Unknown type!";
        return H5Datatype ( H5T_NATIVE_INT );
    }

    template<>
    inline H5Datatype get_hdf5_data_type ( const int& )
    {
        return H5Datatype ( H5T_NATIVE_INT );
    }

    template<>
    inline H5Datatype get_hdf5_data_type ( const long int& )
    {
        return H5Datatype ( H5T_NATIVE_LONG );
    }

    template<>
    inline H5Datatype get_hdf5_data_type ( const double& )
    {
        return H5Datatype ( H5T_NATIVE_DOUBLE );
    }

    template<>
    inline H5Datatype get_hdf5_data_type ( const float& )
    {
        return H5Datatype ( H5T_NATIVE_FLOAT );
    }

    template<>
    inline H5Datatype get_hdf5_data_type ( const char& )
    {
        return H5Datatype ( H5T_C_S1 );
    }

    template<typename VecT>
    H5Datatype get_hdf5_data_type ( const std::vector<VecT>& data )
    {
        return get_hdf5_data_type ( data[0] );
    }

    /// \brief Create compound HDF5 datatype for vector of strings.
    ///
    /// \param data Input vector of strings
    ///
    /// \return Compound string datatype
    ///
    /// The function creates a compound datatype of type H5T_STRING.
    /// Size of Datatype H5T_STRING is length of data. If data is empty
    /// (when reading from a file into an empty buffer) size is variable.

    template<>
    inline H5Datatype get_hdf5_data_type ( const std::vector<std::string>& data )
    {
        int max = get_max_stringlength ( data );
        H5Datatype strtype = H5Tcopy ( H5T_C_S1 );
        if ( max <= 0 ) H5Tset_size ( strtype.id ( ), H5T_VARIABLE );
        else H5Tset_size ( strtype.id ( ), max );
        H5Datatype datatype = H5Tcreate ( H5T_COMPOUND, MaxStrLength );
        H5Tinsert ( datatype.id ( ), "string", 0, strtype.id ( ) );
        return datatype;
    }

    /// \brief Create standard HDF5 datatype for single string.
    ///
    /// \param data Input string
    ///
    /// \return Standard string datatype
    ///
    /// Size of Datatype H5T_STRING is length of data. If data is empty
    /// (when reading from a file into an empty buffer) size is variable.

    template<>
    inline H5Datatype get_hdf5_data_type ( const std::string& data )
    {
        int max = get_stringlength ( data );
        H5Datatype datatype = H5Tcopy ( H5T_C_S1 );
        if ( max <= 0 ) H5Tset_size ( datatype.id ( ), H5T_VARIABLE );
        else H5Tset_size ( datatype.id ( ), max );
        return datatype;
    }

    template<>
    inline H5Datatype get_hdf5_data_type ( const StringContainer& data )
    {
        return get_hdf5_data_type ( std::string ( data.string ) );
    }

    /// \brief Find HDF5 datatype, only used in context with reading/writing maps.
    ///
    /// \param data Vector of pairs
    ///
    /// \return Compound datatype of input pair

    template<typename KeyT, typename ValueT>
    H5Datatype get_hdf5_data_type ( const std::vector< std::pair<KeyT, ValueT> > &data )
    {
        hsize_t size = data.size ( );

        std::vector<KeyT> vfirst;
        std::vector<ValueT> vsecond;
        vfirst.resize ( size );
        vsecond.resize ( size );

        for ( int i = 0; i < size; i++ )
        {
            vfirst[i] = ( data[i] ).first;
            vsecond[i] = ( data[i] ).second;
        }

        H5Datatype keytype = get_hdf5_data_type ( vfirst );
        H5Datatype valuetype = get_hdf5_data_type ( vsecond );

        int key_offset = 0;
        int value_offset = H5Tget_size ( keytype.id ( ) );
        int pair_size = H5Tget_size ( keytype.id ( ) ) + H5Tget_size ( valuetype.id ( ) );

        H5Datatype datatype = H5Tcreate ( H5T_COMPOUND, pair_size );
        H5Tinsert ( datatype.id ( ), "Key", key_offset, keytype.id ( ) );
        H5Tinsert ( datatype.id ( ), "Value", value_offset, valuetype.id ( ) );
        return datatype;
    }

    template<typename KeyT, typename ValueT>
    H5Datatype get_hdf5_data_type ( const std::pair<KeyT, ValueT> &data )
    {
        H5Datatype keytype = get_hdf5_data_type ( data.first );
        H5Datatype valuetype = get_hdf5_data_type ( data.second );

        int key_offset = 0;
        int value_offset = H5Tget_size ( keytype.id ( ) );
        int pair_size = H5Tget_size ( keytype.id ( ) ) + H5Tget_size ( valuetype.id ( ) );

        H5Datatype datatype = H5Tcreate ( H5T_COMPOUND, pair_size );
        H5Tinsert ( datatype.id ( ), "Key", key_offset, keytype.id ( ) );
        H5Tinsert ( datatype.id ( ), "Value", value_offset, valuetype.id ( ) );
        return datatype;
    }

    ///
    /// \brief Class for HDF5 datasets
    ///
    /// An HDF5 dataset is an n-dimensional array of data elements
    /// (currently a one-dimensional array).
    /// To create or access a dataset, one must first obtain
    /// the pointer to a HDF5 group where it is or is being linked to.
    /// Furthermore, one must specify the size of the dimension of the
    /// dataset, its name, the mode (read or write) and a buffer with
    /// data to be written to a file or a buffer for data read from a file.

    class H5Dataset
    {
      public:
        /// \brief Constructor
        ///
        /// \param group_ptr Pointer to HDF5 group the dataset is linked to
        /// \param dim_space Size of dimensions of dataset
        /// \param datasetname Name of dataset
        /// \param mode Read or write
        /// \param data Buffer with data to be written to the file or read from file
        ///
        /// If the dataset already exists, it will be opened. In case of
        /// write mode, data will be overwritten. If it does not exist a
        /// new dataset will be created in write mode. Trying to read from
        /// an inexistent dataset will return an error.

        template<class DataType>
        H5Dataset ( H5GroupPtr group_ptr,
                    hsize_t dim_space,
                    const std::string& datasetname,
                    const std::string& mode,
                    DataType* data )
        {
            group_ptr_ = group_ptr;
            h5dtype = get_hdf5_data_type ( *data );
            // Open Dataset
            // Check if Dataset already exists.
            H5Propertylist lapl ( H5P_DATASET_ACCESS );
            htri_t dataset = H5Lexists ( group_ptr->id ( ),
                                         datasetname.c_str ( ), lapl.id ( ) );

            // Create Dataspace for Dataset.
            space = H5Dataspace ( dim_space );

            // Dataset storage layout is chunked.
            // This is necessary in order to extent Datasets once they are created.
            H5Propertylist chpl ( H5P_DATASET_CREATE );
            H5Pset_chunk ( chpl.id ( ), 1, &dim_space );

            // If Dataset exists, overwrite old Dataset.
            // If Dataset does not exist, create new Dataset.
            if ( mode == "r" )
            {
                if ( dataset == true )
                {
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
                    id_ = H5Dopen ( group_ptr->id ( ), datasetname.c_str ( ) );
#        else
                    id_ = H5Dopen1 ( group_ptr->id ( ), datasetname.c_str ( ) );
#        endif
                }
                else std::cerr << "HDF5 ERROR:  H5Dataset does not exist! \n";
            }
            else
            {
                if ( dataset == true )
                {
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
                    id_ = H5Dopen ( group_ptr->id ( ), datasetname.c_str ( ) );
#        else
                    id_ = H5Dopen1 ( group_ptr->id ( ), datasetname.c_str ( ) );
#        endif

                    H5Dataspace filespace = H5Dataspace ( id_ );
                    H5Datatype filetype = H5Dget_type ( id_ );

                    // If Datatypes don't match, old Dataset will be unlinked from the Group and
                    // new Dataset with new Datatype is linked to Group instead.
                    // Memory of old Dataset can not be released.
                    if ( H5Tequal ( h5dtype.id ( ), filetype.id ( ) ) <= 0 )
                    {
                        std::cerr << "HDF5 WARNING: H5Dataset already exists with a different Datatype! \n";
                        std::cerr << "        The existing H5Dataset will be overwritten. \n";
                        H5Ldelete ( group_ptr->id ( ), datasetname.c_str ( ), lapl.id ( ) );
                        H5Dclose ( id_ );
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
                        id_ = H5Dcreate ( group_ptr->id ( ), datasetname.c_str ( ), h5dtype.id ( ),
                                          space.id ( ), chpl.id ( ) );
#        else
                        id_ = H5Dcreate1 ( group_ptr->id ( ), datasetname.c_str ( ), h5dtype.id ( ),
                                           space.id ( ), chpl.id ( ) );
#        endif

                    }
                    else if ( H5Tequal ( h5dtype.id ( ), filetype.id ( ) ) > 0 )
                    {
                        // If Dataspaces in File is smaller than new Dataspace, Dataset is extended.
                        if ( filespace.number_of_elements ( ) < space.number_of_elements ( ) )
                        {
                            std::cerr << "HDF5 WARNING: H5Dataset already exists with a different Dataspace! \n";
                            std::cerr << "        The existing H5Dataset will be overwritten. \n";
                            H5Dset_extent ( id_, &dim_space );
                        } // If Dataspaces in File is bigger than new Dataspace, old Dataset will be
                            // unlinked from the Group and new Dataset with new Dataspace is linked
                            // to Group instead. Memory of old Dataset can not be released.
                        else if ( filespace.number_of_elements ( ) > space.number_of_elements ( ) )
                        {
                            std::cerr << "HDF5 WARNING: H5Dataset already exists with a different Dataspace! \n";
                            std::cerr << "        The existing H5Dataset will be overwritten. \n";
                            H5Ldelete ( group_ptr->id ( ), datasetname.c_str ( ), lapl.id ( ) );
                            H5Dclose ( id_ );
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
                            id_ = H5Dcreate ( group_ptr->id ( ), datasetname.c_str ( ), h5dtype.id ( ),
                                              space.id ( ), chpl.id ( ) );
#        else
                            id_ = H5Dcreate1 ( group_ptr->id ( ), datasetname.c_str ( ), h5dtype.id ( ),
                                               space.id ( ), chpl.id ( ) );
#        endif
                        }
                    }

                }
                else
                {
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
                    id_ = H5Dcreate ( group_ptr->id ( ), datasetname.c_str ( ), h5dtype.id ( ),
                                      space.id ( ), chpl.id ( ) );
#        else
                    id_ = H5Dcreate1 ( group_ptr->id ( ), datasetname.c_str ( ), h5dtype.id ( ),
                                       space.id ( ), chpl.id ( ) );
#        endif
                }
            }

        }

        /// \brief Simplified constructor for appending attributes.
        ///
        /// \param group_ptr Pointer to HDF5 group
        /// \param datasetname Name of dataset
        ///
        /// In case of only appending an attribute, one can use this simplified
        /// dataset constructor.

        H5Dataset ( H5GroupPtr group_ptr, const std::string& datasetname )
        {
            group_ptr_ = group_ptr;
            H5Propertylist lapl ( H5P_DATASET_ACCESS );
            htri_t dataset = H5Lexists ( group_ptr->id ( ),
                                         datasetname.c_str ( ), lapl.id ( ) );
            if ( dataset == true )
            {
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
                id_ = H5Dopen ( group_ptr->id ( ), datasetname.c_str ( ) );
#        else
                id_ = H5Dopen1 ( group_ptr->id ( ), datasetname.c_str ( ) );
#        endif
            }
            else std::cerr << "HDF5 ERROR:  H5Dataset does not exist! \n";
        }

        /// \brief Destructor

        ~H5Dataset ( )
        {
            H5Dclose ( id_ );
        }

        /// \brief Access the identifier of a dataset.
        ///
        /// \return The identifier

        const hid_t& id ( ) const
        {
            return id_;
        }

        /// \brief Write a dataset from buffer into a file.
        ///
        /// \param dim_memspace Size of dimensions of dataspace describing the dimensionality of dataset
        /// \param offset Offset of processors' ownership of data to be written to the file
        /// \param data Buffer with data to be written to the file

        template<class DataType>
        void write ( hsize_t dim_memspace, hsize_t offset, DataType* data )
        {
            // Create Dataspace in memory.
            H5Dataspace memspace ( dim_memspace );

            // Select Hyperslab in the File.
            H5Sselect_hyperslab ( space.id ( ), H5S_SELECT_SET,
                                  &offset, NULL, &dim_memspace, NULL );

            // Create property list for collective dataset write.
            H5Propertylist plist ( H5P_DATASET_XFER );
            H5Pset_dxpl_mpio ( plist.id ( ), H5FD_MPIO_INDEPENDENT );

            // Write Data
            H5Dwrite ( id_, h5dtype.id ( ), memspace.id ( ),
                       space.id ( ), plist.id ( ), data );
        }

        /// \brief Read a dataset from a file into a buffer.
        ///
        /// \param dim_memspace Size of dimensions of dataspace describing the dimensionality of dataset
        /// \param offset Offset of processors' ownership of data to be written to the file
        /// \param data Buffer to receive data read from file

        template<class DataType>
        void read ( hsize_t dim_memspace, hsize_t offset, DataType* buffer )
        {
            // Create Dataspace in memory.
            H5Dataspace memspace ( dim_memspace );

            // Get Datatype of memory.
            H5Datatype memtype = get_hdf5_data_type ( *buffer );

            // Select Hyperslab in the File.
            H5Sselect_hyperslab ( space.id ( ), H5S_SELECT_SET,
                                  &offset, NULL, &dim_memspace, NULL );
            // Read Data from File.
            H5Dread ( id_, memtype.id ( ), memspace.id ( ),
                      space.id ( ), H5P_DEFAULT, buffer );
        }

        hid_t id_;
        H5Dataspace space;
        H5GroupPtr group_ptr_;
        H5Datatype h5dtype;

    };

    typedef std::tr1::shared_ptr<H5Dataset> H5DatasetPtr;

    ///
    /// \brief Class for HDF5 attributes
    ///
    /// An HDF5 attribute is a small dataset attached to an object
    /// (group or dataset). It can be used to describe the object
    /// and/or the objects intended usage.
    /// To create or access an attribut, one must first obtain a pointer to a
    /// dataset or a group, where the attribute is attached or being attached to.
    /// Furthermore, the size of its dimensions, name and the mode (read or write)
    /// must be specified.

    template<class DataType>
    class H5Attribute
    {
      public:
        /// \brief Constructor

        H5Attribute ( ) : id_ ( -1 )
        {
        }

        /// \brief Constructor for appending an attribute to a dataset.
        ///
        /// \param dataset_ptr Pointer to HDF5 dataset where to append the attribute or where the attribute is appended to
        /// \param size Size of dimensions of attribute
        /// \param att_name Name of attribute
        /// \param mode Read or write
        ///
        /// If the attribute already exists, it will be overwritten.
        /// Otherwise it will be created.

        H5Attribute ( H5DatasetPtr dataset_ptr,
                      hsize_t& size,
                      const std::string& att_name,
                      const std::string& mode )
        {
            if ( mode == "r" )
            {
                open ( dataset_ptr->id ( ), att_name );
                size = get_att_size ( );
            }
            else
            {
                create ( dataset_ptr->id ( ), size, att_name );
            }
        }

        /// \brief Constructor for appending an attribute to a group.
        ///
        /// \param dataset_ptr Pointer to HDF5 group the attribute is appended to
        /// \param size Size of dimensions of attribute
        /// \param att_name Name of attribute
        /// \param mode Read or write
        ///
        /// If the attribute already exists, it will be overwritten.
        /// Otherwise it will be created.

        H5Attribute ( H5GroupPtr group_ptr,
                      hsize_t& size,
                      const std::string& att_name,
                      const std::string& mode )
        {
            if ( mode == "r" )
            {
                open ( group_ptr->id ( ), att_name );
                size = get_att_size ( );
            }
            else
            {
                create ( group_ptr->id ( ), size, att_name );
            }
        }

        /// \brief Destructor

        ~H5Attribute ( )
        {
            if ( id_ > 0 )
            {
                H5Aclose ( id_ );
            }
        }

        /// \brief Access the identifier of an attribute.
        ///
        /// \return The id

        const hid_t& id ( ) const
        {
            return id_;
        }

        /// \brief Write an attribute into a file.
        ///
        /// \param data Attribute to be written

        void write ( DataType* data )
        {
            // Get HDF5 Datatype
            H5Datatype h5dtype;
            h5dtype = get_hdf5_data_type ( DataType ( ) );

            // Write Attribute
            H5Awrite ( id_, h5dtype.id ( ), data );
        }

        /// \brief Read an attribute from a file.
        ///
        /// \param buffer Buffer to receive data of attribute read from file

        void read ( DataType* buffer )
        {
            // Read Attribute
            H5Datatype h5dtype = get_hdf5_data_type ( buffer[0] );

            // Read Attribute from file
            // Special treatment for strings
            if ( ( H5Tget_class ( h5dtype.id ( ) ) == H5T_STRING )
                 && ( H5Tequal ( h5dtype.id ( ), H5T_C_S1 ) <= 0 ) )
            {
                // This test only passes if DataType == std::string, and not char.

                hsize_t space = get_att_size ( );
                if ( space == 0 )
                {
                    return;
                }
                std::vector<char*> string_att ( space );

                H5Aread ( id_, h5dtype.id ( ), &string_att[0] );

                for ( int i = 0; i < space; i++ )
                {
                    ( std::string& )buffer[i] = string_att[i];

                    free ( string_att[i] );
                }
            } // Default treatment for all other datatypes
            else
                H5Aread ( id_, h5dtype.id ( ), buffer );

        }

        /// \brief Get number of elements in an attribute.
        ///
        /// \return Number of elements

        hsize_t get_att_size ( )
        {
            hid_t space_id = H5Aget_space ( id_ );
            return (H5Sget_simple_extent_npoints ( space_id ) );
        }

      private:
        hid_t id_;

        H5Attribute ( const H5Attribute& type )
        {
        }

        H5Attribute& operator= ( const H5Attribute& type )
        {
        }
        /// \brief Create an attribute.
        ///
        /// \param pointer_id Pointer to HDF5 dataset or group where to append the attribute
        /// \param size Size of dimensions of attribute
        /// \param att_name Name of attribute

        void create ( hid_t pointer_id, hsize_t size,
                      const std::string& att_name )
        {
            // Create Attribute property list
            H5Propertylist atplist ( H5P_ATTRIBUTE_CREATE );

            // Check if Attribute exists
            htri_t attr = H5Aexists ( pointer_id, att_name.c_str ( ) );

            // If Attribute exists, delete old Attribute.
            if ( attr == true )
            {
                H5Adelete ( pointer_id, att_name.c_str ( ) );
            }

            // Create Dataspace for Attribute
            H5Dataspace attspace ( size );

            // Get HDF5 Datatype
            H5Datatype h5dtype;
            h5dtype = get_hdf5_data_type ( DataType ( ) );

            // Create Attribute
#        if H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8
            id_ = H5Acreate ( pointer_id, att_name.c_str ( ),
                              h5dtype.id ( ), attspace.id ( ), atplist.id ( ) );
#        else
            id_ = H5Acreate1 ( pointer_id, att_name.c_str ( ),
                               h5dtype.id ( ), attspace.id ( ), atplist.id ( ) );
#        endif
        }

        /// \brief Open an attribute.
        ///
        /// \param pointer_id Pointer to HDF5 dataset or group the attribute is appended to
        /// \param att_name Name of attribute

        void open ( hid_t pointer_id, const std::string& att_name )
        {
            // Create Attribute property list
            H5Propertylist atplist ( H5P_ATTRIBUTE_CREATE );

            // Check if Attribute exists
            htri_t attr = H5Aexists ( pointer_id, att_name.c_str ( ) );

            // If Attribute does not exist tell user
            if ( attr == false )
            {
                std::cerr << "HDF5 ERROR:  H5Attribute does not exist! \n";
            }

            // Open Attribute
            id_ = H5Aopen ( pointer_id, att_name.c_str ( ), atplist.id ( ) );
        }

    };

    /// \brief Return number of dimensions of dataset.
    ///
    /// \param group_ptr Pointer to HDF5 group
    /// \param datasetname Name of dataset
    ///
    /// \returns Number of dimensions od dataset

    inline int get_ndims_dataset ( H5GroupPtr group_ptr,
                                   const std::string& datasetname )
    {
        // Get Dataset pointer
        H5DatasetPtr dataset_ptr;
        dataset_ptr.reset ( new H5Dataset ( group_ptr, datasetname ) );

        // Get Dataspace of Dataset
        H5Dataspace dataspace ( dataset_ptr->id ( ) );

        // Get number of dimensions of Dataspace
        hsize_t ndims = H5Sget_simple_extent_ndims ( dataspace.id ( ) );
        // Return
        return int(ndims );
    }

    /// \brief Return size of dimensions of dataset.
    ///
    /// \param group_ptr Pointer to HDF5 group
    /// \param datasetname Name of dataset
    ///
    /// \returns Size of dimensions of dataset

    inline int get_dataset_size ( H5GroupPtr group_ptr,
                                  const std::string& datasetname )
    {
        hsize_t size;

        // Get Dataset pointer
        H5DatasetPtr dataset_ptr;
        dataset_ptr.reset ( new H5Dataset ( group_ptr, datasetname ) );

        // Get Dataspace of Dataset
        H5Dataspace dataspace ( dataset_ptr->id ( ) );

        // Get size of Dataspace
        H5Sget_simple_extent_dims ( dataspace.id ( ), &size, NULL );

        return int(size );
    }

    /// \brief Write vector into dataset in file.
    ///
    /// \param group_ptr Pointer to HDF5 group
    /// \param datasetname Name of dataset
    /// \param array Buffer with vector of data to be written to the file

    template<typename T>
    void write_array ( H5GroupPtr group_ptr,
                       const std::string& datasetname,
                       const std::vector<T>& array )
    {
        H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, array.size ( ),
                                                   datasetname, "w", &array[0] ) );
        dataset_ptr->write ( array.size ( ), 0, &array[0] );
    }

    /// \brief Write array into dataset in file.
    ///
    /// \param group_ptr Pointer to HDF5 group
    /// \param datasetname Name of dataset
    /// \param array Buffer with array of data to be written to the file
    /// \param size Size of array

    template<typename T>
    void write_array ( H5GroupPtr group_ptr,
                       const std::string& datasetname,
                       const T* array,
                       int size )
    {
        H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, size,
                                                   datasetname, "w", array ) );
        dataset_ptr->write ( size, 0, array );
    }

    //helper function 

    template<typename T>
    inline void write_array_parallel ( H5GroupPtr group_ptr, std::string datasetname, const T* data, int data_size, const MPI_Comm& comm )
    {
        int rank, num_part;
        MPI_Comm_rank ( comm, &rank );
        MPI_Comm_size ( comm, &num_part );

        int local_data_size = data_size;
        int data_offset = 0;
        if ( rank > 0 )
        {
            MPI_Status status;
            MPI_Recv ( &data_offset, 1, MPI_INT, rank - 1, rank - 1, comm, &status );
        }
        int data_next_offset = data_offset + local_data_size;
        if ( rank < num_part - 1 )
        {
            MPI_Send ( &data_next_offset, 1, MPI_INT, ( rank + 1 ), rank, comm );
        }
        int global_data_size = data_next_offset;
        MPI_Bcast ( &global_data_size, 1, MPI_INT, num_part - 1, comm );

        if ( global_data_size > 0 )
        {
            H5DatasetPtr data_dataset_ptr ( new H5Dataset ( group_ptr, global_data_size,
                                                            datasetname, "w", data ) );
            data_dataset_ptr->write ( local_data_size, data_offset, data );
        }
        else
        {
            //std::cerr << "Could not write dataset " << datasetname << " because it has no elements!" <<std::endl;
        }
        std::stringstream data_offset_name;
        data_offset_name << datasetname << "_offset_";
        H5DatasetPtr data_offset_ptr ( new H5Dataset ( group_ptr, num_part,
                                                       data_offset_name.str ( ), "w", &data_next_offset ) );
        data_offset_ptr->write ( 1, rank, &data_next_offset );
    }

    //helper function 

    template<typename T>
    inline void write_array_parallel ( H5GroupPtr group_ptr, std::string datasetname, const std::vector<T>& data, const MPI_Comm& comm )
    {
        write_array_parallel ( group_ptr, datasetname, &data[0], data.size ( ), comm );
    }

    template<typename KeyT, typename ValueT>
    inline void write_map_parallel ( H5GroupPtr group_ptr, std::string datasetname, const std::map< KeyT, ValueT>& map_data, const MPI_Comm& comm )
    {
        int rank, num_part;
        MPI_Comm_rank ( comm, &rank );
        MPI_Comm_size ( comm, &num_part );

        typename std::vector< std::pair<KeyT, ValueT> > temp;
        typename std::map<KeyT, ValueT>::const_iterator it;

        for ( it = map_data.begin ( ); it != map_data.end ( ); it++ )
            temp.push_back ( *it );

        hsize_t size = temp.size ( );

        write_array_parallel ( group_ptr, datasetname, temp, comm );
    }

    //helper function 

    template<typename T>
    void read_array_parallel ( H5GroupPtr group_ptr, std::string datasetname, T *&data, int& data_size, const MPI_Comm& comm )
    {
        int rank, num_part;
        MPI_Comm_rank ( comm, &rank );
        MPI_Comm_size ( comm, &num_part );

        std::stringstream data_offset_name;
        data_offset_name << datasetname << "_offset_";
        int data_next_offset;
        H5DatasetPtr data_offset_ptr ( new H5Dataset ( group_ptr, num_part,
                                                       data_offset_name.str ( ), "r", &data_next_offset ) );
        data_offset_ptr->read ( 1, rank, &data_next_offset );

        int data_offset = 0;
        if ( rank > 0 )
        {
            MPI_Status status;
            MPI_Recv ( &data_offset, 1, MPI_INT, rank - 1, rank - 1, comm, &status );
        }
        if ( rank < num_part - 1 )
        {
            MPI_Send ( &data_next_offset, 1, MPI_INT, ( rank + 1 ), rank, comm );
        }
        int local_data_size = data_next_offset - data_offset;
        int global_data_size = data_next_offset;
        MPI_Bcast ( &global_data_size, 1, MPI_INT, num_part - 1, comm );

        data_size = local_data_size;
        if ( global_data_size > 0 )
        {
            data = ( T * ) malloc ( sizeof (T ) * ( data_size ) );
            H5DatasetPtr data_dataset_ptr ( new H5Dataset ( group_ptr, global_data_size,
                                                            datasetname, "r", data ) );
            data_dataset_ptr->read ( local_data_size, data_offset, data );
        } //else set data = NULL???
    }

    //helper function 

    template<typename T>
    void read_array_parallel ( H5GroupPtr group_ptr, std::string datasetname, std::vector<T>& data, const MPI_Comm& comm )
    {
        int rank, num_part;
        MPI_Comm_rank ( comm, &rank );
        MPI_Comm_size ( comm, &num_part );

        std::stringstream data_offset_name;
        data_offset_name << datasetname << "_offset_";
        int data_next_offset;
        H5DatasetPtr data_offset_ptr ( new H5Dataset ( group_ptr, num_part,
                                                       data_offset_name.str ( ), "r", &data_next_offset ) );
        data_offset_ptr->read ( 1, rank, &data_next_offset );
        int data_offset = 0;
        if ( rank > 0 )
        {
            MPI_Status status;
            MPI_Recv ( &data_offset, 1, MPI_INT, rank - 1, rank - 1, comm, &status );
        }
        if ( rank < num_part - 1 )
        {
            MPI_Send ( &data_next_offset, 1, MPI_INT, ( rank + 1 ), rank, comm );
        }
        int local_data_size = data_next_offset - data_offset;
        int global_data_size = data_next_offset;
        MPI_Bcast ( &global_data_size, 1, MPI_INT, num_part - 1, comm );

        data.resize ( local_data_size );
        if ( global_data_size > 0 )
        {
            H5DatasetPtr data_dataset_ptr ( new H5Dataset ( group_ptr, global_data_size,
                                                            datasetname, "r", &data[0] ) );
            data_dataset_ptr->read ( local_data_size, data_offset, &data[0] );
        }
    }

    /*void read_array_parallel(H5GroupPtr group_ptr, std::string datasetname, std::vector<std::string>& data, MPI_Comm& comm)
    {
        std::vector<StringContainer> tc = convert ( data );
        read_array_parallel(group_ptr, datasetname, tc, comm);
    }*/
    /// \brief Write single value into dataset in file.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of the dataset
    /// \param value Buffer with single value to be written to the file

    template<typename T>
    void write_value ( H5GroupPtr group_ptr,
                       const std::string& datasetname,
                       const T& value )
    {
        write_array ( group_ptr, datasetname, &value, 1 );
    }

    /// \brief Convert vector of strings to vector of character array.
    ///
    /// \param v Vector of strings to be converted
    /// \returns Vector of character arrays
    ///
    /// HDF5 does not support writing std::strings. To do so, they must be
    /// converted into character arrays first. Maximum string size is 512 Byte.

    inline std::vector<StringContainer> convert ( const std::vector<std::string> &v )
    {
        std::vector<StringContainer> tc;
        int max = get_max_stringlength ( v );
        for ( std::vector<std::string>::const_iterator
              i = v.begin ( ), end = v.end ( ); i != end; ++i )
        {
            StringContainer t;
            memset ( &t, 0, sizeof (StringContainer ) );
            strncpy ( t.string, i->c_str ( ), max );
            tc.push_back ( t );
        }
        return tc;
    }

    /// \brief Write vector of strings into dataset in file.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param array Buffer with vector of strings to be written to file
    ///
    /// Strings must be converted to character arrays. Then they can be
    /// written out easily.

    inline void write_array ( H5GroupPtr group_ptr,
                              const std::string& datasetname,
                              const std::vector<std::string>& v )
    {

        H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, v.size ( ),
                                                   datasetname, "w", &v ) );

        std::vector<StringContainer> tc = convert ( v );

        dataset_ptr->write ( v.size ( ), 0, &tc[0] );
    }

    /// \brief Write array of strings into dataset in file.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param array Buffer with array of strings to be written to file
    /// \param size Size of array
    ///
    /// First, array of strings is converted to vector of strings of the same size.
    /// Then write_array function for vectors is called.

    inline void write_array ( H5GroupPtr group_ptr,
                              const std::string& datasetname,
                              const std::string* value,
                              int size )
    {
        std::vector<std::string> v;
        v.resize ( size );
        for ( int i = 0; i < size; i++ )
        {
            v[i] = value[i];
        }

        write_array ( group_ptr, datasetname, v );
    }

    inline void write_array_parallel ( H5GroupPtr group_ptr, std::string datasetname, const std::vector<std::string>& data, const MPI_Comm& comm )
    {
        std::vector<StringContainer> tc = convert ( data );
        write_array_parallel ( group_ptr, datasetname, tc, comm );
    }
    /// \brief Write single string into dataset in file.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param array Buffer with single string to be written to file
    ///
    /// First, single string is converted to vector with strings of size 1.
    /// Then write_array function for vectors is called.

    inline void write_value ( H5GroupPtr group_ptr,
                              const std::string& datasetname,
                              const std::string& value )
    {
        std::vector<std::string> v ( 1, value );
        write_array ( group_ptr, datasetname, v );
    }

    /// \brief Read vector from dataset in file into buffer.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param array Buffer to receive vector of data of dataset read from file

    template<typename T>
    void read_array ( H5GroupPtr group_ptr,
                      const std::string& datasetname,
                      std::vector<T>& array )
    {
        //Resize vector to fit elements
        int size = get_dataset_size ( group_ptr, datasetname );
        array.resize ( size );

        H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, size,
                                                   datasetname, "r", &array[0] ) );
        dataset_ptr->read ( size, 0, &array[0] );

    }

    /// \brief Read array from dataset in file into buffer.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param array Buffer to receive array of data of dataset read from file
    ///
    /// Reading arrays is based on the read_array function for vectors.
    /// A temporary vector is used for reading data from file. Data in
    /// vector is written into array afterwards.

    template<typename T>
    void read_array ( H5GroupPtr group_ptr,
                      const std::string& datasetname,
                      T* array )
    {
        std::vector<T> v;
        read_array ( group_ptr, datasetname, v );
        for ( int i = 0; i < v.size ( ); i++ ) array[i] = v[i];
    }

    /// \brief Read single value from dataset in file into buffer.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param value Buffer to receive data of dataset read from file
    ///
    /// Reading single values is based on the read_array function for arrays.

    template<typename T>
    void read_value ( H5GroupPtr group_ptr,
                      const std::string& datasetname,
                      T* value )
    {
        read_array ( group_ptr, datasetname, value );
    }

    /// \brief Read vector of strings from dataset in file into buffer.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param v Buffer to receive vector of strings of dataset read from file

    inline void read_array ( H5GroupPtr group_ptr,
                             const std::string& datasetname,
                             std::vector<std::string>& v )
    {
        // Set right size for v
        hsize_t size = get_dataset_size ( group_ptr, datasetname );
        v.resize ( size );

        H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, size,
                                                   datasetname, "r", &v[0] ) );

        H5Dataspace memspace ( size );

        std::vector <StringContainer> buffer;
        buffer.resize ( size );

        // Read Data from file.
        H5Datatype datatype = H5Dget_type ( dataset_ptr->id ( ) );

        H5Dataspace filespace ( dataset_ptr->id ( ) );

        H5Dread ( dataset_ptr->id ( ), datatype.id ( ), memspace.id ( ),
                  filespace.id ( ), H5P_DEFAULT, &buffer[0] );

        for ( size_t i = 0; i != size; ++i )
        {
            v[i] = std::string ( &buffer[i].string[0] );
        }

    }

    /// \brief Read array of strings from dataset in file into buffer.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param value Buffer to receive array of strings of dataset read from file
    ///
    /// Reading arrays of strings is based on the read_array function for vectors of strings.
    /// A temporary vector is used for reading data from file. Data in
    /// vector is written into array afterwards.

    inline void read_array ( H5GroupPtr group_ptr,
                             const std::string& datasetname,
                             std::string* value )
    {

        std::vector<std::string> v;
        read_array ( group_ptr, datasetname, v );
        for ( size_t i = 0; i != v.size ( ); ++i ) value[i] = v[i];

    }

    inline void read_array_parallel ( H5GroupPtr group_ptr, std::string datasetname, std::vector<std::string>& data, const MPI_Comm& comm )
    {
        //std::vector<StringContainer> tc(1);
        std::vector<StringContainer> tc ( 4 );
        for ( int i = 0; i < MaxStrLength; ++i )
            tc[0].string[i] = 'a';
        read_array_parallel<StringContainer>( group_ptr, datasetname, tc, comm );
        data.resize ( tc.size ( ) );
        for ( int i = 0; i < tc.size ( ); ++i )
        {
            data[i].assign ( tc[i].string );
        }

    }

    /// \brief Read single string from dataset in file into buffer.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param value Buffer to receive string of dataset read from file
    ///
    /// Reading single values is based on the read_array function for arrays of strings.

    inline void read_value ( H5GroupPtr group_ptr,
                             const std::string& datasetname,
                             std::string& value )
    {
        read_array ( group_ptr, datasetname, &value );
    }

    /// \brief Convert string to character array.
    ///
    /// \param string
    /// \param max Size of character array
    ///
    /// \return Character array

    inline StringContainer string2char ( const std::string &v, int max )
    {
        StringContainer tc;
        memset ( &tc, 0, sizeof (StringContainer ) );
        strncpy ( &tc.string[0], v.c_str ( ), max );

        return tc;
    }

    /// \brief Prepare write map into dataset in file.
    ///
    /// \param group_ptr Pointer to an HDF5 group
    /// \param datasetname Name of dataset
    /// \param m Map to be written out.
    ///
    /// Key and value of the map are transformed into a vector of pairs.
    /// The function write_map_content writes its content to a file.

    template<typename KeyT, typename ValueT>
    void write_map ( H5GroupPtr group_ptr,
                     const std::string& datasetname,
                     const std::map<KeyT, ValueT>& m )
    {

        // Transform map into a vector of pairs.
        typename std::vector< std::pair<KeyT, ValueT> > temp;
        typename std::map<KeyT, ValueT>::const_iterator it;

        for ( it = m.begin ( ); it != m.end ( ); it++ )
            temp.push_back ( *it );

        hsize_t size = temp.size ( );

        H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, size,
                                                   datasetname, "w", &temp ) );
        write_map_content ( dataset_ptr, temp );

    }

    /// \brief Write content of map into dataset in file.
    ///
    /// \param dataset_ptr Pointer to an HDF5 dataset
    /// \param temp Vector of pairs to be written out.
    ///
    /// Standard write_map_content. Not usable for a map containing strings.

    template<typename KeyT, typename ValueT>
    void write_map_content ( H5DatasetPtr dataset_ptr,
                             const std::vector< std::pair<KeyT, ValueT> > &temp )
    {
        hsize_t size = temp.size ( );
        dataset_ptr->write ( size, 0, &temp[0] );
    }

    /// \brief Write content of map into dataset in file with string value.
    ///
    /// \param dataset_ptr Pointer to an HDF5 dataset
    /// \param temp Vector of pairs to be written out.

    template<typename KeyT>
    void write_map_content ( H5DatasetPtr dataset_ptr,
                             const std::vector< std::pair<KeyT, std::string> > &temp )
    {

        typename std::vector< std::pair<KeyT, StringContainer> > temp_string;
        hsize_t size = temp.size ( );
        temp_string.resize ( size );
        int max = 0;
        for ( int i = 0; i < size; i++ )
        {
            if ( max < get_stringlength ( ( temp[i] ).second ) ) max = ( ( temp[i] ).second ).size ( );
        }

        for ( int i = 0; i < size; i++ )
        {
            ( temp_string[i] ).first = ( temp[i] ).first;
            ( temp_string[i] ).second = string2char ( ( temp[i] ).second, max );
        }
        dataset_ptr->write ( size, 0, &temp_string[0] );
    }

    /// \brief Write content of map into dataset in file with string key.
    ///
    /// \param dataset_ptr Pointer to an HDF5 dataset
    /// \param temp Vector of pairs to be written out.

    template<typename ValueT>
    void write_map_content ( H5DatasetPtr dataset_ptr,
                             const std::vector< std::pair<std::string, ValueT> > &temp )
    {
        typename std::vector< std::pair<StringContainer, ValueT> > temp_string;
        hsize_t size = temp.size ( );
        temp_string.resize ( size );

        int max = 0;
        for ( int i = 0; i < size; i++ )
        {
            if ( max < get_stringlength ( ( temp[i] ).first ) ) max = ( ( temp[i] ).first ).size ( );
        }

        for ( int i = 0; i < size; i++ )
        {
            ( temp_string[i] ).first = string2char ( ( temp[i] ).first, max );
            ( temp_string[i] ).second = ( temp[i] ).second;
        }

        dataset_ptr->write ( size, 0, &temp_string[0] );

    }

    /// \brief Write content of map into dataset in file with string key and value.
    ///
    /// \param dataset_ptr Pointer to an HDF5 dataset
    /// \param temp Vector of pairs to be written out.

    inline void write_map_content ( H5DatasetPtr dataset_ptr,
                                    const std::vector< std::pair<std::string, std::string> > &temp )
    {
        std::vector< std::pair<StringContainer, StringContainer> > temp_string;
        hsize_t size = temp.size ( );
        temp_string.resize ( size );

        std::size_t max1 = 0;
        std::size_t max2 = 0;
        for ( size_t i = 0; i != size; ++i )
        {
            max1 = std::max ( max1, temp[i].first.size ( ) );
            max2 = std::max ( max2, temp[i].second.size ( ) );
        }

        for ( size_t i = 0; i != size; ++i )
        {
            ( temp_string[i] ).first = string2char ( ( temp[i] ).first, max1 );
            ( temp_string[i] ).second = string2char ( ( temp[i] ).second, max2 );
        }
        dataset_ptr->write ( size, 0, &temp_string[0] );
    }

    /// \brief Prepare read map from dataset in file into buffer.
    ///
    /// \param group_ptr Pointer to HDF5 group
    /// \param datasetname Name of dataset
    /// \param m Buffer for map to be read from file
    ///
    /// Key and value of the map must be read into a vector of pairs.
    /// Then they are transformed into map.
    /// The function read_map_content reads the map from file into a vector of pairs.

    template<typename KeyT, typename ValueT>
    void read_map ( H5GroupPtr group_ptr,
                    const std::string& datasetname,
                    std::map<KeyT, ValueT>& m )
    {

        // Prepare Buffer
        typename std::vector< std::pair<KeyT, ValueT> > buffer;
        hsize_t size = get_dataset_size ( group_ptr, datasetname );
        buffer.resize ( size );

        // Create DatasetPtr
        H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, size,
                                                   datasetname, "r", &buffer ) );

        read_map_content ( dataset_ptr, size, buffer );

        // Write vector<pair> into map
        for ( int i = 0; i < size; i++ )
            m[buffer[i].first] = buffer[i].second;
    }

    template<typename KeyT, typename ValueT>
    void read_map_parallel ( H5GroupPtr group_ptr,
                             const std::string& datasetname,
                             std::map<KeyT, ValueT>& map_data, const MPI_Comm& comm )
    {
        int rank, num_part;
        MPI_Comm_rank ( comm, &rank );
        MPI_Comm_size ( comm, &num_part );

        typename std::vector< std::pair<KeyT, ValueT> > temp;

        read_array_parallel ( group_ptr, datasetname, temp, comm );

        typename std::vector< std::pair<KeyT, ValueT> >::iterator it;
        for ( it = temp.begin ( ); it != temp.end ( ); ++it )
        {
            map_data.insert ( *it );
        }
    }

    /// \brief Read content of map from dataset in file into vector of pairs.
    ///
    /// \param dataset_ptr Pointer to HDF5 dataset
    /// \param temp Buffer for vector of pairs to be read from file.
    ///
    /// Standard read_map_content. Not usable for a map containing strings.

    template<typename KeyT, typename ValueT>
    void read_map_content ( H5DatasetPtr dataset_ptr, hsize_t size,
                            std::vector< std::pair<KeyT, ValueT> > &temp )
    {

        // Read data from file into vector<pair>
        H5Datatype datatype = H5Dget_type ( dataset_ptr->id ( ) );
        H5Dataspace memspace ( size );
        H5Dataspace filespace ( dataset_ptr->id ( ) );

        H5Dread ( dataset_ptr->id ( ), datatype.id ( ), memspace.id ( ),
                  filespace.id ( ), H5P_DEFAULT, &temp[0] );
    }

    /// \brief Read content of map from dataset in file into vector of pairs with string key.
    ///
    /// \param dataset_ptr Pointer to HDF5 dataset
    /// \param temp Buffer for vector of pairs to be read from file.

    template<typename KeyT>
    void read_map_content ( H5DatasetPtr dataset_ptr, hsize_t size,
                            std::vector< std::pair<KeyT, std::string> > &temp )
    {

        H5Datatype datatype = H5Dget_type ( dataset_ptr->id ( ) );
        H5Dataspace memspace ( size );
        H5Dataspace filespace ( dataset_ptr->id ( ) );

        typename std::vector< std::pair<KeyT, StringContainer> > buffer;
        buffer.resize ( size );

        // Read Data from File.
        H5Dread ( dataset_ptr->id ( ), datatype.id ( ), memspace.id ( ),
                  filespace.id ( ), H5P_DEFAULT, &buffer[0] );

        for ( int i = 0; i < size; i++ )
        {
            ( temp[i] ).first = ( buffer[i] ).first;
            ( temp[i] ).second = std::string ( &( ( buffer[i] ).second ).string[0] );
        }

    }

    /// \brief Read content of map from dataset in file into vector of pairs with string value.
    ///
    /// \param dataset_ptr Pointer to HDF5 dataset
    /// \param temp Buffer for vector of pairs to be read from file.

    template<typename ValueT>
    void read_map_content ( H5DatasetPtr dataset_ptr, hsize_t size,
                            std::vector< std::pair<std::string, ValueT> > &temp )
    {

        H5Datatype datatype = H5Dget_type ( dataset_ptr->id ( ) );
        H5Dataspace memspace ( size );
        H5Dataspace filespace ( dataset_ptr->id ( ) );

        typename std::vector< std::pair<StringContainer, ValueT> > buffer;
        buffer.resize ( size );

        // Read Data from File.
        H5Dread ( dataset_ptr->id ( ), datatype.id ( ), memspace.id ( ),
                  filespace.id ( ), H5P_DEFAULT, &buffer[0] );

        for ( int i = 0; i < size; i++ )
        {
            ( temp[i] ).first = std::string ( &( ( buffer[i] ).first ).string[0] );
            ( temp[i] ).second = ( buffer[i] ).second;
        }

    }

    /// \brief Read content of map from dataset in file into vector of pairs with string key and value.
    ///
    /// \param dataset_ptr Pointer to HDF5 dataset
    /// \param temp Buffer for vector of pairs to be read from file.

    inline void read_map_content ( H5DatasetPtr dataset_ptr, hsize_t size,
                                   std::vector< std::pair<std::string, std::string> > &temp )
    {
        H5Datatype datatype = H5Dget_type ( dataset_ptr->id ( ) );
        H5Dataspace memspace ( size );
        H5Dataspace filespace ( dataset_ptr->id ( ) );

        std::vector< std::pair<StringContainer, StringContainer> > buffer;
        buffer.resize ( size );

        // Read Data from File.
        H5Dread ( dataset_ptr->id ( ), datatype.id ( ), memspace.id ( ),
                  filespace.id ( ), H5P_DEFAULT, &buffer[0] );

        for ( size_t i = 0; i != size; ++i )
        {
            ( temp[i] ).first = std::string ( &( ( buffer[i] ).first ).string[0] );
            ( temp[i] ).second = std::string ( &( ( buffer[i] ).second ).string[0] );
        }

    }

}

#    endif // WITH_HDF5
#endif // include guards
