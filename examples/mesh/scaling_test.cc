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

/// \author Philipp Gerstner

#include "scaling_test.h"

const int DEBUG_LEVEL = 1;
const int PRINT_PROC = 1;
static const char* PARAM_FILENAME = "scaling_test.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class ScalingTest
{
  public:

    ScalingTest ( const std::string& param_filename,
                  const std::string& path_mesh,
                  const std::string& csv_path )
    : path_mesh ( path_mesh ),
    csv_path_ ( csv_path ),
    comm_ ( MPI_COMM_WORLD ),
    rank_ ( -1 ),
    num_partitions_ ( -1 ),
    params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
    is_done_ ( false ),
    refinement_level_ ( 0 )
    {
        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_partitions_ );
    }

    // Main algorithm

    void run ( )
    {
        if ( csv_path_.empty ( ) )
        {
            csv_path_ = ".";
        }
        // Construct / read in the initial mesh.
        this->csv_names_.clear ( );
        this->csv_names_.push_back ( "Read time" );
        this->csv_names_.push_back ( "Read  in cells" );
        this->csv_names_.push_back ( "Sequential refine level " );
        this->csv_names_.push_back ( "Initial cells " );
        this->csv_names_.push_back ( "Sequential refine time" );
        this->csv_names_.push_back ( "Partition time" );

        this->read_mesh ( );

        this->partition_mesh ( );

        if ( this->rank_ == MASTER_RANK )
        {
            std::stringstream csv_file_name;
#ifdef USE_MESH_P4EST
            csv_file_name << csv_path_ << "/pXest_mesh_" << DIMENSION << "D_" << this->num_partitions_ << "_init_time.csv";
#else
            csv_file_name << csv_path_ << "/dbview_mesh_" << DIMENSION << "D_" << this->num_partitions_ << "_init_time.csv";
#endif
            this->csv_writer_.InitFilename ( csv_file_name.str ( ) );
            this->csv_writer_.Init ( this->csv_names_ );
            this->csv_writer_.write ( this->csv_quantities_ );
        }

        // Prepare headers of CSV file
        this->csv_names_.clear ( );
        this->csv_names_.push_back ( "Level" );
        this->csv_names_.push_back ( "Number of cells " );
        this->csv_names_.push_back ( "Refine time " );
        this->csv_names_.push_back ( "Compute ghost time " );
        this->csv_names_.push_back ( "Number of ghost cells " );
        this->csv_names_.push_back ( "Build Iflist time " );
        this->csv_names_.push_back ( "Copy mesh time " );
        this->csv_names_.push_back ( "Compute parents time " );

        if ( this->rank_ == MASTER_RANK )
        {
            std::stringstream csv_file_name;
#ifdef USE_MESH_P4EST
            csv_file_name << csv_path_ << "/pXest_mesh_" << DIMENSION << "D_" << this->num_partitions_ << "_run_time.csv";
#else
            csv_file_name << csv_path_ << "/dbview_mesh_" << DIMENSION << "D_" << this->num_partitions_ << "_run_time.csv";
#endif
            this->csv_writer_.InitFilename ( csv_file_name.str ( ) );
            this->csv_writer_.Init ( this->csv_names_ );
        }

        const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
        const int ref_type = params_["Mesh"]["Refinement"].get<int>( 1 );

        // Main adaptation loop.
        while ( !is_done_ )
        {
            if ( this->refinement_level_ >= final_ref_level )
            {
                break;
            }
            // Prepare data for statistics in CSV file
            this->csv_quantities_.clear ( );

            // Add current refinement level
            this->csv_quantities_.push_back ( this->refinement_level_ + 1 );

            LOG_INFO ( "Refinement Level", " ===== " << refinement_level_ + 1 << "\n=====" );

            if ( ref_type == 1 )
            {
                this->refine_mesh_uniform ( );
            }
            else
            {
                this->refine_mesh_local ( );
            }

            this->compute_ghosts ( );

            this->build_interface_list ( );

            this->copy_mesh ( );

            this->compute_parents ( );

            // Write CSV output
            if ( this->rank_ == MASTER_RANK )
            {
                this->csv_writer_.write ( this->csv_quantities_ );
            }
        }
    }

    ~ScalingTest ( )
    {
    }

  private:
    // Member functions

    std::string path_mesh;
    std::string csv_path_;

    void read_mesh ( );
    void partition_mesh ( );
    void refine_mesh_uniform ( );
    void refine_mesh_local ( );
    void compute_ghosts ( );
    void copy_mesh ( );
    void build_interface_list ( );
    void compute_parents ( );

    // method of adaptive refinement
    int refinement_;
    // member variables
    // MPI communicator.
    MPI_Comm comm_;
    // Local process rank and number of processes.
    int rank_, num_partitions_;

    // Parameter data read in from file.
    PropertyTree params_;

    // Local mesh and mesh on master process.
    MeshPtr mesh_, master_mesh_, mesh_without_ghost_;
    // Solution space.
    VectorSpace<double> space_;

    // Flag for stopping adaptive loop.
    bool is_done_;
    // Current refinement level.
    int refinement_level_;
    int num_initial_cells_;

    /// CSV writer
    CSVWriter<double> csv_writer_;
    /// Headers of columns in CSV file
    std::vector<std::string> csv_names_;
    /// Data  in CSV file (filled and written to file after
    /// every level)
    std::vector<double> csv_quantities_;

}; // end class ScalingTest

// Program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    // Set default parameter file
    std::string param_filename ( PARAM_FILENAME );
    std::string path_mesh;
    std::string csv_path;
    // If set take parameter file specified on console
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }
    // If set take mesh following path specified on console
    if ( argc > 2 )
    {
        path_mesh = std::string ( argv[2] );
    }
    // If set take mesh following path specified on console
    if ( argc > 3 )
    {
        csv_path = std::string ( argv[3] );
    }

    try
    {
        // Create log files for INFO and DEBUG output
        //    std::ofstream info_log("poisson_tutorial_info_log");
        LogKeeper::get_log ( "info" ).set_target ( &( std::cout ) );
        //std::ofstream debug_log("poisson_tutorial_debug_log");
        LogKeeper::get_log ( "debug" ).set_target ( &( std::cout ) );
        //     std::ofstream error_log ( "poisson_tutorial_error_log" );
        LogKeeper::get_log ( "error" ).set_target ( &( std::cout ) );

        // Create application object and run it
        ScalingTest app ( param_filename, path_mesh, csv_path );

        app.run ( );

    }
    catch ( std::exception& e )
    {
        std::cerr << "\nProgram ended with uncaught exception.\n";
        std::cerr << e.what ( ) << "\n";
        return -1;
    }
    MPI_Finalize ( );
    return 0;
}

void ScalingTest::read_mesh ( )
{
    // Read in the mesh.
    //The mesh is chosen according to the dimension of the problem.
    std::string mesh_name;
    Timer timer;

    if ( DIMENSION == 2 )
    {
        mesh_name = params_["Mesh"]["Filename2"].get<std::string>( "unit_square_inner_square.inp" );
    }
    if ( DIMENSION == 3 )
    {
        mesh_name = params_["Mesh"]["Filename3"].get<std::string>( "unit_cube.inp" );
    }

    std::string mesh_filename;
    if ( path_mesh.empty ( ) )
    {
        mesh_filename = std::string ( DATADIR ) + mesh_name;
    }
    else
    {
        mesh_filename = path_mesh + mesh_name;
    }

    std::vector<MasterSlave> period ( 0, MasterSlave ( 0., 0., 0., 0 ) );
#ifdef PARALLEL_READ
    timer.start ( );
#    ifdef USE_MESH_P4EST
    master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0, period, mesh::IMPL_P4EST );
#    endif

    timer.stop ( );
    if ( rank_ == MASTER_RANK )
    {
        this->csv_quantities_.push_back ( timer.get_duration ( ) );
        int num_read_cells = master_mesh_->num_entities ( DIMENSION );
        this->csv_quantities_.push_back ( num_read_cells );
    }

    timer.reset ( );
    timer.start ( );

    // Refine the mesh until the initial refinement level is reached.
    const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( 3 );
    if ( initial_ref_lvl > 0 )
    {
        LOG_DEBUG ( 1, "Initial refinement level " << initial_ref_lvl );
        master_mesh_ = master_mesh_->refine_uniform_seq ( initial_ref_lvl );
    }
    refinement_level_ = initial_ref_lvl;
    timer.stop ( );
    double ref_time = timer.get_duration ( );
    if ( rank_ == MASTER_RANK )
    {
        this->csv_quantities_.push_back ( initial_ref_lvl );
        this->num_initial_cells_ = master_mesh_->num_entities ( DIMENSION );
        this->csv_quantities_.push_back ( this->num_initial_cells_ );
        this->csv_quantities_.push_back ( ref_time );
    }
#else
    if ( rank_ == MASTER_RANK )
    {
        timer.start ( );
#    ifdef USE_MESH_P4EST
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0, period, mesh::IMPL_P4EST );
#    else
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0, period, mesh::IMPL_DBVIEW );
#    endif

        timer.stop ( );
        this->csv_quantities_.push_back ( timer.get_duration ( ) );
        int num_read_cells = master_mesh_->num_entities ( DIMENSION );
        this->csv_quantities_.push_back ( num_read_cells );

        timer.reset ( );
        timer.start ( );
        // Refine the mesh until the initial refinement level is reached.
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( 3 );
        if ( initial_ref_lvl > 0 )
        {
            LOG_DEBUG ( 1, "Initial refinement level " << initial_ref_lvl );
            master_mesh_ = master_mesh_->refine_uniform_seq ( initial_ref_lvl );
        }
        refinement_level_ = initial_ref_lvl;
        timer.stop ( );
        double ref_time = timer.get_duration ( );
        this->csv_quantities_.push_back ( initial_ref_lvl );
        this->num_initial_cells_ = master_mesh_->num_entities ( DIMENSION );
        this->csv_quantities_.push_back ( this->num_initial_cells_ );
        this->csv_quantities_.push_back ( ref_time );
    }
    // Broadcast information from master to slaves.
    MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );
    MPI_Bcast ( &num_initial_cells_, 1, MPI_INT, MASTER_RANK, comm_ );
#endif
}

void ScalingTest::partition_mesh ( )
{
    int uniform_ref_steps;
    Timer timer;
    timer.start ( );
#ifdef USE_MESH_P4EST
    mesh_ = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &uniform_ref_steps, mesh::IMPL_P4EST );
#else
    mesh_without_ghost_ = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &uniform_ref_steps, mesh::IMPL_DBVIEW );
#endif

    timer.stop ( );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );

    refinement_level_ += uniform_ref_steps;

    /*
    PVtkWriter writer ( comm_ );
    std::ostringstream name;
    name << "scaling_test_initial_mesh.pvtu";
    std::string output_file = name.str ( );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
     */
}

void ScalingTest::copy_mesh ( )
{
    Timer timer;
    timer.start ( );
#ifdef USE_MESH_P4EST
    MeshPtr inter_mesh ( new MeshPXest ( DIMENSION, DIMENSION ) );
#else
    MeshPtr inter_mesh ( new MeshDbView ( DIMENSION, DIMENSION ) );
#endif
    inter_mesh->deep_copy_from ( mesh_ );

    timer.stop ( );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );
}

void ScalingTest::compute_parents ( )
{
    Timer timer;
    timer.start ( );
    std::map<Id, std::vector<Id> > descendants;

    // Loop over all cells in fine mesh
    for ( EntityIterator cell_it = mesh_->begin ( DIMENSION ); cell_it != mesh_->end ( DIMENSION ); cell_it++ )
    {
        Id child_id = cell_it->id ( );
        EntityNumber child_index = cell_it->index ( );
        Id parent_id = mesh_->get_parent_cell_id ( child_index );

        if ( parent_id >= 0 )
        {
            std::map<Id, std::vector<Id> >::iterator it = descendants.find ( parent_id );
            if ( it == descendants.end ( ) )
            {
                std::vector<Id> tmp_id ( 1, child_id );
                descendants[parent_id] = tmp_id;
            }
            else
            {
                it->second.push_back ( child_id );
            }
        }
    }

    timer.stop ( );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );
}

void ScalingTest::build_interface_list ( )
{
    Timer timer;
    timer.start ( );

    InterfaceList if_list = InterfaceList::create ( mesh_ );

    timer.stop ( );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );
}

void ScalingTest::compute_ghosts ( )
{
    Timer timer;
    timer.start ( );

    SharedVertexTable shared_verts;
#ifdef USE_MESH_P4EST
    mesh_ = compute_ghost_cells ( *mesh_, comm_, shared_verts, mesh::IMPL_P4EST, 1 );
    timer.stop ( );
#else
    mesh_ = compute_ghost_cells ( *mesh_without_ghost_, comm_, shared_verts, mesh::IMPL_DBVIEW, 1 );
    timer.stop ( );

#endif
    const int init_ref_level = params_["Mesh"]["InitialRefLevel"].get<int>( 6 );
    int num_cells = static_cast < int > ( std::pow ( static_cast < double > ( 2 ), DIMENSION * ( refinement_level_ - init_ref_level ) ) ) * this->num_initial_cells_;

    //assert ( num_cells == mesh_->num_global_cells ( comm_ ) );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );

    int num_ghost_local = mesh_->num_ghost_cells ( );
    int num_ghost_global = 0;
    MPI_Allreduce ( &num_ghost_local, &num_ghost_global, 1, MPI_INT, MPI_SUM, comm_ );
    this->csv_quantities_.push_back ( num_ghost_global );

    PVtkWriter writer ( comm_ );
    std::ostringstream name;
    name << "scaling_test_mesh_" << refinement_level_ << ".pvtu";
    std::string output_file = name.str ( );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );

}

void ScalingTest::refine_mesh_uniform ( )
{
    const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
    const int init_ref_level = params_["Mesh"]["InitialRefLevel"].get<int>( 6 );
    if ( refinement_level_ >= final_ref_level )
    {
        is_done_ = true;
    }
    else
    {
        ++refinement_level_;
        Timer timer;
        timer.start ( );
#ifdef USE_MESH_P4EST
        mesh_ = mesh_->refine ( );
        timer.stop ( );
        int num_cells = mesh_->num_global_cells ( comm_ );
#else
        mesh_without_ghost_ = mesh_without_ghost_->refine ( );
        timer.stop ( );
        int num_cells = mesh_->num_global_cells ( comm_ );
#endif

        this->csv_quantities_.push_back ( num_cells );
        this->csv_quantities_.push_back ( timer.get_duration ( ) );
    }
}

void ScalingTest::refine_mesh_local ( )
{
    const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
    const int init_ref_level = params_["Mesh"]["InitialRefLevel"].get<int>( 6 );
    double refine_frac = params_["Mesh"]["FractionToRefine"].get<double>( 0.2 );
    double coarsen_frac = params_["Mesh"]["FractionToCoarsen"].get<double>( 0.1 );
    int threshold = params_["Mesh"]["CoarsenThreshold"].get<int>( 100 );
    int coarsen_marker = params_["Mesh"]["CoarsenFlag"].get<int>( 1 );
    int balance_mode = params_["Mesh"]["BalanceConnectionMode"].get<int>( 1 );

    std::vector<int> adapt_markers;
    std::vector<double> indicator;

    int num_cells = mesh_->num_entities ( DIMENSION );
    indicator.resize ( num_cells, 0. );
    for ( int c = 0; c < num_cells; ++c )
    {
        if ( this->mesh_->cell_is_local ( c ) )
        {
            indicator[c] = ( double ) c;
        }
        else
        {
            indicator[c] = 0.;
        }
    }

    local_fixed_fraction_strategy ( refine_frac, coarsen_frac, threshold, coarsen_marker, indicator, adapt_markers );

    if ( refinement_level_ >= final_ref_level )
    {
        is_done_ = true;
    }
    else
    {
        ++refinement_level_;
        Timer timer;
        timer.start ( );
#ifdef USE_MESH_P4EST
        boost::intrusive_ptr<MeshPXest> mesh_pXest = boost::static_pointer_cast<MeshPXest> ( mesh_ );
        mesh_pXest->set_connection_mode ( balance_mode );

        mesh_ = mesh_->refine ( adapt_markers );
        timer.stop ( );
        int num_cells = mesh_->num_global_cells ( comm_ );
#else
        mesh_without_ghost_ = mesh_without_ghost_->refine ( adapt_markers );
        timer.stop ( );
        int num_cells = static_cast < int > ( std::pow ( static_cast < double > ( 2 ), DIMENSION * ( refinement_level_ - init_ref_level ) ) ) * this->num_initial_cells_;
        assert ( num_cells == mesh_without_ghost_->num_global_cells ( comm_ ) );
#endif

        this->csv_quantities_.push_back ( num_cells );
        this->csv_quantities_.push_back ( timer.get_duration ( ) );
    }
}
