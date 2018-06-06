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

#include "cell_type.h"

#include "common/vector_algebra.h"
#include "types.h"

namespace hiflow
{
    namespace mesh
    {

        //////////////// Point ////////////////

        Point::Point ( )
        : CellType ( CellType::POINT, 0, 1 )
        {
        }

        void Point::create_regular_entities ( )
        {
        }

        void Point::create_refined_vertices ( )
        {
        }

        void Point::create_refined_cells ( )
        {
        }

        void Point::create_refinements ( )
        {
        }

        //////////////// Line ////////////////
        namespace
        {

            const int EDGES_OF_LINE[1][2] = {
                { 0, 1 }
            };
            const int LINE_SUB_CELLS[2][2] = {
                { 0, 2 },
                { 2, 1 }
            };
            const int LINE_REFINEMENT[2] = { 1, 2 };
        }

        Line::Line ( )
        : CellType ( CellType::LINE, 1, 2 )
        {
        }

        void Line::create_regular_entities ( )
        {
        }

        void Line::create_refined_vertices ( )
        {
            add_refined_vertex ( std::vector<int>( &EDGES_OF_LINE[0][0],
                                 &EDGES_OF_LINE[0][2] ) );
        }

        void Line::create_refined_cells ( )
        {
            add_refined_cell ( std::vector<int>( &LINE_SUB_CELLS[0][0],
                               &LINE_SUB_CELLS[0][2] ) );
            add_refined_cell ( std::vector<int>( &LINE_SUB_CELLS[1][0],
                               &LINE_SUB_CELLS[1][2] ) );
        }

        void Line::create_refinements ( )
        {
            add_refinement ( std::vector<int>( &LINE_REFINEMENT[0], &LINE_REFINEMENT[2] ) );
        }

        //////////////// Triangle ////////////////
        namespace
        { // anonymous namespace used for file-local variables

            const int EDGES_OF_TRIANGLE[3][2] = {
                { 0, 1 },
                { 1, 2 },
                { 2, 0 }
            };

            const int TRI_4TRI_SUB_CELLS[4][3] = {
                { 0, 3, 5 },
                { 1, 4, 3 },
                { 2, 5, 4 },
                { 3, 4, 5 }
            };

            const int TRI_BISECTION_CELLS[3][2][3] = {
                {
                    { 0, 1, 4 },
                    { 0, 4, 2 }
                },
                {
                    { 0, 1, 5 },
                    { 5, 1, 2 }
                },
                {
                    { 0, 3, 2 },
                    { 3, 1, 2 }
                }
            };

            const int TRI_4TRI_REFINEMENT[4] = { 1, 2, 3, 4 };
            const int TRI_BISECTION_REFINEMENTS[3][2] = {
                { 5, 6 },
                { 7, 8 },
                { 9, 10 }
            };

        }

        Triangle::Triangle ( )
        : CellType ( CellType::TRIANGLE, 2, 3 )
        {
        }

        void Triangle::create_regular_entities ( )
        {
            // add edges
            for ( int i = 0; i < 3; ++i )
            {
                add_regular_entity ( 1, std::vector<int>( &EDGES_OF_TRIANGLE[i][0],
                                     &EDGES_OF_TRIANGLE[i][2] ) );
            }
        }

        void Triangle::create_refined_vertices ( )
        {
            // mid-edge vertices
            for ( int i = 0; i < 3; ++i )
            {
                add_refined_vertex ( std::vector<int>( &EDGES_OF_TRIANGLE[i][0],
                                     &EDGES_OF_TRIANGLE[i][2] ) );
            }
        }

        void Triangle::create_refined_cells ( )
        {
            // tri -> 4 tri subcells (similar sub-triangles)
            for ( int i = 0; i < 4; ++i )
            {
                add_refined_cell ( std::vector<int>( &TRI_4TRI_SUB_CELLS[i][0],
                                   &TRI_4TRI_SUB_CELLS[i][3] ) );
            }

            // bisection by vertex 0 - 2
            for ( int v = 0; v < 3; ++v )
            {
                for ( int i = 0; i < 2; ++i )
                {
                    add_refined_cell ( std::vector<int>( &TRI_BISECTION_CELLS[v][i][0], &TRI_BISECTION_CELLS[v][i][3] ) );
                }
            }
        }

        void Triangle::create_refinements ( )
        {
            add_refinement ( std::vector<int>( &TRI_4TRI_REFINEMENT[0], &TRI_4TRI_REFINEMENT[4] ) );
            for ( int v = 0; v < 3; ++v )
            {
                add_refinement ( std::vector<int>( &TRI_BISECTION_REFINEMENTS[v][0], &TRI_BISECTION_REFINEMENTS[v][2] ) );
            }
        }

        bool Triangle::check_cell_geometry ( const std::vector<Coordinate>& coords, GDim gdim ) const
        {
            // Check if edge vectors obey "right-hand" rule.
            assert ( gdim == 2 );

            Vec<3, double> pt[3];
            int k = 0;
            for ( int i = 0; i < 3; ++i )
            {
                for ( int j = 0; j < gdim; ++j )
                {
                    pt[i][j] = coords[k++];
                }
                pt[2][2] = 0.;
            }

            Vec<3, double> v[2];
            for ( int i = 0; i < 2; ++i )
            {
                v[i] = pt[i + 1] - pt[0];
            }

            const double orientation = cross ( v[0], v[1] )[2];

            return orientation > 0.0;
        }

        //////////////// Quadrilateral ////////////////
        namespace
        {
            const int EDGES_OF_QUAD[4][2] = {
                {0, 1 },
                {1, 2 },
                {2, 3 },
                {3, 0 }
            };

            const int MID_QUAD_VERTEX[4] = { 0, 1, 2, 3 };

            const int QUAD_4QUAD_SUB_CELLS[4][4] = {
                { 0, 4, 8, 7 },
                { 4, 1, 5, 8 },
                { 8, 5, 2, 6 },
                { 7, 8, 6, 3 }
            };

            // horizontal quad refinment
            const int QUAD_2QUAD_H_SUB_CELLS[2][4] = {
                { 0, 1, 5, 7 },
                { 7, 5, 2, 3 }
            };

            // vertical quad refinement
            const int QUAD_2QUAD_V_SUB_CELLS[2][4] = {
                { 0, 4, 6, 3 },
                { 4, 1, 2, 6 }
            };

            // diagonal quad refinements
            const int QUAD_2TRI_02_SUB_CELLS[2][3] = {
                { 0, 1, 2 },
                { 3, 0, 2 }
            };

            const int QUAD_2TRI_13_SUB_CELLS[2][3] = {
                { 0, 1, 3 },
                { 3, 1, 2 }
            };

            const int QUAD_3TRI_SUB_CELLS[4][3][3] = {
                { // corners 0 and 1
                    { 0, 1, 6 },
                    { 6, 1, 2 },
                    { 3, 0, 6 }
                },
                { // corners 1 and 2
                    { 0, 1, 7 },
                    { 7, 1, 2 },
                    { 3, 7, 2 }
                },
                { // corners 2 and 3
                    { 0, 4, 3 },
                    { 3, 4, 2 },
                    { 4, 1, 2 }
                },
                { // corners 3 and 0
                    { 0, 1, 5 },
                    { 0, 5, 3 },
                    { 3, 5, 2 }
                }
            };

            const int QUAD_4QUAD_REFINEMENT[4] = { 1, 2, 3, 4 };
            const int QUAD_2QUAD_H_REFINEMENT[2] = { 5, 6 };
            const int QUAD_2QUAD_V_REFINEMENT[2] = { 7, 8 };
            const int QUAD_2TRI_02_REFINEMENT[2] = { 9, 10 };
            const int QUAD_2TRI_13_REFINEMENT[2] = { 11, 12 };
            const int QUAD_3_TRI_REFINEMENTS[4][3] = {
                { 13, 14, 15 },
                { 16, 17, 18 },
                { 19, 20, 21 },
                { 22, 23, 24 }
            };
        }

        Quadrilateral::Quadrilateral ( )
        : CellType ( CellType::QUADRILATERAL, 2, 4 )
        {
        }

        void Quadrilateral::create_regular_entities ( )
        {
            // add edges
            for ( int i = 0; i < 4; ++i )
            {
                add_regular_entity ( 1, std::vector<int>( &EDGES_OF_QUAD[i][0],
                                     &EDGES_OF_QUAD[i][2] ) );
            }
        }

        void Quadrilateral::create_refined_vertices ( )
        {
            // mid-edge vertices
            for ( int i = 0; i < 4; ++i )
            {
                add_refined_vertex ( std::vector<int>( &EDGES_OF_QUAD[i][0],
                                     &EDGES_OF_QUAD[i][2] ) );

            }

            // mid-cell vertex
            add_refined_vertex ( std::vector<int>( &MID_QUAD_VERTEX[0], &MID_QUAD_VERTEX[4] ) );
        }

        void Quadrilateral::create_refined_cells ( )
        {
            // quad -> 4 quad subcells (all squares)
            for ( int i = 0; i < 4; ++i )
            {
                add_refined_cell ( std::vector<int>( &QUAD_4QUAD_SUB_CELLS[i][0],
                                   &QUAD_4QUAD_SUB_CELLS[i][4] ) );
            }

            // quad -> 2 quad horizontal subcells
            for ( int i = 0; i < 2; ++i )
            {
                add_refined_cell ( std::vector<int>( &QUAD_2QUAD_H_SUB_CELLS[i][0],
                                   &QUAD_2QUAD_H_SUB_CELLS[i][4] ) );
            }

            // quad -> 2 quad vertical subcells
            for ( int i = 0; i < 2; ++i )
            {
                add_refined_cell ( std::vector<int>( &QUAD_2QUAD_V_SUB_CELLS[i][0],
                                   &QUAD_2QUAD_V_SUB_CELLS[i][4] ) );
            }

            // quad -> 2 triangles diagonal 0-2
            for ( int i = 0; i < 2; ++i )
            {
                add_refined_cell ( std::vector<int>( &QUAD_2TRI_02_SUB_CELLS[i][0],
                                   &QUAD_2TRI_02_SUB_CELLS[i][3] ) );
            }

            // quad -> 2 triangles diagonal 1-3
            for ( int i = 0; i < 2; ++i )
            {
                add_refined_cell ( std::vector<int>( &QUAD_2TRI_13_SUB_CELLS[i][0],
                                   &QUAD_2TRI_13_SUB_CELLS[i][3] ) );
            }

            // quad -> 3 triangles for each midpoint
            for ( int r = 0; r < 4; ++r )
            {
                for ( int i = 0; i < 3; ++i )
                {
                    add_refined_cell ( std::vector<int>( &QUAD_3TRI_SUB_CELLS[r][i][0], &QUAD_3TRI_SUB_CELLS[r][i][3] ) );
                }
            }
        }

        void Quadrilateral::create_refinements ( )
        {
            add_refinement ( std::vector<int>( &QUAD_4QUAD_REFINEMENT[0], &QUAD_4QUAD_REFINEMENT[4] ) );
            add_refinement ( std::vector<int>( &QUAD_2QUAD_H_REFINEMENT[0], &QUAD_2QUAD_H_REFINEMENT[2] ) );
            add_refinement ( std::vector<int>( &QUAD_2QUAD_V_REFINEMENT[0], &QUAD_2QUAD_V_REFINEMENT[2] ) );
            add_refinement ( std::vector<int>( &QUAD_2TRI_02_REFINEMENT[0], &QUAD_2TRI_02_REFINEMENT[2] ) );
            add_refinement ( std::vector<int>( &QUAD_2TRI_13_REFINEMENT[0], &QUAD_2TRI_13_REFINEMENT[2] ) );
            for ( int r = 0; r < 4; ++r )
            {
                add_refinement ( std::vector<int>( &QUAD_3_TRI_REFINEMENTS[r][0], &QUAD_3_TRI_REFINEMENTS[r][3] ) );
            }
        }

        //////////////// Tetrahedron ////////////////
        namespace
        {
            const int MIDPOINT_OF_TET[4] = { 0, 1, 2, 3 };

            const int EDGES_OF_TET[6][2] = {
                {0, 1 }, // 0
                {1, 2 }, // 1
                {2, 0 }, // 2
                {1, 3 }, // 3
                {3, 0 }, // 4
                {2, 3 } // 5
            };

            const int FACES_OF_TET[4][3] = {
                {0, 2, 1 }, // bottom 0
                {0, 1, 3 }, // front  1
                {0, 3, 2 }, // left   2
                {1, 2, 3 } // right  3
            };

            const int TET_8TET_SUB_CELLS[8][4] = {
                {0, 4, 6, 8 }, // 0 bottom left
                {4, 1, 5, 7 }, // 1 bottom right
                {6, 5, 2, 9 }, // 2 bottom back
                {4, 5, 6, 7 }, // 3 middle bottom
                {4, 7, 6, 8 }, // 4 middle front
                {6, 7, 9, 8 }, // 5 middle left
                {6, 7, 5, 9 }, // 6 middle right
                {8, 7, 9, 3 } // 7 top
            };

            const int TET_4HEX_SUB_CELLS[4][8] = {
                {0, 4, 10, 6, 8, 11, 14, 12 }, // left front
                {4, 1, 5, 10, 11, 7, 13, 14 }, // right front
                {10, 5, 2, 6, 14, 13, 9, 12 }, // back
                {11, 7, 13, 14, 8, 3, 9, 12 } // top
            };

            const int TET_8TET_REFINEMENT[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
            const int TET_4HEX_REFINEMENT[4] = { 9, 10, 11, 12 };
        }

        Tetrahedron::Tetrahedron ( )
        : CellType ( CellType::TETRAHEDRON, 3, 4 )
        {
        }

        void Tetrahedron::create_regular_entities ( )
        {
            // edges
            for ( int i = 0; i < 6; ++i )
            {
                add_regular_entity ( 1, std::vector<int>( &EDGES_OF_TET[i][0], &EDGES_OF_TET[i][2] ) );
            }

            // faces
            for ( int i = 0; i < 4; ++i )
            {
                add_regular_entity ( 2, std::vector<int>( &FACES_OF_TET[i][0], &FACES_OF_TET[i][3] ) );
            }
        }

        void Tetrahedron::create_refined_vertices ( )
        {
            // mid-edge vertices
            for ( int i = 0; i < 6; ++i )
            {
                add_refined_vertex ( std::vector<int>( &EDGES_OF_TET[i][0], &EDGES_OF_TET[i][2] ) );
            }

            // mid-face vertices
            for ( int i = 0; i < 4; ++i )
            {
                add_refined_vertex ( std::vector<int>( &FACES_OF_TET[i][0], &FACES_OF_TET[i][3] ) );
            }

            // mid vertex
            add_refined_vertex ( std::vector<int>( &MIDPOINT_OF_TET[0], &MIDPOINT_OF_TET[4] ) );
        }

        void Tetrahedron::create_refined_cells ( )
        {
            // tet -> 8 tet subcells
            for ( int i = 0; i < 8; ++i )
            {
                add_refined_cell ( std::vector<int>( &TET_8TET_SUB_CELLS[i][0], &TET_8TET_SUB_CELLS[i][4] ) );
            }

            // tet -> 4 hex subcells
            for ( int i = 0; i < 4; ++i )
            {
                add_refined_cell ( std::vector<int>( &TET_4HEX_SUB_CELLS[i][0], &TET_4HEX_SUB_CELLS[i][8] ) );
            }
        }

        void Tetrahedron::create_refinements ( )
        {
            add_refinement ( std::vector<int>( &TET_8TET_REFINEMENT[0], &TET_8TET_REFINEMENT[8] ) );
            add_refinement ( std::vector<int>( &TET_4HEX_REFINEMENT[0], &TET_4HEX_REFINEMENT[4] ) );
        }

        bool Tetrahedron::check_cell_geometry ( const std::vector<Coordinate>& coords, GDim gdim ) const
        {
            // TODO Can we do this in higher geom. dimensions?
            assert ( gdim == 3 );

            Vec<3, double> pt[4];
            int k = 0;
            for ( int i = 0; i < 4; ++i )
            {
                for ( int j = 0; j < 3; ++j )
                {
                    pt[i][j] = coords[k++];
                }
            }

            Vec<3, double> v[3];
            for ( int i = 0; i < 3; ++i )
            {
                v[i] = pt[i + 1] - pt[0];
            }

            const double orientation = dot ( v[2], cross ( v[0], v[1] ) );

            return orientation > 0.0;
        }

        //////////////// HEXAHEDRON ////////////////
        namespace
        {
            const int EDGES_OF_HEX[12][2] = {
                {0, 1 },
                {1, 2 },
                {2, 3 },
                {3, 0 },
                {4, 5 },
                {5, 6 },
                {6, 7 },
                {7, 4 },
                {0, 4 },
                {1, 5 },
                {2, 6 },
                {3, 7 }
            };

            const int FACES_OF_HEX[6][4] = {
                {0, 3, 2, 1 }, // bottom 0
                {0, 1, 5, 4 }, // front  1
                {0, 4, 7, 3 }, // left   2
                {1, 2, 6, 5 }, // right  3
                {2, 3, 7, 6 }, // back   4
                {4, 5, 6, 7 } // top    5
            };

            const int MID_CELL_VERTEX[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

            const int HEX_8HEX_SUB_CELLS[8][8] = {
                {0, 8, 20, 11, 16, 21, 26, 22 },
                {8, 1, 9, 20, 21, 17, 23, 26 },
                {20, 9, 2, 10, 26, 23, 18, 24 },
                {11, 20, 10, 3, 22, 26, 24, 19 },
                {16, 21, 26, 22, 4, 12, 25, 15 },
                {21, 17, 23, 26, 12, 5, 13, 25 },
                {26, 23, 18, 24, 25, 13, 6, 14 },
                {22, 26, 24, 19, 15, 25, 14, 7 }
            };

            const int HEX_8HEX_REFINEMENT[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
        }

        Hexahedron::Hexahedron ( )
        : CellType ( CellType::HEXAHEDRON, 3, 8 )
        {
        }

        void Hexahedron::create_regular_entities ( )
        {
            // edges
            for ( int i = 0; i < 12; ++i )
            {
                add_regular_entity ( 1, std::vector<int>( &EDGES_OF_HEX[i][0], &EDGES_OF_HEX[i][2] ) );
            }

            // faces
            for ( int i = 0; i < 6; ++i )
            {
                add_regular_entity ( 2, std::vector<int>( &FACES_OF_HEX[i][0], &FACES_OF_HEX[i][4] ) );
            }
        }

        void Hexahedron::create_refined_vertices ( )
        {
            // mid-edge vertices
            for ( int i = 0; i < 12; ++i )
            {
                add_refined_vertex ( std::vector<int>( &EDGES_OF_HEX[i][0], &EDGES_OF_HEX[i][2] ) );
            }

            // mid-face vertices
            for ( int i = 0; i < 6; ++i )
            {
                add_refined_vertex ( std::vector<int>( &FACES_OF_HEX[i][0], &FACES_OF_HEX[i][4] ) );
            }

            // cell mid-point
            add_refined_vertex ( std::vector<int>( &MID_CELL_VERTEX[0], &MID_CELL_VERTEX[8] ) );
        }

        void Hexahedron::create_refined_cells ( )
        {
            // hex -> 8hex sub-cells
            for ( int i = 0; i < 8; ++i )
            {
                add_refined_cell ( std::vector<int>( &HEX_8HEX_SUB_CELLS[i][0], &HEX_8HEX_SUB_CELLS[i][8] ) );
            }
        }

        void Hexahedron::create_refinements ( )
        {
            add_refinement ( std::vector<int>( &HEX_8HEX_REFINEMENT[0], &HEX_8HEX_REFINEMENT[8] ) );
        }

        //////////////// PYRAMID ////////////////
        namespace
        {
            const int EDGES_OF_PYR[8][2] = {
                {0, 1 },
                {1, 2 },
                {2, 3 },
                {3, 0 },
                {0, 4 },
                {1, 4 },
                {2, 4 },
                {3, 4 }
            };

            const int FACES_QUAD_OF_PYR[4] = { 0, 3, 2, 1 }; //bottom 0
            const int FACES_TET_OF_PYR[4][3] = {
                {0, 1, 4 }, //front 1
                {1, 2, 4 }, //right 2
                {2, 3, 4 }, //back 3
                {0, 4, 3 } //left 4
            };

            const int MID_QUAD_IN_PYR_VERTEX[4] = { 0, 1, 2, 3 };

            const int PYR_6PYR_SUB_CELLS[6][5] = {
                {9, 10, 11, 12, 4 }, //top
                {9, 12, 11, 10, 13 }, //center
                {0, 5, 13, 8, 9 }, //bottom-front-left
                {1, 6, 13, 5, 10 }, //bottom-front-rigth
                {2, 7, 13, 6, 11 }, //bottom-back-right
                {3, 8, 13, 7, 12 } //bottom-back-left
            };

            const int PYR_4TET_SUB_CELLS[4][4] = {
                {9, 13, 10, 5 }, //front
                {10, 13, 11, 6 }, //right
                {11, 13, 12, 7 }, //back
                {12, 13, 9, 8 } //left
            };

            const int PYR_6PYR4TET_REFINEMENT[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

        }

        Pyramid::Pyramid ( )
        : CellType ( CellType::PYRAMID, 3, 5 )
        {
        }

        void Pyramid::create_regular_entities ( )
        {
            // edges
            for ( int i = 0; i < 8; ++i )
            {
                add_regular_entity ( 1, std::vector<int>( &EDGES_OF_PYR[i][0], &EDGES_OF_PYR[i][2] ) );
            }

            // faces
            add_regular_entity ( 2, std::vector<int>( &FACES_QUAD_OF_PYR[0], &FACES_QUAD_OF_PYR[4] ) );
            for ( int i = 0; i < 4; ++i )
            {
                add_regular_entity ( 2, std::vector<int>( &FACES_TET_OF_PYR[i][0], &FACES_TET_OF_PYR[i][3] ) );
            }
        }

        void Pyramid::create_refined_vertices ( )
        {
            // mid-edge vertices
            for ( int i = 0; i < 8; ++i )
            {
                add_refined_vertex ( std::vector<int>( &EDGES_OF_PYR[i][0], &EDGES_OF_PYR[i][2] ) );
            }

            // mid-face vertices
            add_refined_vertex ( std::vector<int>( &MID_QUAD_IN_PYR_VERTEX[0], &MID_QUAD_IN_PYR_VERTEX[4] ) );
        }

        void Pyramid::create_refined_cells ( )
        {
            // 6 pyramids
            for ( int i = 0; i < 6; ++i )
            {
                add_refined_cell ( std::vector<int>( &PYR_6PYR_SUB_CELLS[i][0], &PYR_6PYR_SUB_CELLS[i][5] ) );
            }

            // 4 tetrahedrons
            for ( int i = 0; i < 4; ++i )
            {
                add_refined_cell ( std::vector<int>( &PYR_4TET_SUB_CELLS[i][0], &PYR_4TET_SUB_CELLS[i][4] ) );
            }
        }

        void Pyramid::create_refinements ( )
        {
            add_refinement ( std::vector<int>( &PYR_6PYR4TET_REFINEMENT[0], &PYR_6PYR4TET_REFINEMENT[10] ) );
        }

    } // namespace mesh
} // namespace hiflow
