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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-10-07

#ifndef SRC_LINEAR_SOLVER_AMG_STANDARDCOARSENINGSETTINGS_H_
#    define SRC_LINEAR_SOLVER_AMG_STANDARDCOARSENINGSETTINGS_H_

/// Settings of the StandardCoarseningMethod

struct StandardCoarseningSettings
{
    /// Default constructor

    StandardCoarseningSettings ( )
    : strength_threshold ( 0.25 )
    {
    }

    /// Threshold defining the strongly coupled variables.
    /// Typically in the range of [0.25, 0.5].
    double strength_threshold;
};

#endif /* SRC_LINEAR_SOLVER_AMG_STANDARDCOARSENINGSETTINGS_H_ */
