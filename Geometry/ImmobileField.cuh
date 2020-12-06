/**
* @file ImmobileField.cuh
* @brief Header file for cellfield.
*        A cellfield is a field that is defined at the center of
*        every cells of the grid.
*
* @author Calogero B. Rizzo
*
* @copyright This file is part of the PAR2 software.
*            Copyright (C) 2018 Calogero B. Rizzo
*
* @license This program is free software: you can redistribute it and/or modify
*          it under the terms of the GNU General Public License as published by
*          the Free Software Foundation, either version 3 of the License, or
*          (at your option) any later version.
*
*          This program is distributed in the hope that it will be useful,
*          but WITHOUT ANY WARRANTY; without even the implied warranty of
*          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*          GNU General Public License for more details.
*
*          You should have received a copy of the GNU General Public License
*          along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PAR2_IMMOBILEFIELD_CUH
#define PAR2_IMMOBILEFIELD_CUH

#include <fstream>
#include <iostream>
#include <algorithm>

#include "Vector.cuh"
#include "CartesianGrid.cuh"

namespace par2
{

    namespace immobillefield
    {
        /**
        * @brief Build the vector containing a cellfield.
        * @param g Grid where the field is defined
        * @param data Vector that contains the cellfield values
        * @param v Init value for data
        * @tparam T Float number precision
        * @tparam Vector Container for data vectors
        */
        template<typename T, class Vector>
        void build(const grid::Grid<T>& g, Vector& data, T v = 0)
        {
            int size = g.nx*g.ny*g.nz;
            data.resize(size, v);
        }

    };
};

#endif //PAR2_CELLFIELD_CUH
