/**
* @file MassField.cuh
* @brief Header file for massfield.
*        A massfield is a field that is defined at the center of
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

#ifndef PAR2_MASSFIELD_CUH
#define PAR2_MASSFIELD_CUH

#include <fstream>
#include <iostream>
#include <algorithm>

#include "Vector.cuh"
#include "CartesianGrid.cuh"

namespace par2
{

    namespace massfield
    {
        /**
        * @brief Build the vector containing a massfield.
        * @param g Grid where the field is defined
        * @param data Vector that contains the massfield values
        * @param m Init value for data
        * @tparam T Float number precision
        * @tparam Vector Container for data vectors
        */
        template<typename T, class Vector>
        void build(const grid::Grid<T>& g, Vector& mmb1, Vector& mmb2, Vector& mimb1, Vector& mimb2, T m = 0)
        {
            int size = g.nx*g.ny*g.nz;
            mmb1.resize(size, m);
            mmb2.resize(size, m);
            mimb1.resize(size, m);
            mimb2.resize(size, m);
        }
    };
};

#endif //PAR2_MASSFIELD_CUH
