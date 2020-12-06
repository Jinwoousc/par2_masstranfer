/**
* @file PParticles.cu
* @brief MoveParticle functor.
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

#ifndef PAR2_MOVEPARTICLE_CUH
#define PAR2_MOVEPARTICLE_CUH

#include <thrust/tuple.h>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>
#include <curand_kernel.h>
#include "../Geometry/CartesianGrid.cuh"
#include "../Geometry/FaceField.cuh"
#include "../Geometry/CornerField.cuh"
#include "../Geometry/MassField.cuh"
#include "../Geometry/Vector.cuh"

namespace par2
{
    /**
    * @struct MoveParticle
    * @brief Thrust functor for one step of the particle tracking method.
    * @tparam T Float number precision
    */
    template<typename T>
    struct MoveParticle
    {
        // Raw pointer to velocity vectors (facefield)
        T* datax;
        T* datay;
        T* dataz;
        // Raw pointer to velocity vectors (cornerfield or cellfield)
        T* cdatax;
        T* cdatay;
        T* cdataz;
        // Raw pointer to mass field
        int* mmb1;
        int* mmb2;
        int* mimb1;
        int* mimb2;
        // Grid and physical variables
        grid::Grid<T> grid;
        T dt;
        unsigned int nParticles;
        int massType;
        T mass1alpha;
        T mass1beta;
        T mass2alpha;
        T mass2beta;
        T expvalue1;
        T expvalue2;
        T molecularDiffusion;
        T alphaL;
        T alphaT;
        // Mass transfer transition probability
        T P100, P101, P110, P111, P200, P201, P210, P211;
        // Raw pointer to curand state vector
        curandState_t* states;
        // If true, use trilinear interpolation
        bool useTrilinearCorrection;

        /**
        * @brief Constructor.
        * @param _grid Grid where the velocity is defined
        */
        MoveParticle(const grid::Grid<T>& _grid) : grid(_grid), dt(0.0) {};

        /**
        * @brief Initialize the functor.
        * @param _datax Velocity vector (x-direction)
        * @param _datay Velocity vector (y-direction)
        * @param _dataz Velocity vector (z-direction)
        * @param _molecularDiffusion Effective molecular diffusion
        * @param _alphaL Longitudinal dispersivity
        * @param _alphaT Transverse dispersivity
        * @param _states Curand states
        * @param _useTrilinearCorrection True if trilinear correction is used
        */
        void initialize(thrust::device_vector<T> &_datax,
                        thrust::device_vector<T> &_datay,
                        thrust::device_vector<T> &_dataz,
                        thrust::device_vector<T> &_cdatax,
                        thrust::device_vector<T> &_cdatay,
                        thrust::device_vector<T> &_cdataz,
                        thrust::device_vector<int> &_mmb1,
                        thrust::device_vector<int> &_mmb2,
                        thrust::device_vector<int> &_mimb1,
                        thrust::device_vector<int> &_mimb2,
                        T _nParticles,
                        T _molecularDiffusion,
                        T _alphaL,
                        T _alphaT,
                        int _massType,
                        T _mass1alpha,
                        T _mass1beta,
                        T _mass2alpha,
                        T _mass2beta,
                        thrust::device_vector<curandState_t> &_states,
                        bool _useTrilinearCorrection)
        {
            datax = thrust::raw_pointer_cast(_datax.data());
            datay = thrust::raw_pointer_cast(_datay.data());
            dataz = thrust::raw_pointer_cast(_dataz.data());

            cdatax = thrust::raw_pointer_cast(_cdatax.data());
            cdatay = thrust::raw_pointer_cast(_cdatay.data());
            cdataz = thrust::raw_pointer_cast(_cdataz.data());

            mmb1 = thrust::raw_pointer_cast(_mmb1.data());
            mmb2 = thrust::raw_pointer_cast(_mmb2.data());
            mimb1 = thrust::raw_pointer_cast(_mimb1.data());
            mimb2 = thrust::raw_pointer_cast(_mimb2.data());

            nParticles = _nParticles;
            molecularDiffusion = _molecularDiffusion;
            alphaL = _alphaL;
            alphaT = _alphaT;
            massType = _massType;
            mass1alpha =  _mass1alpha;
            mass1beta = _mass1beta;
            mass2alpha =  _mass2alpha;
            mass2beta = _mass2beta;

            states = thrust::raw_pointer_cast(_states.data());

            useTrilinearCorrection = _useTrilinearCorrection;
        }

        /**
        * @brief Set time step.
        * @param _dt Time step
        */
        void setTimeStep(T _dt)
        {
            dt = _dt;

            // mass transition matrix
            if (massType == 1)
            {
                expvalue1 = exp(-(1+mass1beta)*mass1alpha*dt);
                P100 = ( 1 + mass1beta*expvalue1 ) / ( 1+mass1beta );
                P101 = ( 1 - expvalue1 ) / ( 1+mass1beta );
                P110 = ( mass1beta - mass1beta*expvalue1 ) / ( 1+mass1beta );
                P111 = ( mass1beta + expvalue1 ) / ( 1+mass1beta );
            }
            else if (massType == 2)
            {
                expvalue1 = exp(-(1+mass1beta)*mass1alpha*dt);
                P100 = ( 1 + mass1beta*expvalue1 ) / ( 1+mass1beta );
                P101 = ( 1 - expvalue1 ) / ( 1+mass1beta );
                P110 = ( mass1beta - mass1beta*expvalue1 ) / ( 1+mass1beta );
                P111 = ( mass1beta + expvalue1 ) / ( 1+mass1beta );
                expvalue2 = exp(-(1+mass2beta)*mass2alpha*dt);
                P200 = ( 1 + mass2beta*expvalue2 ) / ( 1+mass2beta );
                P201 = ( 1 - expvalue2 ) / ( 1+mass2beta );
                P210 = ( mass2beta - mass2beta*expvalue2 ) / ( 1+mass2beta );
                P211 = ( mass2beta + expvalue2 ) / ( 1+mass2beta );
            }
        }

        using Position = thrust::tuple<T, T, T, T, T, unsigned int>;

        /**
        * @brief Execute one step of the particle tracking method
        *        on one particle.
        * @param p Initial position of the particle
        * @return Final position of the particle
        */
        __device__
        Position operator()(Position p) const
        {
            int idx, idy, idz;
            grid::idPoint(grid,
                          thrust::get<0>(p),
                          thrust::get<1>(p),
                          thrust::get<2>(p),
                          &idx, &idy, &idz);

            int id = grid::mergeId(grid, idx, idy, idz);

            bool idValid = grid::validId<T>(grid, idx, idy, idz);

            // Mass transfer consideration
            bool mobilephase = true;
            if (massType == 1)
            {
                // calculate next time step mass distribution & mass transition probability
                // mmb; mimb;
                // heterogeneous masstransfer
                T nmmb = P100*mmb2[id] + P101*mimb2[id];
                T nmimb = P110*mmb2[id] + P111*mimb2[id];
                // homogeneous masstransfer
                // T nmmb = P100*mmb2[0] + P101*mimb2[0];
                // T nmimb = P110*mmb2[0] + P111*mimb2[0];

                T pmmb = nmmb/(nmmb+nmimb);
                thrust::get<4>(p) = pmmb;

                // apply mass transition probability to particles
                T randomY = curand_uniform_double(&states[thrust::get<5>(p)]);

                if (thrust::get<3>(p) == 0) // from mobile phase
                {
                    if (randomY <= pmmb) // to mobile phase
                    {
                        // heterogeneous masstransfer
                        atomicAdd( &mmb1[id], -1 );
                        // homogeneous masstransfer
                        // atomicAdd( &mmb1[0], -1 );
                    }
                    else // to immobile phase
                    {
                        mobilephase = false;
                        // heterogeneous masstransfer
                        atomicAdd( &mmb1[id], -1 );
                        atomicAdd( &mimb1[id], 1 );
                        // homogeneous masstransfer
                        // atomicAdd( &mmb1[0], -1 );
                        // atomicAdd( &mimb1[0], 1 );
                        thrust::get<3>(p) = 1;
                    }
                }
                else if (thrust::get<3>(p) == 1) // from immobile phase
                {
                    if (randomY <= pmmb) // to mobile phase
                    {
                        // heterogeneous masstransfer
                        atomicAdd( &mimb1[id], -1 );
                        // homogeneous masstransfer
                        // atomicAdd( &mimb1[0], -1 );
                        thrust::get<3>(p) = 0;
                    }
                    else // to immobile phase
                    {
                        mobilephase = false;
                    }
                }
            }
            else if (massType == 2)
            {
                // calculate next time step mass distribution & mass transition probability
                // mmb; mimb;

                // apply mass transition probability to particles
                T randomY = curand_uniform_double(&states[thrust::get<5>(p)]);
                mobilephase = true;
                // heterogeneous masstransfer
                atomicAdd( &mmb1[id], -1 );
                // homogeneous masstransfer
                // atomicAdd( &mmb1[0], -1 );
            }
            else
            {
                // heterogeneous masstransfer
                atomicAdd( &mmb1[id], -1 );
                // homogeneous masstransfer
                // atomicAdd( &mmb1[0], -1 );
            }

            if (mobilephase)
            {
                // Velocity term (linear interpolation)
                T vlx, vly, vlz;
                facefield::in<T>
                            (datax, datay, dataz, grid, idx, idy, idz,
                             idValid,
                             thrust::get<0>(p),
                             thrust::get<1>(p),
                             thrust::get<2>(p),
                             &vlx, &vly, &vlz);

                T vcx, vcy, vcz;
                if (useTrilinearCorrection)
                {
                    // Velocity correction div(D) (trilinear interpolation)
                    cornerfield::velocityCorrection<T>
                                (cdatax, cdatay, cdataz, grid, idx, idy, idz,
                                 idValid,
                                 thrust::get<0>(p),
                                 thrust::get<1>(p),
                                 thrust::get<2>(p),
                                 molecularDiffusion, alphaL, alphaT,
                                 &vcx, &vcy, &vcz);
                }
                else
                {
                    // Velocity correction div(D) (block-centered finite difference)
                    vcx = idValid ? cdatax[id] : 0;
                    vcy = idValid ? cdatay[id] : 0;
                    vcz = idValid ? cdataz[id] : 0;
                }

                // Displacement Matrix (trilinear interpolation)
                T B00, B11, B22, B01, B02, B12;
                cornerfield::displacementMatrix<T>
                            (cdatax, cdatay, cdataz, grid, idx, idy, idz,
                             idValid,
                             thrust::get<0>(p),
                             thrust::get<1>(p),
                             thrust::get<2>(p),
                             molecularDiffusion, alphaL, alphaT, dt,
                             &B00, &B11, &B22, &B01, &B02, &B12);

                // Random Displacement
                T xi0 = curand_normal_double(&states[thrust::get<5>(p)]);
                T xi1 = curand_normal_double(&states[thrust::get<5>(p)]);
                T xi2 = curand_normal_double(&states[thrust::get<5>(p)]);

                // Update positions
                T dpx, dpy, dpz;

                dpx = (idValid ? ((vlx + vcx)*dt + (B00*xi0 + B01*xi1 + B02*xi2)) : 0);
                dpy = (idValid ? ((vly + vcy)*dt + (B01*xi0 + B11*xi1 + B12*xi2)) : 0);
                if (grid.nz != 1)
                    dpz = (idValid ? ((vlz + vcz)*dt + (B02*xi0 + B12*xi1 + B22*xi2)) : 0);

                // Closed condition: particles cannot exit the domain
                thrust::get<0>(p) += (grid::validX(grid, thrust::get<0>(p) + dpx) ? dpx : 0);
                thrust::get<1>(p) += (grid::validY(grid, thrust::get<1>(p) + dpy) ? dpy : 0);
                if (grid.nz != 1)
                    thrust::get<2>(p) += (grid::validZ(grid, thrust::get<2>(p) + dpz) ? dpz : 0);

                // Update mobile mass field
                // heterogeneous masstransfer
                grid::idPoint(grid,
                              thrust::get<0>(p),
                              thrust::get<1>(p),
                              thrust::get<2>(p),
                              &idx, &idy, &idz);
                int id = grid::mergeId(grid, idx, idy, idz);
                atomicAdd( &mmb1[id], 1 );

                // homogeneous masstransfer
                // atomicAdd( &mmb1[0], 1 );
            }
            return p;
        }
    };
}

#endif //PAR2_MOVEPARTICLE_CUH
