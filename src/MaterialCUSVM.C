// -*-c++-*-
//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// # Produced at the Lawrence Livermore National Laboratory.
// #
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// #
// # LLNL-CODE-643337
// #
// # All rights reserved.
// #
// # This file is part of SW4, Version: 1.0
// #
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// #
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991.
// #
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details.
// #
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

#include "Require.h"

#include <cstring>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <fcntl.h>
#include "EW.h"
#include "MaterialCUSVM.h"
#include "Byteswapper.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif

using namespace std;


//-----------------------------------------------------------------------
MaterialCUSVM::MaterialCUSVM( EW* a_ew, const string a_file, const string a_directory):
    mEW(a_ew),
    m_model_file(a_file),
    m_model_dir(a_directory),
    m_use_attenuation(false)
{
    mCoversAllPoints = false;
    // Check that the depths make sense
    if (a_ew != NULL) {
        m_use_attenuation = a_ew->usingAttenuation();
        read_CUSVM();
    }
}

//-----------------------------------------------------------------------
MaterialCUSVM::~MaterialCUSVM()
{
}

//-----------------------------------------------------------------------
void MaterialCUSVM::set_material_properties(std::vector<Sarray> & rho,
        std::vector<Sarray> & cs,
        std::vector<Sarray> & cp,
        std::vector<Sarray> & xis,
        std::vector<Sarray> & xip )
{
// Assume attenuation arrays defined on all grids if they are defined on grid zero.
    bool use_q = m_use_attenuation && xis[0].is_defined() && xip[0].is_defined();
    bool is_debug = true;

    bool bulldoze = false;
    bool even_stretch = false;

    uint64_t outside=0, material=0;
    double lon, lat, depth, topo_elev, min_depth;
    char inname[128], outname[128], cmd[2048];

    FILE *fptr;

    double mylon, mylat, myz, vp, vs, density, vp_min, vs_min, density_min;
    int containunit;

    // Find the relative dimension size of upper and lower interface for each grid patch
    int i1, j1;
    int g_fac[16];

    // No grid size reduction at the first curvilinear grid
    g_fac[mEW->mNumberOfGrids-1] = 1;
    for(int g=mEW->mNumberOfGrids-2; g >= 0; g--) {
        if(g != mEW->mNumberOfCartesianGrids-1)
            g_fac[g] = g_fac[g+1] * 2;
        else
            g_fac[g] = g_fac[g+1];
    }

    // Get topo
    /* fprintf(stderr, "Start to query topo\n"); */
    /* mEW->extractTopographyFromCUSVM("cusvm"); */
    /* fprintf(stderr, "Done query topo\n"); */

    m_zminloc = 0;
    vp_min = 200;
    vs_min = 0;
    density_min = 800;
    // Note, due to the abnormal surface velocity, we replace it with depth=10 
    min_depth = 10;

    for(int g=0; g < mEW->mNumberOfGrids; g++) {
        sprintf(inname, "/tmp/cusvm.in.%d.%d", g, mEW->getRank());
        sprintf(outname, "/tmp/cusvm.out.%d.%d", g, mEW->getRank());
        sprintf(cmd, "mkdir -p %s", outname);
        system(cmd);

        bool curvilinear = mEW->topographyExists() && g >= mEW->mNumberOfCartesianGrids;
        int npts = (mEW->m_jEnd[g] - mEW->m_jStart[g] + 1) * (mEW->m_iEnd[g] - mEW->m_iStart[g] + 1) * (mEW->m_kEnd[g] -  mEW->m_kStart[g] + 1);

        // Query Vp, Vs, Rho
        fptr = fopen(inname, "w");
        fprintf(fptr, "%d\n", npts);
        for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i) {
            for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j) {
                // get lat lon for current grid g
                float_sw4 x = (i-1)*mEW->mGridSize[g];
                float_sw4 y = (j-1)*mEW->mGridSize[g];
                float_sw4 z = 0;
                // (x, y, z) is the coordinate of current grid point
                mEW->computeGeographicCoord(x, y, lon, lat);

                /* i1 = i * g_fac[g]; */
                /* j1 = j * g_fac[g]; */
                /* topo_elev = mEW->mTopo(i1, j1, 1); */

                for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k) {
                    if( curvilinear )
                        z = mEW->mZ[g](i,j,k);
                    else
                        z = mEW->m_zmin[g] + (k-1)*mEW->mGridSize[g];

                    depth = z;
                    if (depth < min_depth)
                        depth = min_depth;

                    fprintf(fptr, "%f %f %f\n", lon, lat, depth);
                    material++;
                } // End for k
            } // End for i
        } // End for j
        fclose(fptr);

        // query CUSVM
        printf("Rank %d: query CUSVM with %d pts\n", mEW->getRank(), npts);
        sprintf(cmd, "/pscratch/sd/h/houhun/cusvm.tang/geodataquery 0 2 /pscratch/sd/h/houhun/cusvm.tang/DATABASE %s %s", inname, outname);
        /* sprintf(cmd, "/pscratch/sd/h/houhun/cusvm.tang/geodataquery 0 1 /pscratch/sd/h/houhun/cusvm.tang/DATABASE %s %s", inname, outname); */
        system(cmd);

        if (is_debug)
            fprintf(stderr, "Rank %d: Grid %d, queried %d points\n", mEW->getRank(), g, material);

        material = 0;
        // read the output
        sprintf(outname, "/tmp/cusvm.out.%d.%d/output.out", g, mEW->getRank());
        fptr = fopen(outname, "r");
        int read_npts;
        fscanf(fptr, "%d", &read_npts);

        for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i) {
            for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j) {
                float_sw4 x = (i-1)*mEW->mGridSize[g];
                float_sw4 y = (j-1)*mEW->mGridSize[g];
                float_sw4 z;

                i1 = i * g_fac[g];
                j1 = j * g_fac[g];
                /* topo_elev = mEW->mTopo(i1, j1, 1); */

                mEW->computeGeographicCoord(x, y, lon, lat);

                for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k) {
                    if( curvilinear )
                        z = mEW->mZ[g](i,j,k);
                    else
                        z = mEW->m_zmin[g] + (k-1)*mEW->mGridSize[g];

                    // There is a value discontinuity issue with the USGS data due to stiching different models
                    // current workaround is using depth=10 instead of actual surface values
                    depth = z;
                    if (depth < min_depth)
                        depth = min_depth;

                    /* if (i == 1200 && k == 1) { */
                    /*     fprintf(stderr, "Rank %d grid %d: (%d, %d, %d) [%lf, %lf, %lf] : %lf, %lf, %lf\n", */
                    /*                      mEW->getRank(), g, i, j, k, mylon, mylat, myz, vp, vs, density); */
                    /* } */

                    material++;
                    fscanf(fptr, "%lf %lf %lf %lf %lf %lf %d", &mylon, &mylat, &myz, &vp, &vs, &density, &containunit);

                    if (vp < vp_min) {
                        /* fprintf(stderr, "Rank %d grid %d: (%d, %d, %d) [%f, %f, %f] vp = %lf, adjust to %lf\n", */
                        /*                  mEW->getRank(), g, i, j, k, mylon, mylat, myz, vp, vp_min); */
                        vp = vp_min;
                    }
                    if (vs < vs_min) {
                        fprintf(stderr, "Rank %d grid %d: (%d, %d, %d) [%f, %f, %f] vs = %lf, adjust to %lf\n",
                                         mEW->getRank(), g, i, j, k, mylon, mylat, myz, vs, vs_min);
                        vs = vs_min;
                    }
                    if (density < density_min) {
                        fprintf(stderr, "Rank %d grid %d: (%d, %d, %d) [%f, %f, %f] rho = %lf, adjust to %lf\n",
                                         mEW->getRank(), g, i, j, k, mylon, mylat, myz, density, density_min);
                        density = density_min;
                    }

                    rho[g](i, j, k) = density;
                    cp[g](i, j, k)  = vp;
                    cs[g](i, j, k)  = vs;
                    if( use_q ) {
                        // Qs=36               for VS < 360;
                        // Qs=0.1*Vs           for 360 ≤ VS < 500;
                        // Qs=0.3* VS − 100    for 500 ≤ VS < 1000;
                        // Qs=167 ·*VS+  33    for 1000 ≤ VS < 4000
                        // Qs=700              for VS ≥ 4000
                        // Qp=2*Qs
                        if (vs <= 360)
                            xis[g](i, j, k)  = 36;
                        else if (vs > 360 && vs < 500)
                            xis[g](i, j, k)  = vs * 0.1;
                        else if (vs >= 500 && vs < 1000)
                            xis[g](i, j, k)  = vs * 0.3 - 100.0;
                        else if (vs >= 1000 && vs < 4000)
                            xis[g](i, j, k)  = vs * 167.0 + 33.0;
                        else if (vs >= 4000)
                            xis[g](i, j, k)  = 700;

                        xip[g](i, j, k)  = xis[g](i, j, k) * 2.0;
                    }


                    if (fabs(lon - mylon) > 1e-5 )
                        printf("x=%.1ff, y=%.1f, sw4_lon=%f does not match ucvm_lon=%f!\n", x, y, lon, mylon);
                    if (fabs(lat - mylat) > 1e-5 )
                        printf("x=%.1f, y=%.1f, sw4_lat=%f does not match ucvm_lat=%f!\n", x, y, lat, mylat);
                    if (fabs(depth - myz) > 1e-5 )
                        printf("x=%.1f, y=%.1f, sw4_z=%f does not match ucvm_z=%f!\n", x, y, depth, myz);

                    if (density < 10 || density > 10000) {
                        fprintf(stderr, "Rank %d grid %d: (%d, %d, %d) [%f, %f, %f] density = %f\n",
                                         mEW->getRank(), g, i, j, k, mylon, mylat, myz, density);
                    }
                    if (vp < 10 || vp > 10000) {
                        fprintf(stderr, "Rank %d grid %d: (%d, %d, %d) [%f, %f, %f] vp = %f\n",
                                         mEW->getRank(), g, i, j, k, mylon, mylat, myz, vp);
                    }
                    if (vs < 10 || vs > 10000) {
                        fprintf(stderr, "Rank %d grid %d: (%d, %d, %d) [%f, %f, %f] vs = %f\n",
                                         mEW->getRank(), g, i, j, k, mylon, mylat, myz, vs);
                    }

                } // End for k
            } // End for j
        } // End for i
        fclose(fptr);
        if (is_debug) {
            printf("Read %d points\n", material);
            printf("Rank %d grid %d: rho min = %.2f, max = %.2f\n", mEW->getRank(), g, rho[g].minimum(), rho[g].maximum());
            printf("Rank %d grid %d: cp min = %.2f, max = %.2f\n", mEW->getRank(), g, cp[g].minimum(), cp[g].maximum());
            printf("Rank %d grid %d: cs min = %.2f, max = %.2f\n", mEW->getRank(), g, cs[g].minimum(), cs[g].maximum());
        }

    } // end for g...

    mEW->communicate_arrays( rho );
    mEW->communicate_arrays( cs );
    mEW->communicate_arrays( cp );
    mEW->material_ic( rho );
    mEW->material_ic( cs );
    mEW->material_ic( cp );
    if( use_q ) {
        mEW->communicate_arrays( xis );
        mEW->communicate_arrays( xip );
        mEW->material_ic( xis);
        mEW->material_ic( xip );
    }

    size_t materialSum, outsideSum;
    int mpisizelong, mpisizelonglong, mpisizeint;
    MPI_Type_size(MPI_LONG,&mpisizelong );
    MPI_Type_size(MPI_LONG_LONG,&mpisizelonglong );
    MPI_Type_size(MPI_INT,&mpisizeint );
    if( sizeof(size_t) == mpisizelong ) {
        MPI_Reduce(&material, &materialSum, 1, MPI_LONG, MPI_SUM, 0, mEW->m_1d_communicator );
        MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG, MPI_SUM, 0, mEW->m_1d_communicator );
    }
    else if( sizeof(size_t) == mpisizelonglong ) {
        MPI_Reduce(&material, &materialSum, 1, MPI_LONG_LONG, MPI_SUM, 0, mEW->m_1d_communicator );
        MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG_LONG, MPI_SUM, 0, mEW->m_1d_communicator );
    }
    else if( sizeof(size_t) == mpisizeint ) {
        MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0, mEW->m_1d_communicator );
        MPI_Reduce(&outside,   &outsideSum, 1, MPI_INT, MPI_SUM, 0, mEW->m_1d_communicator );
    }
    else {
        int materialsumi, outsidesumi, materiali=material, outsidei=outside;
        MPI_Reduce(&materiali, &materialsumi, 1, MPI_INT, MPI_SUM, 0, mEW->m_1d_communicator );
        MPI_Reduce(&outsidei,   &outsidesumi, 1, MPI_INT, MPI_SUM, 0, mEW->m_1d_communicator );
        materialSum=materialsumi;
        outsideSum=outsidesumi;
    }
    if (mEW->getRank() == 0)
        //      cout << endl
        //           << "--------------------------------------------------------------\n"
        //           << "CUSVM Initialized Node Types: " << endl
        //           << "   Material:        " << materialSum << endl
        //           << endl
        //           << "*Outside Domain:    " << outsideSum << endl
        //           << endl
        //           << "--------------------------------------------------------------\n"
        //           << endl;
        cout << endl
             << "CUSVM command: outside = " << outsideSum << ", material = " << materialSum << endl;

    /* material_check(false); */
} // End of set_material_properties


//-----------------------------------------------------------------------
void MaterialCUSVM::read_CUSVM()
{
    // Timers
    double time_start, time_end;
    double intf_start, intf_end, mat_start, mat_end;
    time_start = MPI_Wtime();

    /* fill_in_fluids(); */

    time_end = MPI_Wtime();
    if (mEW->getRank() == 0) {
        cout << "MaterialCUSVM::read_CUSVM, time to read material file: " << time_end - time_start << " seconds." << endl;
    }
    cout.flush();
}

//-----------------------------------------------------------------------
void MaterialCUSVM::fill_in_fluids()
{
// Start from p=0
// start from the last (bottom) block and progress upwards
    if( !m_outside ) {
        for( int p=m_npatches-1 ; p >= 0; p-- ) {
            if( !m_isempty[p] ) {
                #pragma omp parallel for
                for( int j=mMaterial_cs[p].m_jb ; j <= mMaterial_cs[p].m_je ; j++ ) {
                    for( int i=mMaterial_cs[p].m_ib ; i <= mMaterial_cs[p].m_ie ; i++ ) {
                        int k0 = mMaterial_cs[p].m_kb;
                        while( mMaterial_cs[p](1,i,j,k0) < 0 && k0 < mMaterial_cs[p].m_ke )
                            k0++;
                        // consider the case where the top block is all water. Then k0 = mMaterial[p].m_ke and mMaterial[p](3,i,j,k0)=-999
                        // k0 is now the first k with cs > 0.
                        if (mMaterial_cs[p](1,i,j,k0) < 0) {
                            // get value from block p+1
                            if (p<m_npatches-1) {
                                int pd=p+1, id, jd, kd; // index of donor block
                                float_sw4 xm=(i-1)*m_hh[p];
                                float_sw4 ym=(j-1)*m_hh[p];
                                // get closest (id,jd) index on patch pd
                                id = static_cast<int>( 1 + trunc(xm/m_hh[pd]) );
                                jd = static_cast<int>( 1 + trunc(ym/m_hh[pd]) );
                                kd = mMaterial_cs[pd].m_kb; // get value from top of block pd

                                if (! (id >= mMaterial_cs[pd].m_ib && id <= mMaterial_cs[pd].m_ie &&
                                        jd >= mMaterial_cs[pd].m_jb && jd <= mMaterial_cs[pd].m_je )) {
                                    // out of bounds: find nearest interior point
                                    if (id < mMaterial_cs[pd].m_ib) id=mMaterial_cs[pd].m_ib;
                                    if (id > mMaterial_cs[pd].m_ie) id=mMaterial_cs[pd].m_ie;
                                    if (jd < mMaterial_cs[pd].m_jb) jd=mMaterial_cs[pd].m_jb;
                                    if (jd > mMaterial_cs[pd].m_je) jd=mMaterial_cs[pd].m_je;

                                    printf("WARNING: nearest grid point to (%e,%e) was outside local part of block pd=%i\n"
                                           " using id=%i, jd=%i, at (%e, %e)\n", xm, ym, pd, id, jd, (id-1)*m_hh[pd], (jd-1)*m_hh[pd]);

                                }
                                // get values from block 'pd'
                                mMaterial_rho[p](1,i,j,k0)= mMaterial_rho[pd](1,id,jd,kd);
                                mMaterial_cp[p](1,i,j,k0)= mMaterial_cp[pd](1,id,jd,kd);
                                mMaterial_cs[p](1,i,j,k0)= mMaterial_cs[pd](1,id,jd,kd);
                                if (m_use_attenuation) {
                                    mMaterial_qp[p](1,i,j,k0)= mMaterial_qp[pd](1,id,jd,kd);
                                    mMaterial_qs[p](1,i,j,k0)= mMaterial_qs[pd](1,id,jd,kd);
                                }
                            }
                            else {
                                printf("ERROR: found undefined material properties in last material block\n"
                                       " patch p=%i, i=%i, j=%i, k0=%i\n", p, i, j, k0);
                            }
                        }

                        for( int k=mMaterial_cs[p].m_kb ; k < k0 ; k++ ) {
                            mMaterial_rho[p](1,i,j,k) = mMaterial_rho[p](1,i,j,k0);
                            mMaterial_cp[p](1,i,j,k)  = mMaterial_cp[p](1,i,j,k0);
                            mMaterial_cs[p](1,i,j,k)  = mMaterial_cs[p](1,i,j,k0);
                            if( m_use_attenuation ) {
                                mMaterial_qp[p](1,i,j,k) = mMaterial_qp[p](1,i,j,k0);
                                mMaterial_qs[p](1,i,j,k) = mMaterial_qs[p](1,i,j,k0);
                            }
                        } // End for k

                    } // End for i
                } // End for j
            } // End if !m_isempty
        } // End for p
    } // End if !outside
}

//-----------------------------------------------------------------------
void MaterialCUSVM::material_check( bool water )
{
    bool printsmallcpcs=false;
    for( int p=1 ; p < m_npatches ; p++ )
    {
        double csmin=1e38,cpmin=1e38,cratmin=1e38,csmax=-1e38,cpmax=-1e38,cratmax=-1e38;
        double rhomin=1e38, rhomax=-1e38;
        for( int k=mMaterial_cs[p].m_kb ; k<= mMaterial_cs[p].m_ke ; k++ )
            for( int j=mMaterial_cs[p].m_jb ; j<= mMaterial_cs[p].m_je ; j++ )
                for( int i=mMaterial_cs[p].m_ib ; i<= mMaterial_cs[p].m_ie ; i++ )
                {
                    if( water || mMaterial_cs[p](1,i,j,k) != -999 )
                    {
                        if( mMaterial_rho[p](1,i,j,k) < rhomin )
                            rhomin = mMaterial_rho[p](1,i,j,k);
                        if( mMaterial_rho[p](1,i,j,k) > rhomax )
                            rhomax = mMaterial_rho[p](1,i,j,k);
                        if( mMaterial_cs[p](1,i,j,k) < csmin )
                            csmin = mMaterial_cs[p](1,i,j,k);
                        if( mMaterial_cs[p](1,i,j,k) > csmax )
                            csmax = mMaterial_cs[p](1,i,j,k);
                        if( mMaterial_cp[p](1,i,j,k) < cpmin )
                            cpmin = mMaterial_cp[p](1,i,j,k);
                        if( mMaterial_cp[p](1,i,j,k) > cpmax )
                            cpmax = mMaterial_cp[p](1,i,j,k);
                        double crat = mMaterial_cp[p](1,i,j,k)/mMaterial_cs[p](1,i,j,k);
                        if( crat < cratmin ) {
                            cratmin = crat;
                            if( printsmallcpcs && crat < 1.41 ) {
                                cout << "crat= " << crat << " at " << i << " " <<  j << " " << k << endl;
                                cout << " material is " << mMaterial_rho[p](1,i,j,k) << " " << mMaterial_cp[p](1,i,j,k) << " "
                                     << mMaterial_cs[p](1,i,j,k) << " " << mMaterial_qp[p](1,i,j,k) << " " << mMaterial_qs[p](1,i,j,k) << endl;
                            }
                        }
                        if( crat > cratmax )
                            cratmax = crat;
                    }
                }
        double cmins[4]= {csmin,cpmin,cratmin,rhomin}, cmaxs[4]= {csmax,cpmax,cratmax,rhomax};
        double cminstot[4], cmaxstot[4];
        MPI_Reduce(cmins, cminstot, 4, MPI_DOUBLE, MPI_MIN, 0, mEW->m_1d_communicator );
        MPI_Reduce(cmaxs, cmaxstot, 4, MPI_DOUBLE, MPI_MAX, 0, mEW->m_1d_communicator );
        int myid;
        MPI_Comm_rank(mEW->m_1d_communicator,&myid);
        if( myid == 0 )
            //	 if( mEW->getRank()==0 )
        {
            if( p== 1 && !water )
                cout << "S-file limits, away from water: " << endl;
            else if( p== 1 )
                cout << "S-file limits : " << endl;
            cout << "  Patch no " << p << " : " << endl;
            cout << "    cp    min and max " << cminstot[1] << " " << cmaxstot[1] << endl;
            cout << "    cs    min and max " << cminstot[0] << " " << cmaxstot[0] << endl;
            cout << "    cp/cs min and max " << cminstot[2] << " " << cmaxstot[2] << endl;
            cout << "    rho   min and max " << cminstot[3] << " " << cmaxstot[3] << endl;
        }
    }
}


