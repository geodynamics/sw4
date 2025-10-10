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
#include "MaterialUCVM.h"
#include "Byteswapper.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif

using namespace std;


//-----------------------------------------------------------------------
MaterialUCVM::MaterialUCVM( EW* a_ew, const string a_file, const string a_directory):
    mEW(a_ew),
    m_model_file(a_file),
    m_model_dir(a_directory),
    m_use_attenuation(false)
{
    mCoversAllPoints = false;
    // Check that the depths make sense
    if (a_ew != NULL) {
        m_use_attenuation = a_ew->usingAttenuation();
        read_UCVM();
    }
}

//-----------------------------------------------------------------------
MaterialUCVM::~MaterialUCVM()
{
}

//-----------------------------------------------------------------------
void MaterialUCVM::set_material_properties(std::vector<Sarray> & rho,
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
    double lon, lat, elev, topo_elev;
    int nrow, nfile = 0;
    char inname[128], outname[128], cmd[2048];
    double squash_bottom = 7000;
    double squash_power = 1.1;

    FILE *fptr;

    double mylon, mylat, myz, myelev, myvs30, cvm_vp, cvm_vs, cvm_rho, gtl_vp, gtl_vs, gtl_rho, comb_vp, comb_vs, comb_rho;
    char cvm_name[128], gtl_name[128], comb_name[128];

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
    /* mEW->extractTopographyFromUCVM("ucvm"); */
    /* fprintf(stderr, "Done query topo\n"); */

    // The default floors for taper interpolation is 500m, 1700m, 1700m for vs,vp,density.
    const char *ucvm_cmd = getenv("UCVM_QUERY_CMD");
    if (ucvm_cmd == NULL || strstr(ucvm_cmd, "ucvm_query") == NULL) {
        if (mEW->getRank() == 0)
            fprintf(stderr, "UCVM_QUERY_CMD env variable not set correctly [%s], exiting...\n", ucvm_cmd);
        exit(-1);
    }

    // const char *ucvm_cmd = "ucvm_query -f /global/cfs/cdirs/m3354/tang/ucvm/install.25.7/conf/ucvm.conf -m cvmsi ";
    // const char *ucvm_cmd = "ucvm_query -f /global/cfs/cdirs/m3354/tang/ucvm/install.25.7/conf/ucvm.conf -m cvmsi,elygtl:taper -L 750,1700,1700 ";
    // const char *ucvm_cmd = "ucvm_query -f /global/cfs/cdirs/m3354/tang/ucvm/install.25.7/conf/ucvm.conf -m cvmsi,elygtl:taper -L 1000,1700,1700 ";
    if (mEW->getRank() == 0)
        fprintf(stderr, "Using command: %s\n", ucvm_cmd);

    for(int g=0; g < mEW->mNumberOfGrids; g++) {
        sprintf(inname, "/tmp/ucvm.in.%d.%d", g, mEW->getRank());
        sprintf(outname, "/tmp/ucvm.out.%d.%d", g, mEW->getRank());
        remove(outname);

        bool curvilinear = mEW->topographyExists() && g >= mEW->mNumberOfCartesianGrids;
        int npts = (mEW->m_jEnd[g] - mEW->m_jStart[g] + 1) * (mEW->m_iEnd[g] - mEW->m_iStart[g] + 1) * (mEW->m_kEnd[g] -  mEW->m_kStart[g] + 1);
        nrow = 0;
        nfile = 0;
        int total_batch = (int)npts / 20000;
        if (is_debug)
            fprintf(stderr, "Rank %d: Grid %d, %d points, %d batches, grid factor %d\n", mEW->getRank(), g, npts, (int)npts / 20000, g_fac[g]);

        // Query Vp, Vs, Rho
        fptr = fopen(inname, "w");
        for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i) {
            for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j) {
                // get lat lon for current grid g
                float_sw4 x = (i-1)*mEW->mGridSize[g];
                float_sw4 y = (j-1)*mEW->mGridSize[g];
                float_sw4 z = 0;
                // (x, y, z) is the coordinate of current grid point
                mEW->computeGeographicCoord(x, y, lon, lat);

                i1 = i * g_fac[g];
                j1 = j * g_fac[g];
                /* topo_elev = mEW->mTopo(i1, j1, 1); */

                for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k) {
                    if( curvilinear )
                        z = mEW->mZ[g](i,j,k);
                    else
                        z = mEW->m_zmin[g] + (k-1)*mEW->mGridSize[g];

                    elev = z;
                    if (elev < 0)
                        elev = 0;

                    /* // Bulldoze by removing anything above sea-level, filling constant values from vertical surface below sea-level */
                    /* if (bulldoze) { */
                    /*     // Adjust with topo */
                    /*     elev = z + topo_elev; */

                    /*     if (elev  < 0) */
                    /*         elev = 0; */

                    /*     /1* fprintf(stderr, "Rank %d grid %d: (%d, %d, %d) [%f, %f, %f] topo(%d, %d)=%f \n", *1/ */
                    /*     /1*         mEW->getRank(), g, i, j, k, lon, lat, z, i1, j1, topo_elev); *1/ */
                    /* } */
                    /* // Stretch & Squash by squashing above sea-level and stretching below sea-level with a bottom (7km) */
                    /* else { */
                    /*     if (z <= squash_bottom) { */
                    /*         if (z < 0) */
                    /*             z = 0; */
                    /*         // Even stretching */
                    /*         if (even_stretch) */
                    /*             elev = z / squash_bottom * (squash_bottom + topo_elev); */
                    /*         else */
                    /*             elev = pow(z / squash_bottom, squash_power) * (squash_bottom + topo_elev); */
                    /*     } */
                    /*     else */
                    /*         elev = z + topo_elev; */

                    /*     // Debug */
                    /*     if (nrow % 10000 == 0) { */
                    /*         fprintf(stderr, "Rank %d grid %d: (%d, %d, %d) [%f, %f, %f] topo(%d, %d)=%f, z=%f, squashed z=%f \n", */
                    /*                 mEW->getRank(), g, i, j, k, lon, lat, z, i1, j1, topo_elev, z, elev); */
                    /*     } */
                    /* } */

                    nrow++;
                    // Batch the processing to every 20000 points, the limit of ucvm_query
                    if (nrow >= 20000) {
                        fflush(fptr);
                        fclose(fptr);
                        if (is_debug)
                            fprintf(stderr, "Query batch %d / %d\n", nfile, total_batch);
                        // query UCVM and append to output file
                        sprintf(cmd, "%s < %s >> %s", ucvm_cmd, inname, outname);
                        system(cmd);

                        fptr = fopen(inname, "w");
                        nrow = 0;
                        nfile++;
                    }
                    fprintf(fptr, "%f %f %f\n", lon, lat, elev);
                    material++;
                } // End for k
            } // End for i
        } // End for j
        fclose(fptr);

        // query UCVM
        if (nrow > 0) {
            printf("Query last batch %d\n", nfile);
            sprintf(cmd, "%s < %s >> %s", ucvm_cmd, inname, outname);
            system(cmd);
        }

        if (is_debug)
            fprintf(stderr, "Rank %d: Grid %d, queried %d points\n", mEW->getRank(), g, material);

        material = 0;
        // read the output
        fptr = fopen(outname, "r");
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

                    elev = z;
                    if (elev < 0)
                        elev = 0;

                    /* if (bulldoze) { */
                    /*     // Bulldoze by removing anything above sea-level, fill constant values from vertical surface below sea-level */
                    /*     elev = z + topo_elev; */
                    /*     if (elev < 0) */
                    /*         elev = 0; */
                    /* } */
                    /* else { */
                    /*     // Squash by squashing above sea-level and stretching below sea-level with a bottom (7km) */
                    /*     if (z <= squash_bottom) { */
                    /*         if (z < 0) */
                    /*             z = 0; */
                    /*         // Even stretching */
                    /*         if (even_stretch) */
                    /*             elev = z / squash_bottom * (squash_bottom + topo_elev); */
                    /*         else */
                    /*             elev = pow(z / squash_bottom, squash_power) * (squash_bottom + topo_elev); */
                    /*     } */
                    /*     else */
                    /*         elev = z + topo_elev; */
                    /* } */


                    material++;
                    fscanf(fptr, " %lf %lf %lf %lf %lf %s %lf %lf %lf %s %lf %lf %lf %s %lf %lf %lf",
                                   &mylon, &mylat, &myz, &myelev, &myvs30, cvm_name, &cvm_vp, &cvm_vs, &cvm_rho,
                                   gtl_name, &gtl_vp, &gtl_vs, &gtl_rho, comb_name, &comb_vp, &comb_vs, &comb_rho);

                    rho[g](i, j, k) = comb_rho;
                    cp[g](i, j, k)  = comb_vp;
                    cs[g](i, j, k)  = comb_vs;
                    if( use_q ) {
                        // if (z <= 50)
                        //     xis[g](i, j, k) = 10.0;
                        // else if (z <= 100)
                        //     xis[g](i, j, k) = 20.0;
                        // else if (z < 200)
                        //     xis[g](i, j, k) = 30.0;
                        // else
                        //     xis[g](i, j, k)  = comb_vs / 1000.0 * 100.0;

                        // xip[g](i, j, k)  = xis[g](i, j, k) * 2.0;

			if (z <= 50)
			    xis[g](i, j, k) = 40.0;   //(damp  = 0.0125 ) 
			else if (z <= 100)
			    xis[g](i, j, k) = 90.0;  //(damp =  0.0056)
			else if (z < 200)
			    xis[g](i, j, k) = 140.0;    //(damp = 0.0036)
			else
			    xis[g](i, j, k) = fmax(comb_vs/1000.0*250, 140); //<== to avoid sudden Q reduction at 200 m) 

                        xip[g](i, j, k)  = xis[g](i, j, k) * 2.0;
                    }

                    if (fabs(lon - mylon) > 1e-3 )
                        printf("x=%.1ff, y=%.1f, sw4_lon=%f does not match ucvm_lon=%f!\n", x, y, lon, mylon);
                    if (fabs(lat - mylat) > 1e-3 )
                        printf("x=%.1f, y=%.1f, sw4_lat=%f does not match ucvm_lat=%f!\n", x, y, lat, mylat);
                    if (fabs(elev - myz) > 1e-3 )
                        printf("x=%.1f, y=%.1f, sw4_z=%f does not match ucvm_z=%f!\n", x, y, elev, myz);

                    if (comb_rho < 10 || comb_rho > 10000) {
                        fprintf(stderr, "Rank %d grid %d (%d, %d, %d) [%f, %f, %f] rho = %f\n",
                                         mEW->getRank(), g, i, j, k, mylon, mylat, myz, comb_rho);
                    }
                    if (comb_vp < 10 || comb_vp > 10000) {
                        fprintf(stderr, "Rank %d grid %d (%d, %d, %d) [%f, %f, %f] vp = %f\n",
                                         mEW->getRank(), g, i, j, k, mylon, mylat, myz, comb_vp);
                    }
                    if (comb_vs < 10 || comb_vs > 10000) {
                        fprintf(stderr, "Rank %d grid %d (%d, %d, %d) [%f, %f, %f] vs = %f\n",
                                         mEW->getRank(), g, i, j, k, mylon, mylat, myz, comb_vs);
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
        //           << "UCVM Initialized Node Types: " << endl
        //           << "   Material:        " << materialSum << endl
        //           << endl
        //           << "*Outside Domain:    " << outsideSum << endl
        //           << endl
        //           << "--------------------------------------------------------------\n"
        //           << endl;
        cout << endl
             << "UCVM command: outside = " << outsideSum << ", material = " << materialSum << endl;

    /* material_check(false); */
} // End of set_material_properties


//-----------------------------------------------------------------------
void MaterialUCVM::read_UCVM()
{
    // Timers
    double time_start, time_end;
    double intf_start, intf_end, mat_start, mat_end;
    time_start = MPI_Wtime();

    /* fill_in_fluids(); */

    time_end = MPI_Wtime();
    if (mEW->getRank() == 0) {
        cout << "MaterialUCVM::read_UCVM, time to read material file: " << time_end - time_start << " seconds." << endl;
    }
    cout.flush();
}

//-----------------------------------------------------------------------
void MaterialUCVM::fill_in_fluids()
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
void MaterialUCVM::material_check( bool water )
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


