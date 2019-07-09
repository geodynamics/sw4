#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <byteswap.h>
#include <hdf5.h>

#define sw4_float float
#define VALID_VALUE -998

sw4_float get_min_max(sw4_float *data, int n, sw4_float *min, sw4_float *max)
{
    int i;
    *max = *min = data[0];
    for (i = 1; i < n; i++) {
        if (data[i] > *max) 
            *max = data[i];
        else if (data[i] < *min) 
            *min = data[i];
    }
}

void smooth_interface(sw4_float *data, int maxIter, int imin, int imax, int jmin, int jmax)
{
    sw4_float rf=0.2; // rf<0.25 for stability
    int i, j, iter;
    int isize = imax - imin, jsize = jmax - jmin;
    int size = isize * jsize;
    sw4_float *tmp  = (sw4_float*)calloc(size, sizeof(sw4_float));

    // Laplacian filter
    for (iter=0; iter < maxIter; iter++) {
        for (i = imin+1; i < imax-1; ++i) {
            for (j = jmin+1; j < jmax-1; ++j) {
                tmp[i*jsize + j] = data[i*jsize + j] + rf*(data[(i+1)*jsize + j] + 
                                   data[(i-1)*jsize + j] + data[i*jsize + (j+1)] + 
                                   data[i*jsize + (j-1)] - 4.*data[i*jsize + j]);
            }
        }

        // Neumann boundary conditions
        for (j = jmin+1; j < jmax-1; ++j) {
            i = imin;
            tmp[i*jsize + j] = tmp[(i+1)*jsize + j];
            i = imax-1;
            tmp[i*jsize + j] = tmp[(i-1)*jsize + j];
        }

        for (i = imin+1; i < imax-1; ++i) {
            j = jmin;
            tmp[i*jsize + j] = tmp[i*jsize + (j+1)];
            j = jmax-1;
            tmp[i*jsize + j] = tmp[i*jsize + (j-1)];
        }

        // Corners
        i = imin;
        j = jmin;
        tmp[i*jsize + j] = tmp[(i+1)*jsize + (j+1)];

        i = imax - 1;
        j = jmin;
        tmp[i*jsize + j] = tmp[(i-1)*jsize + (j+1)];

        i = imin;
        j = jmax - 1;
        tmp[i*jsize + j] = tmp[(i+1)*jsize + (j-1)];

        i = imax - 1;
        j = jmax - 1;
        tmp[i*jsize + j] = tmp[(i-1)*jsize + (j-1)];

        // update solution
        memcpy(data, tmp, sizeof(sw4_float)*size);
    }// end for iter

    free(tmp);
}

sw4_float interpolate(sw4_float up_val,   sw4_float down_val, 
                      sw4_float up_depth, sw4_float down_depth, sw4_float cur_z)
{
    sw4_float res;

    if (down_val < VALID_VALUE) {
        /* fprintf(stderr, "Error with interpolation, down value is invalid, " */
        /*        "up_depth=%.2f, down_depth=%.2f, cur_depth=%.2f, up_val=%.2f, down_val=%.2f!\n", */
        /*         up_depth, down_depth, cur_z, up_val, down_val); */
        return VALID_VALUE - 1.0;
    }

    if (up_val < VALID_VALUE) {
        /* printf("Use down value only, up_depth=%.2f, down_depth=%.2f, cur_depth=%.2f, " */
        /*        "up_val=%.2f, down_val=%.2f!\n", up_depth, down_depth, cur_z, up_val, down_val); */
        return down_val;
    }

    if (up_depth == down_depth) {
        return down_val;
    }
    res = up_val + (down_val - up_val) * ((cur_z - up_depth) / (down_depth - up_depth));

    return res;
}


int main(int argc, char *argv[])
{
    FILE *fp = NULL;
    char *rfilename, *sfilename, *mercstr;
    int magic, prec, att, mlen, nb, onesw = 1, need_swap, is_debug;
    uint64_t i, j, elem_cnt, *data_cnt = NULL;

    double az, lon0, lat0;
    size_t rsize, *data_sizes = NULL, file_size;
    void **data = NULL;
    double *hh, *hv, *z0;
    int *nc, *ni, *nj, *nk;
    sw4_float *float_ptr, *new_float_ptr;
    int grid_nz[10];
    int ngrid;
    sw4_float tmp_data[10];
    sw4_float cur_z, mul_factor;
    int r_level, cnt = 0, sfile_idx, up_k, down_k, up_idx, down_idx, look_down_cnt, tmpk, debug;
    sw4_float top_v, bot_v, up_depth, down_depth;
    sw4_float up_rho, down_rho, up_p, down_p, up_s, down_s, up_qp, down_qp, up_qs, down_qs;
    hid_t h5_dtype;
    hid_t file_id, group_id, dataset_id, attribute_id, memtype, str_dataspace, string_type, attr_space, attr_space2, attr_id;
    hid_t dataset_space, attr_space3, attr_space_n;
    char group_name[128];
    hsize_t ndim, dims[3];
    herr_t status;
    char interface_name[32];
    int ii, jj, kk, ref_ratio = 2;
    sw4_float **top_data_f;
    sw4_float **bot_data_f;
    sw4_float rfile_depth[10];
    sw4_float sfile_depth[10] = {0};
    sw4_float bottom_depth;
    sw4_float sfile_zdiff[10] = {0};
    int sfile_nk[10];
    int n_processed;
    sw4_float **sfile_rho;
    sw4_float **sfile_p  ;
    sw4_float **sfile_s  ;
    sw4_float **sfile_qp ;
    sw4_float **sfile_qs ;
    hid_t m_group_id;


    is_debug = 1;

    rfilename = argv[1];

    fp = fopen(rfilename, "r");
    if (fp == NULL) {
        fprintf(stderr,"Error opening file [%s]\n", rfilename);
        return -1;
    }

    // Read header
    rsize = fread(&magic, sizeof(int), 1, fp);
    if (rsize != 1) {
        fprintf(stderr, "Error reading magic\n");
        return -1;
    }
    bswap_32(onesw);
    if (magic == 1) {
        need_swap = 0;
    }
    else if (magic == onesw) {
        need_swap = 1;
    }
    else {
        fprintf(stderr, "Cannot determine byte order from magic number on file\n");
        return -1;
    }

    rsize = fread(&prec, sizeof(int), 1, fp);
    if (rsize != 1) {
        fprintf(stderr, "Error reading prec\n");
        return -1;
    }
    if (need_swap == 1) 
        bswap_32(prec);
    
    if (prec != 4 && prec != 8) {
        fprintf(stderr, "Error with prec=%d, expecting 4 or 8\n", prec);
        return -1;
    }

    rsize = fread(&att, sizeof(int), 1, fp);
    if (rsize != 1) {
        fprintf(stderr, "Error reading att\n");
        return -1;
    }
    if (need_swap == 1) 
        bswap_32(att);

    rsize = fread(&az, sizeof(double), 1, fp);
    if (rsize != 1) {
        fprintf(stderr, "Error reading az\n");
        return -1;
    }
    if (need_swap == 1) 
        bswap_64(az);

    rsize = fread(&lon0, sizeof(double), 1, fp);
    if (rsize != 1) {
        fprintf(stderr, "Error reading lon0\n");
        return -1;
    }
    if (need_swap == 1) 
        bswap_64(lon0);

    rsize = fread(&lat0, sizeof(double), 1, fp);
    if (rsize != 1) {
        fprintf(stderr, "Error reading lat0\n");
        return -1;
    }
    if (need_swap == 1) 
        bswap_64(lat0);

    rsize = fread(&mlen, sizeof(int), 1, fp);
    if (rsize != 1) {
        fprintf(stderr, "Error reading mlen\n");
        return -1;
    }
    if (need_swap == 1) 
        bswap_32(mlen);

    mercstr = (char*)calloc(1, mlen+1);
    if (mercstr == NULL) {
        fprintf(stderr, "Error allocating space for mercstr\n");
        return -1;
    }
    rsize   = fread(mercstr, mlen, 1, fp);
    if (rsize != 1) {
        fprintf(stderr, "Error reading mercstr\n");
        return -1;
    }

    rsize = fread(&nb, sizeof(int), 1, fp);
    if (rsize != 1) {
        fprintf(stderr, "Error reading nb\n");
        return -1;
    }
    if (need_swap == 1) 
        bswap_32(nb);


    printf("magic=%d, prec=%d, att=%d, az=%f, lon0=%f, lat0=%f, mlen=%d, mercstr=[%s], nb=%d\n", 
            magic, prec, att, az, lon0, lat0, mlen, mercstr, nb);


    hh = (double*)malloc(nb * sizeof(double));
    hv = (double*)malloc(nb * sizeof(double));
    z0 = (double*)malloc(nb * sizeof(double));

    nc = (int*)malloc(nb * sizeof(int));
    ni = (int*)malloc(nb * sizeof(int));
    nj = (int*)malloc(nb * sizeof(int));
    nk = (int*)malloc(nb * sizeof(int));

    // Read block header
    for (i = 0; i < nb; i++) {
        rsize = fread(&hh[i], sizeof(double), 1, fp);
        if (rsize != 1) {
            fprintf(stderr, "Error reading hh\n");
            return -1;
        }
        if (need_swap == 1) 
            bswap_64(hh[i]);

        rsize = fread(&hv[i], sizeof(double), 1, fp);
        if (rsize != 1) {
            fprintf(stderr, "Error reading hv\n");
            return -1;
        }
        if (need_swap == 1) 
            bswap_64(hv[i]);

        rsize = fread(&z0[i], sizeof(double), 1, fp);
        if (rsize != 1) {
            fprintf(stderr, "Error reading z0\n");
            return -1;
        }
        if (need_swap == 1) 
            bswap_64(z0[i]);

        rsize = fread(&nc[i], sizeof(int), 1, fp);
        if (rsize != 1) {
            fprintf(stderr, "Error reading nc\n");
            return -1;
        }
        if (need_swap == 1) 
            bswap_32(nc[i]);

        if (i == 0 && nc[i] != 1) {
            fprintf(stderr, "Error with block %d nc=%d, expecting 1\n", i, nc[i]);
            return -1;
        }
        else if (i > 0 && att == 0 && nc[i] != 3) {
            fprintf(stderr, "Error with block %d nc=%d, expecting 3\n", i, nc[i]);
            return -1;
        }
        else if (i > 0 && att == 1 && nc[i] != 5) {
            fprintf(stderr, "Error with block %d nc=%d, expecting 5\n", i, nc[i]);
            return -1;
        }


        rsize = fread(&ni[i], sizeof(int), 1, fp);
        if (rsize != 1) {
            fprintf(stderr, "Error reading ni\n");
            return -1;
        }
        if (need_swap == 1) 
            bswap_32(ni[i]);

        rsize = fread(&nj[i], sizeof(int), 1, fp);
        if (rsize != 1) {
            fprintf(stderr, "Error reading nj\n");
            return -1;
        }
        if (need_swap == 1) 
            bswap_32(nj[i]);

        rsize = fread(&nk[i], sizeof(int), 1, fp);
        if (rsize != 1) {
            fprintf(stderr, "Error reading nk\n");
            return -1;
        }
        if (need_swap == 1) 
            bswap_32(nk[i]);

        if (i == 0 && nk[i] != 1) {
            fprintf(stderr, "Error with block %d nk=%d, expecting 1\n", i, nk[i]);
            return -1;
        }

        printf("block %d: hh=%6.2f, hv=%6.2f, z0=%8.2f, nc=%d, ni=%d, nj=%d, nk=%d\n" ,
                i, hh[i], hv[i], z0[i], nc[i], ni[i], nj[i], nk[i]);
    }

    for (j = 0; j < nb-1; j++) {
        rfile_depth[j] = z0[j+1];
    }
    rfile_depth[nb-1] = z0[nb-1] + (nk[nb-1]-1)*hv[nb-1];


    sw4_float topo_min, topo_max, topo_avg = 0;
    sw4_float all_p_min, all_p_max, p_min, p_max;
    sw4_float all_s_min, all_s_max, s_min, s_max;
    sw4_float all_qp_min, all_qp_max, qp_min, qp_max;
    sw4_float all_qs_min, all_qs_max, qs_min, qs_max;
    sw4_float all_rho_min, all_rho_max, rho_min, rho_max;
    int np_gt_2600 = 0, nqs_gt_1000 = 0, nqp_gt_2000, nqp_63000, nqs_40000;
    all_rho_min = all_p_min = all_s_min = all_qp_min = all_qs_min =  1000000;
    all_rho_max = all_p_max = all_s_max = all_qp_max = all_qs_max = -1000000;

    // Read block data
    // Read topography first
    data       = (void**)malloc(nb * sizeof(void*));
    data_cnt   = (uint64_t*)malloc(nb * sizeof(uint64_t));
    data_sizes = (size_t*)malloc(nb * sizeof(size_t));

    for (i = 0; i < nb; i++) {

        data_cnt[i]   = ni[i] * nj[i] * nk[i] * nc[i];
        data_sizes[i] = prec * data_cnt[i];
        data[i]       = calloc(1, data_sizes[i]);
        if (data[i] == NULL) {
            fprintf(stderr, "Error with block %d data allocation\n", i);
            return -1;
        }

        rsize = fread(data[i], data_sizes[i], 1, fp);
        if (rsize != 1) {
            fprintf(stderr, "Error reading topography\n");
            return -1;
        }

        float_ptr = (sw4_float*)data[i];

        if (need_swap == 1) {
            for (j = 0; j < data_cnt[i]; j++) {
                if (prec == 4) 
                    bswap_32(float_ptr[j]);
                else if (prec == 8) 
                    bswap_64(float_ptr[j]);
            }
        }
        printf("\n");

        printf("Block %d data (total_cnt=%d):\n", i, data_cnt[i]);
        for (j = 0; j < 10; j++) {
            if (prec == 4) 
                printf(" ,%.2f", float_ptr[j]);
            else if (prec == 8) 
                printf(" ,%.2f", float_ptr[j]);
        }
        printf(" ... ");
        for (j = data_cnt[i] - 10; j < data_cnt[i]; j++) {
            if (prec == 4) 
                printf(" ,%.2f", float_ptr[j]);
            else if (prec == 8) 
                printf(" ,%.2f", float_ptr[j]);
        }
        printf("\n");

        elem_cnt = ni[i] * nj[i] * nk[i] * nc[i]; 
        rho_min = p_min = s_min = qp_min = qs_min =  1000000;
        rho_max = p_max = s_max = qp_max = qs_max = -1000000;

        if (i == 0) {
            topo_min = topo_max = float_ptr[0];
            for (j = 1; j < elem_cnt; j++) {
                if (topo_min > float_ptr[j]) 
                    topo_min = float_ptr[j];
                if (topo_max < float_ptr[j]) 
                    topo_max = float_ptr[j];
                topo_avg += float_ptr[j];
            }
            topo_avg /= elem_cnt;
            if (is_debug) 
                printf("topo (min, max, avg) = (%.2f, %.2f, %.2f)\n", topo_min, topo_max, topo_avg);

        }
        else {
            j = 0;
            np_gt_2600 = 0;
            nqp_gt_2000 = 0;
            nqs_gt_1000 = 0;
            nqp_63000 = 0;
            nqs_40000 = 0;
            for (j = 0; j < elem_cnt; j += nc[i]) {
                if (float_ptr[j] > 2600) 
                    np_gt_2600++;
                
                if (float_ptr[j] > VALID_VALUE && rho_min > float_ptr[j]) 
                    rho_min = float_ptr[j];
                if (float_ptr[j] > VALID_VALUE && rho_max < float_ptr[j]) 
                    rho_max = float_ptr[j];
                if (float_ptr[j+1] > VALID_VALUE && p_min > float_ptr[j+1]) 
                    p_min = float_ptr[j+1];
                if (float_ptr[j+1] > VALID_VALUE && p_max < float_ptr[j+1]) 
                    p_max = float_ptr[j+1];
                if (float_ptr[j+2] > VALID_VALUE && s_min > float_ptr[j+2]) 
                    s_min = float_ptr[j+2];
                if (float_ptr[j+2] > VALID_VALUE && s_max < float_ptr[j+2]) 
                    s_max = float_ptr[j+2];
                if (nc[i] == 5) {
                    if (float_ptr[j+3] == 63000) 
                        nqp_63000++;
                    if (float_ptr[j+4] == 40000) 
                        nqs_40000++;

                    if (float_ptr[j+3] > 2000 && float_ptr[j+3] < 63000) {
                        nqp_gt_2000++;
                    }
                    if (float_ptr[j+4] > 1000 && float_ptr[j+4] < 40000) {
                        nqs_gt_1000++;
                    }
                    if (float_ptr[j+3] > VALID_VALUE && qp_min > float_ptr[j+3]) 
                        qp_min = float_ptr[j+3];
                    if (float_ptr[j+3] > VALID_VALUE && qp_max < float_ptr[j+3] && float_ptr[j+3] != 63000) 
                        qp_max = float_ptr[j+3];
                    if (float_ptr[j+4] > VALID_VALUE && qs_min > float_ptr[j+4]) 
                        qs_min = float_ptr[j+4];
                    if (float_ptr[j+4] > VALID_VALUE && qs_max < float_ptr[j+4] && float_ptr[j+4] != 40000) 
                        qs_max = float_ptr[j+4];
                }
            }
            if (rho_min < all_rho_min) 
                all_rho_min = rho_min;
            if (p_min < all_p_min) 
                all_p_min = p_min;
            if (s_min < all_s_min) 
                all_s_min = s_min;
            if (qp_min < all_qp_min) 
                all_qp_min = qp_min;
            if (qs_min < all_qs_min) 
                all_qs_min = qs_min;

            if (rho_max > all_rho_max) 
                all_rho_max = rho_max;
            if (p_max > all_p_max) 
                all_p_max = p_max;
            if (s_max > all_s_max) 
                all_s_max = s_max;
            if (qp_max > all_qp_max) 
                all_qp_max = qp_max;
            if (qs_max > all_qs_max) 
                all_qs_max = qs_max;

            // Debug print
            if (is_debug) {
                printf("rho (min, max) = (%.2f, %.2f)\n", rho_min, rho_max);
                printf("p   (min, max) = (%.2f, %.2f)\n", p_min, p_max);
                printf("s   (min, max) = (%.2f, %.2f)\n", s_min, s_max);
                if (nc[i] == 5) {
                    printf("qp  (min, max) = (%.2f, %.2f)\n", qp_min, qp_max);
                    printf("qs  (min, max) = (%.2f, %.2f)\n", qs_min, qs_max);
                }
                /* printf("%d points out of %d has p > 2600\n", np_gt_2600, elem_cnt / nc[i]); */
                printf("%d points out of %d has qp == 63000\n", nqp_63000, elem_cnt / nc[i]);
                printf("%d points out of %d has qs == 40000\n", nqs_40000, elem_cnt / nc[i]);
                printf("%d points out of %d has qp > 2000 && < 63000\n", nqp_gt_2000, elem_cnt / nc[i]);
                printf("%d points out of %d has qs > 1000 && < 40000\n", nqs_gt_1000, elem_cnt / nc[i]);
            }
        } // else
        printf("\n\n");
       
    } // For i < nb
    fflush(stdout);

    // Check if we have reach end of file
    rsize = ftell(fp);
    fseek(fp, 0, SEEK_END);
    file_size = ftell(fp);
    if (rsize != file_size) {
        fprintf(stderr, "Did not read the entire file, possible error?\n");
    }

    // Write to an HDF5 file
    attr_space = H5Screate(H5S_SCALAR);

    char out_fname[128];
    sprintf(out_fname, "%s.h5", rfilename);
    file_id = H5Fcreate(out_fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Z-interfaces
    group_id = H5Gcreate(file_id, "Z_interfaces", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0) {
        fprintf(stderr, "Error with H5Gcreate %d!", __LINE__);
        return -1;
    }

    if (prec == 4) 
        h5_dtype = H5T_NATIVE_FLOAT;
    else if (prec == 8) 
        h5_dtype = H5T_NATIVE_DOUBLE;

    top_data_f = (sw4_float**)calloc(nb, sizeof(sw4_float*));
    bot_data_f = (sw4_float**)calloc(nb, sizeof(sw4_float*));
    bottom_depth = z0[nb-1] + (nk[nb-1]-1)*hv[nb-1];

    /* sw4_float sfile_depth[10] = {0, 617.7, 2117.7, 6505.2, 45105.2}; */
    sfile_depth[0] = 0;
    sfile_depth[1] = 500 + topo_avg;
    sfile_depth[2] = 3000 + topo_avg;
    sfile_depth[3] = 6387.5 + topo_avg;
    sfile_depth[4] = 44987.5 + topo_avg;

    for (i = 1; i < nb; i++) 
        sfile_zdiff[i] = sfile_depth[i] - sfile_depth[i-1];
    
    printf("Bottom depth = %.1f\n", bottom_depth);

    /* make z values be the depth instead of elevation */
    float_ptr = (sw4_float*)data[0];
    for (i = 0; i < ni[0] * nj[0]; i++) {
        float_ptr[i] *= -1;
    }


    // Write topo fist
    sprintf(interface_name, "z_values_0");
    dims[0] = ni[0];
    dims[1] = nj[0];
    dataset_space = H5Screate_simple(2, dims, NULL);
    dataset_id    = H5Dcreate(group_id, interface_name, h5_dtype, dataset_space, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status        = H5Dwrite(dataset_id, h5_dtype, H5S_ALL, dataset_space, H5P_DEFAULT, data[0]);
    if (status < 0) {
        fprintf(stderr, "Error with H5Dwrite %d!", __LINE__);
        return -1;
    }
    H5Dclose(dataset_id);
    H5Sclose(dataset_space);


    for (i = 1; i < nb; i++) {
        elem_cnt = ni[i] * nj[i]; 
        top_data_f[i]  = (sw4_float*)calloc(elem_cnt, sizeof(sw4_float));
        bot_data_f[i]  = (sw4_float*)calloc(elem_cnt, sizeof(sw4_float));
        if (i == 1) {
            memcpy(top_data_f[i], data[0], elem_cnt * sizeof(sw4_float));
        }
        else {
            for (ii = 0; ii < ni[i]; ii++) 
                for (jj = 0; jj < nj[i]; jj++) 
                    top_data_f[i][ii*nj[i] + jj] = bot_data_f[i-1][ii*ref_ratio*nj[i-1] + jj*ref_ratio];

        }

        if (i == nb - 1) {
            for (ii = 0; ii < ni[i]; ii++) 
                for (jj = 0; jj < nj[i]; jj++) 
                        bot_data_f[i][ii*nj[i] + jj] = bottom_depth;
        }
        else {
            for (ii = 0; ii < ni[i]; ii++) 
                for (jj = 0; jj < nj[i]; jj++) 
                    bot_data_f[i][ii*nj[i] + jj] = top_data_f[i][ii*nj[i] + jj] + sfile_zdiff[i];

            smooth_interface(bot_data_f[i], 10, 0, ni[i], 0, nj[i]);
            /* // Make the bottom interface match the top interface of the patch below */
            /* for (ii = 0; ii < ni[i]; ii++) { */
            /*     for (jj = 0; jj < nj[i]; jj++) { */
            /*         if (ii % 2 == 0 && jj % 2 == 0) */ 
            /*             continue; */

            /*         if (ii % 2 != 0 && jj % 2 != 0) */ 
            /*             bot_data_f[i][ii*nj[i] + jj] = 0.25 * (bot_data_f[i][(ii-1)*nj[i] + jj-1] + bot_data_f[i][(ii+1)*nj[i] + jj+1] + */
            /*                                                    bot_data_f[i][(ii-1)*nj[i] + jj-1] + bot_data_f[i][(ii+1)*nj[i] + jj+1]); */
            /*         else if (ii % 2 != 0) */ 
            /*             bot_data_f[i][ii*nj[i] + jj] =  0.5 * (bot_data_f[i][(ii-1)*nj[i] + jj] + bot_data_f[i][(ii+1)*nj[i] + jj]); */
            /*         else */
            /*             bot_data_f[i][ii*nj[i] + jj] =  0.5 * (bot_data_f[i][ii*nj[i] + jj-1] + bot_data_f[i][ii*nj[i] + jj+1]); */
            /*     } */
            /* } */
        }

    }

    int s_level;
    for (s_level = 1; s_level < nb; s_level++) {
        sprintf(interface_name, "z_values_%d", s_level);

        dims[0] = ni[s_level];
        dims[1] = nj[s_level];
        dataset_space = H5Screate_simple(2, dims, NULL);

        dataset_id =  H5Dcreate(group_id, interface_name, h5_dtype, dataset_space, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, h5_dtype,H5S_ALL,dataset_space, H5P_DEFAULT, bot_data_f[s_level]);
        if (status < 0) {
            fprintf(stderr, "Error with H5Dwrite %d!", __LINE__);
            return -1;
        }

        H5Dclose(dataset_id);
        H5Sclose(dataset_space);
    }

    H5Gclose(group_id);

    sfile_rho = (sw4_float**)calloc(nb, sizeof(sw4_float*));
    sfile_p   = (sw4_float**)calloc(nb, sizeof(sw4_float*));
    sfile_s   = (sw4_float**)calloc(nb, sizeof(sw4_float*));
    sfile_qp  = (sw4_float**)calloc(nb, sizeof(sw4_float*));
    sfile_qs  = (sw4_float**)calloc(nb, sizeof(sw4_float*));

    m_group_id = H5Gcreate(file_id, "Material_model", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (m_group_id < 0) {
        fprintf(stderr, "Error with H5Gcreate %d!", __LINE__);
        return -1;
    }
 
    dims[0] = 2;
    attr_space2 = H5Screate_simple(1, dims, NULL);

    // Block 0
    dims[0] = 3;
    attr_space3 = H5Screate_simple(1, dims, NULL);
    for (s_level = nb - 1; s_level >= 1; s_level--) {
    /* for (s_level = 1; s_level < nb; s_level++) { */
        cnt = 0;

        sfile_nk[s_level] = sfile_zdiff[s_level] / hv[s_level] + 1;
        /* sfile_nk[s_level] = ceil(1.0*sfile_zdiff[s_level] / hv[s_level]) + 1; */
        printf("Sfile nk[%d]=%d\n", s_level, sfile_nk[s_level]);

        sfile_rho[s_level] = (sw4_float*)calloc(ni[s_level]*nj[s_level]*sfile_nk[s_level], sizeof(sw4_float));
        sfile_p[s_level]   = (sw4_float*)calloc(ni[s_level]*nj[s_level]*sfile_nk[s_level], sizeof(sw4_float));
        sfile_s[s_level]   = (sw4_float*)calloc(ni[s_level]*nj[s_level]*sfile_nk[s_level], sizeof(sw4_float));
        sfile_qp[s_level]  = (sw4_float*)calloc(ni[s_level]*nj[s_level]*sfile_nk[s_level], sizeof(sw4_float));
        sfile_qs[s_level]  = (sw4_float*)calloc(ni[s_level]*nj[s_level]*sfile_nk[s_level], sizeof(sw4_float));

        sprintf(group_name, "grid_%d", s_level-1); 
        group_id = H5Gcreate(m_group_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (group_id < 0) {
            fprintf(stderr, "Error with H5Gcreate %d!", __LINE__);
            return -1;
        }

        attr_id = H5Acreate(group_id, "Horizontal grid size", H5T_NATIVE_DOUBLE, 
                               attr_space, H5P_DEFAULT, H5P_DEFAULT);
        status  = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &hh[s_level]);
        if (status < 0) {
            fprintf(stderr, "Error with H5Awrite %d!", __LINE__);
            return -1;
        }
        H5Aclose(attr_id);

        // nc
        attribute_id = H5Acreate(group_id, "Number of components", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        if (attribute_id < 0) {
            fprintf(stderr, "Error with H5Acreate %d!", __LINE__);
            return -1;
        }

        status = H5Awrite(attribute_id, H5T_NATIVE_INT, &nc[s_level]);
        if (status < 0) {
            fprintf(stderr, "Error with H5Awrite %d!", __LINE__);
            return -1;
        }
        status = H5Aclose(attribute_id);

        // s_level starts from 1
        n_processed = 0;
        for (ii = 0; ii < ni[s_level]; ii++) {
            for (jj = 0; jj < nj[s_level]; jj++) {
                for (kk = 0; kk < sfile_nk[s_level]; kk++) {

                    /* // debug */
                    /* if (ii == 2499 && jj == 1257 ) { */
                    /*     debug = 1; */
                    /* } */

                    // Get the top and bottom depth of the point to calculate current point's depth
                    top_v = top_data_f[s_level][ii*nj[s_level] + jj];
                    bot_v = bot_data_f[s_level][ii*nj[s_level] + jj];
                    cur_z = (bot_v - top_v) * (1.0 * kk / (sfile_nk[s_level]-1)) + top_v;

                    // Find the rfile block number that the current point maps to
                    for (r_level = 1; r_level < nb; r_level++) {
                        if (cur_z <= rfile_depth[r_level] + 0.000001) {
                            break;
                        }
                    }

                    // rfile and sfile have the same ni[] and nj[] on same level blocks, but different nk
                    sfile_idx   = ii*nj[s_level]*sfile_nk[s_level] + jj*sfile_nk[s_level] + kk;
                    if (sfile_idx < 0 || sfile_idx >= ni[s_level]*nj[s_level]*sfile_nk[s_level]) {
                        printf("Error with sfile idx = %d\n!", sfile_idx);
                    }

                    float_ptr  = (sw4_float*)data[r_level];

                    int level_diff = r_level - s_level;
                    if (abs(level_diff) >= 2) {
                        fprintf(stderr, "Error! Block %d data mapped to block 2 level up/down in rfile\n", s_level);
                    }

                    if (r_level == s_level) {
                        // Current point in sfile is on the same block of rfile
                        // Find the two rfile points that are on top and bootom of this point
                        // rfile's top (base) depth is z0[s_level], bottom depth is z0[s_level]+(nk[s_level]-1)*hv[s_level]
                        mul_factor = 1;
                    }
                    else if (r_level > s_level) {
                        mul_factor = 0.5;
                        // Data is mapped to a coarser level in the rfile
                        /* printf("sfile block %d (%d, %d, %d) %.2f is in coarser block of rfile (%f to %f)\n", */ 
                        /*         s_level, ii, jj, kk, cur_z, rfile_depth[r_level-1], rfile_depth[r_level]); */
                        cnt++;
                    }
                    else {
                        mul_factor = 2;
                        // Data is mapped to a finer level in the rfile
                        if (s_level == 1) {
                            fprintf(stderr, "Error! Block 1 data cannot map to a finer block in rfile\n");
                        }
                        /* printf("sfile block %d (%d, %d, %d) %.2f is in finer block of rfile (%f to %f)\n", */ 
                        /*         s_level, ii, jj, kk, cur_z, rfile_depth[r_level-1], rfile_depth[r_level]); */
                        cnt++;
                    }

                    up_k       = floor((cur_z - z0[r_level]) / hv[r_level]);
                    down_k     =  ceil((cur_z - z0[r_level]) / hv[r_level]);
                    /* if (down_k == up_k) */ 
                    /*     down_k += 1; */
                    up_depth   = z0[r_level] + up_k * hv[r_level];
                    down_depth = z0[r_level] + down_k * hv[r_level];

                    // Sfile grid may be mapped to a different level in rfile
                    int tmp_ii = ii*mul_factor;
                    int tmp_jj = jj*mul_factor;
                    up_idx     = tmp_ii*nj[r_level]*nk[r_level]*nc[r_level] + tmp_jj*nk[r_level]*nc[r_level] + up_k*nc[r_level];
                    down_idx = tmp_ii*nj[r_level]*nk[r_level]*nc[r_level] + tmp_jj*nk[r_level]*nc[r_level] + down_k*nc[r_level];

                    if (up_idx < 0 || up_idx > ni[r_level]*nj[r_level]*nk[r_level]*nc[r_level]) 
                        printf("Error, up_idx = %d out of bound %d!\n", up_idx, ni[r_level]*nj[r_level]*nk[r_level]*nc[r_level]);
                    if (down_idx < 0 || down_idx > ni[r_level]*nj[r_level]*nk[r_level]*nc[r_level]) 
                        printf("Error, down_idx = %d out of bound %d!\n", up_idx, ni[r_level]*nj[r_level]*nk[r_level]*nc[r_level]);

                    up_rho     = float_ptr[up_idx];
                    down_rho   = float_ptr[down_idx];
                    up_p       = float_ptr[up_idx+1];
                    down_p     = float_ptr[down_idx+1];
                    up_s       = float_ptr[up_idx+2];
                    down_s     = float_ptr[down_idx+2];
                    if (nc[s_level] > 3) {
                        up_qp      = float_ptr[up_idx+3];
                        down_qp    = float_ptr[down_idx+3];
                        up_qs      = float_ptr[up_idx+4];
                        down_qs    = float_ptr[down_idx+4];
                    }

                    // Top level, when up value is significantly different from down value, use down value only
                    sw4_float max_diff = 50.0;
                    int use_down_only = 0;
                    if (s_level == 1 ) {
                    /* if (s_level == 1 && kk == 0) { */
                        if (up_rho > VALID_VALUE && (up_rho > max_diff * down_rho || up_rho < (1.0/max_diff) * down_rho)) {
                            printf("rho: sfile gr=%d, (%d, %d, %d), z=%.2f, intf=%.2f, up=%.2f, down=%.2f\n", 
                                    s_level, ii, jj, kk, cur_z, top_v, up_rho, down_rho);
                            use_down_only = 1;
                            up_rho = down_rho;
                        }
                        if (up_rho > VALID_VALUE && (up_p > max_diff * down_p || up_p < (1.0/max_diff) * down_p)) {
                            printf("cp : sfile gr=%d, (%d, %d, %d), z=%.2f, intf=%.2f, up=%.2f, down=%.2f\n", 
                                    s_level, ii, jj, kk, cur_z, top_v, up_p, down_p);
                            use_down_only = 1;
                            up_p   = down_p;
                        }
                        if (up_rho > VALID_VALUE && (up_s > max_diff * down_s || up_s < (1.0/max_diff) * down_s)) {
                            printf("cs : sfile gr=%d, (%d, %d, %d), z=%.2f, intf=%.2f, up=%.2f, down=%.2f\n", 
                                    s_level, ii, jj, kk, cur_z, top_v, up_s, down_s);
                            use_down_only = 1;
                            up_s   = down_s;
                        }
                        if (nc[s_level] > 3) {
                            if (up_rho > VALID_VALUE && (up_qp > max_diff * down_qp || up_qp < (1.0/max_diff) * down_qp)) {
                                printf("qp : sfile gr=%d, (%d, %d, %d), z=%.2f, intf=%.2f, up=%.2f, down=%.2f\n", 
                                        s_level, ii, jj, kk, cur_z, top_v, up_qp, down_qp);
                                use_down_only = 1;
                                up_qp = down_qp;
                            }
                            if (up_rho > VALID_VALUE && (up_qs > max_diff * down_qs || up_qs < (1.0/max_diff) * down_qs)) {
                                printf("qs : sfile gr=%d, (%d, %d, %d), z=%.2f, intf=%.2f, up=%.2f, down=%.2f\n", 
                                        s_level, ii, jj, kk, cur_z, top_v, up_qs, down_qs);
                                use_down_only = 1;
                                up_qs = down_qs;
                            }
                        }
                    }

                    if (use_down_only) {
                        up_rho = down_rho;
                        up_p   = down_p;
                        up_s   = down_s;
                        if (nc[s_level] > 3) {
                            up_qp = down_qp;
                            up_qs = down_qs;
                        }
                    }

                    // Debug
                    if (ii == 5 && jj == 60 && kk == 0) {
                        debug = 1;
                    }

                    look_down_cnt = 0;
                    while (down_s < VALID_VALUE) {
                        /* printf("(%d, %d, %d) z=%.2f, s=%.2f\n", ii, jj, kk, cur_z, down_s); */
                        look_down_cnt++;
                        tmpk = down_k + look_down_cnt;

                        new_float_ptr = float_ptr;
                        if (tmpk >= nk[s_level]) {
                            /* printf("Need to go to next sfile level\n"); */
                            tmpk = 0;
                            // Assume 2x coarser next level
                            down_idx = (ii/2)*nj[s_level+1]*sfile_nk[s_level+1] + (jj/2)*sfile_nk[s_level+1];
                            down_s   = sfile_s[s_level+1][down_idx];
                            up_s   = VALID_VALUE - 1;
                            if (down_s < 0) {
                                printf("Error with s\n");
                            }

                            down_rho = sfile_rho[s_level+1][down_idx];
                            up_rho = VALID_VALUE - 1;
                            down_p   = sfile_p[s_level+1][down_idx];
                            up_p   = VALID_VALUE - 1;

                            if (nc[s_level] > 3) {
                                up_qp   = VALID_VALUE - 1;
                                down_qp = sfile_qp[s_level+1][down_idx];
                                up_qs   = VALID_VALUE - 1;
                                down_qs = sfile_qs[s_level+1][down_idx];
                            }
                            break;
                        }

                        down_idx = tmp_ii*nj[r_level]*nk[r_level]*nc[r_level] + tmp_jj*nk[r_level]*nc[r_level] + tmpk*nc[r_level];
                        down_s   = new_float_ptr[down_idx+2];
                        up_s   = VALID_VALUE - 1;

                        down_rho = new_float_ptr[down_idx];
                        up_rho = VALID_VALUE - 1;
                        down_p   = new_float_ptr[down_idx+1];
                        up_p   = VALID_VALUE - 1;

                        if (nc[s_level] > 3) {
                            up_qp   = VALID_VALUE - 1;
                            down_qp = new_float_ptr[down_idx+3];
                            up_qs   = VALID_VALUE - 1;
                            down_qs = new_float_ptr[down_idx+4];
                        }
                    } // End while down_s < -998

                    /* if (look_down_cnt > 1) */ 
                    /*     printf("s (%d, %d, %d), z = %.2f, looked down %d, = %.2f\n", ii, jj, kk, cur_z, look_down_cnt, down_s); */
 
                    sfile_rho[s_level][sfile_idx] = interpolate(up_rho, down_rho, up_depth, down_depth, cur_z);
                    if ( sfile_rho[s_level][sfile_idx] < all_rho_min || sfile_rho[s_level][sfile_idx] > all_rho_max) 
                        printf("Patch %d, idx = %d (%d, %d, %d), z = %.2f, Rho = %.2f exceeds range (%.2f, %.2f)!\n", 
                                s_level, sfile_idx, ii, jj, kk, cur_z, sfile_rho[s_level][sfile_idx], all_rho_min, all_rho_max);

                    sfile_p[s_level][sfile_idx]   = interpolate(up_p, down_p, up_depth, down_depth, cur_z);
                    if ( sfile_p[s_level][sfile_idx] < all_p_min || sfile_p[s_level][sfile_idx] > all_p_max) 
                        printf("Patch %d, idx = %d (%d, %d, %d), z = %.2f, p = %.2f exceeds range (%.2f, %.2f)!\n", 
                                s_level, sfile_idx, ii, jj, kk, cur_z, sfile_p[s_level][sfile_idx], all_p_min, all_p_max);

                    sfile_s[s_level][sfile_idx]   = interpolate(up_s, down_s, up_depth, down_depth, cur_z);
                    if ( sfile_s[s_level][sfile_idx] < all_s_min || sfile_s[s_level][sfile_idx] > all_s_max) 
                        printf("Patch %d, idx = %d (%d, %d, %d), z = %.2f, s = %.2f exceeds range (%.2f, %.2f)!\n", 
                                s_level, sfile_idx, ii, jj, kk, cur_z, sfile_s[s_level][sfile_idx], all_s_min, all_s_max);

                    if (nc[s_level] > 3) {
                        sfile_qp[s_level][sfile_idx] = interpolate(up_qp, down_qp, up_depth, down_depth, cur_z);
                        if ( sfile_qp[s_level][sfile_idx] < all_qp_min || sfile_qp[s_level][sfile_idx] > all_qp_max)
                            printf("Patch %d, idx = %d (%d, %d, %d), z = %.2f, qp = %.2f exceeds range (%.2f, %.2f)!\n", 
                                    s_level, sfile_idx, ii, jj, kk, cur_z, sfile_qp[s_level][sfile_idx], all_qp_min, all_qp_max);

                        sfile_qs[s_level][sfile_idx] = interpolate(up_qs, down_qs, up_depth, down_depth, cur_z);
                        if ( sfile_qs[s_level][sfile_idx] < all_qs_min || sfile_qs[s_level][sfile_idx] > all_qs_max)
                            printf("Patch %d, idx = %d (%d, %d, %d), z = %.2f, qs = %.2f exceeds range (%.2f, %.2f)!\n", 
                                    s_level, sfile_idx, ii, jj, kk, cur_z, sfile_qs[s_level][sfile_idx], all_qs_min, all_qs_max);

                        if (sfile_qp[s_level][sfile_idx] > 60000) {
                            debug = 1;
                        }
                    }

                    /* if (cur_z > -24.9 && cur_z < 25.1 && sfile_qp[s_level][sfile_idx] > 62000) { */
                    /*     printf("(%d, %d, %d): Sfile z = %.2f qp = %.2f, up_qp = %.2f, down_qp = %.2f\n", */ 
                    /*             ii, jj, kk, cur_z, sfile_qp[s_level][sfile_idx], up_qp, down_qp); */
                    /* } */

                    n_processed++; 

                } //For kk
            } // For jj
        } // For ii

        printf("Processed %d/%d elements\n", n_processed, ni[s_level]*nj[s_level]*sfile_nk[s_level]);
        printf("Block mismatch count = %d/%d\n", cnt, ni[s_level]*nj[s_level]*sfile_nk[s_level]);

        // Write out Data
        ndim = 3;
        dims[0] = ni[s_level];
        dims[1] = nj[s_level];
        dims[2] = sfile_nk[s_level];

        dataset_space = H5Screate_simple(ndim, dims, NULL);

        dataset_id = H5Dcreate(group_id, "Rho", h5_dtype, dataset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, h5_dtype, H5S_ALL, dataset_space, H5P_DEFAULT, sfile_rho[s_level]);
        if (status < 0) {
            fprintf(stderr, "Error with H5Dwrite %d!", __LINE__);
            return -1;
        }
        H5Dclose(dataset_id);

        get_min_max(sfile_rho[s_level], ni[s_level]*nj[s_level]*sfile_nk[s_level], &rho_min, &rho_max);
        printf("rho min=%.2f, max=%.2f\n", rho_min, rho_max);

        dataset_id = H5Dcreate(group_id, "Cp", h5_dtype, dataset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, h5_dtype, H5S_ALL, dataset_space, H5P_DEFAULT, sfile_p[s_level]);
        if (status < 0) {
            fprintf(stderr, "Error with H5Dwrite %d!", __LINE__);
            return -1;
        }
        H5Dclose(dataset_id);
        get_min_max(sfile_p[s_level], ni[s_level]*nj[s_level]*sfile_nk[s_level], &p_min, &p_max);
        printf("cp min=%.2f, max=%.2f\n", p_min, p_max);

        dataset_id = H5Dcreate(group_id, "Cs", h5_dtype, dataset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, h5_dtype, H5S_ALL, dataset_space, H5P_DEFAULT, sfile_s[s_level]);
        if (status < 0) {
            fprintf(stderr, "Error with H5Dwrite %d!", __LINE__);
            return -1;
        }
        H5Dclose(dataset_id);
        get_min_max(sfile_s[s_level], ni[s_level]*nj[s_level]*sfile_nk[s_level], &s_min, &s_max);
        printf("cs min=%.2f, max=%.2f\n", s_min, s_max);

        if (nc[s_level] > 3) {
            dataset_id = H5Dcreate(group_id, "Qp", h5_dtype, dataset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataset_id, h5_dtype, H5S_ALL, dataset_space, H5P_DEFAULT, sfile_qp[s_level]);
            if (status < 0) {
                fprintf(stderr, "Error with H5Dwrite %d!", __LINE__);
                return -1;
            }
            H5Dclose(dataset_id);
            get_min_max(sfile_qp[s_level], ni[s_level]*nj[s_level]*sfile_nk[s_level], &qp_min, &qp_max);
            printf("qp min=%.2f, max=%.2f\n", qp_min, qp_max);

            dataset_id = H5Dcreate(group_id, "Qs", h5_dtype, dataset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataset_id, h5_dtype, H5S_ALL, dataset_space, H5P_DEFAULT, sfile_qs[s_level]);
            if (status < 0) {
                fprintf(stderr, "Error with H5Dwrite %d!", __LINE__);
                return -1;
            }
            H5Dclose(dataset_id);
            get_min_max(sfile_qs[s_level], ni[s_level]*nj[s_level]*sfile_nk[s_level], &qs_min, &qs_max);
            printf("qs min=%.2f, max=%.2f\n", qs_min, qs_max);
        }

        /* np_gt_2600 = 0; */
        /* for (jj = 0; jj < ni[s_level]*nj[s_level]*nk[s_level]; jj++) { */
        /*     if (sfile_p[s_level][jj] > 2600) { */
        /*         np_gt_2600++; */
        /*     } */
        /* } */
        /* printf("sfile: %d points out of %d has p > 2600\n", np_gt_2600, ni[s_level]*nj[s_level]*nk[s_level]); */

        H5Sclose(dataset_space);
        H5Gclose(group_id);


    }// End for s_level = 1 to nb
    H5Gclose(m_group_id);

    for (i = 1; i < nb; i++) {
        free(sfile_rho[i]);
        free(sfile_p[i]);
        free(sfile_s[i]);
        free(sfile_qp[i]);
        free(sfile_qs[i]);
    }
    /* // prec */
    /* attribute_id = H5Acreate(file_id, "Precision", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT); */
    /* if (attribute_id < 0) { */
    /*     fprintf(stderr, "Error with H5Acreate %d!", __LINE__); */
    /*     return -1; */
    /* } */

    /* status = H5Awrite(attribute_id, H5T_NATIVE_INT, &prec); */
    /* if (status < 0) { */
    /*     fprintf(stderr, "Error with H5Awrite %d!", __LINE__); */
    /*     return -1; */
    /* } */
    /* status = H5Aclose(attribute_id); */


    // att
    attribute_id = H5Acreate(file_id, "Attenuation", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attribute_id < 0) {
        fprintf(stderr, "Error with H5Acreate %d!", __LINE__);
        return -1;
    }

    status = H5Awrite(attribute_id, H5T_NATIVE_INT, &att);
    if (status < 0) {
        fprintf(stderr, "Error with H5Awrite %d!", __LINE__);
        return -1;
    }
    status = H5Aclose(attribute_id);

    // Origin lon lat az
    double longlataz[3];
    longlataz[0] = lon0;
    longlataz[1] = lat0;
    longlataz[2] = az;
    dims[0] = 3;
    attr_id = H5Acreate(file_id, "Origin longitude, latitude, azimuth", H5T_NATIVE_DOUBLE, attr_space3, 
                        H5P_DEFAULT, H5P_DEFAULT);
    status  = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, longlataz);
    if (status < 0) {
        fprintf(stderr, "Error with H5Awrite %d!", __LINE__);
        return -1;
    }
    H5Aclose(attr_id);

    // ngrid
    ngrid = (sw4_float)(nb - 1);
    attr_id = H5Acreate(file_id, "ngrids", H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    status  = H5Awrite(attr_id, H5T_NATIVE_INT, &ngrid);
    if (status < 0) {
        fprintf(stderr, "Error with H5Awrite %d!", __LINE__);
        return -1;
    }
    H5Aclose(attr_id);

    /* // grid nz */
    /* grid_nz[10]; */
    /* for (i = 0; i < nb-1; i++) { */
    /*     grid_nz[i] = sfile_nk[i+1]; */
    /* } */
    /* dims[0] = nb - 1; */
    /* attr_space_n = H5Screate_simple(1, dims, NULL); */
    /* attr_id = H5Acreate(file_id, "grid nz", H5T_NATIVE_INT, attr_space_n, H5P_DEFAULT, H5P_DEFAULT); */
    /* status  = H5Awrite(attr_id, H5T_NATIVE_INT, grid_nz); */
    /* if (status < 0) { */
    /*     fprintf(stderr, "Error with H5Awrite %d!", __LINE__); */
    /*     return -1; */
    /* } */
    /* H5Sclose(attr_space_n); */
    /* H5Aclose(attr_id); */

    attr_id = H5Acreate(file_id, "Coarsest horizontal grid spacing", H5T_NATIVE_DOUBLE, attr_space, 
                        H5P_DEFAULT, H5P_DEFAULT);
    status  = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &hh[0]);
    if (status < 0) {
        fprintf(stderr, "Error with H5Awrite %d!", __LINE__);
        return -1;
    }
    H5Aclose(attr_id);


    double minmax_depth[2];
    minmax_depth[0] = (double)-topo_max;
    minmax_depth[1] = (double)bottom_depth;
    attr_id = H5Acreate(file_id, "Min, max depth", H5T_NATIVE_DOUBLE, attr_space2, H5P_DEFAULT, H5P_DEFAULT);
    status  = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, minmax_depth);
    if (status < 0) {
        fprintf(stderr, "Error with H5Awrite %d!", __LINE__);
        return -1;
    }
    H5Aclose(attr_id);

    status = H5Sclose(attr_space);
    status = H5Sclose(attr_space2);
    status = H5Sclose(attr_space3);
    status = H5Fclose(file_id);

    if (fp)   fclose(fp);
    if (hh)   free(hh); 
    if (hv)   free(hv); 
    if (z0)   free(z0); 
    if (nc)   free(nc); 
    if (ni)   free(ni); 
    if (nj)   free(nj); 
    if (nk)   free(nk); 

    if (data_cnt)   free(data_cnt); 
    if (data_sizes) free(data_sizes); 
    if (data) {
        for (i = 0; i < nb; i++) 
            if (data[i]) free(data[i]);
        free(data);
    }

    for (i = 0; i < nb; i++) {
        free(top_data_f[i]);
        free(bot_data_f[i]);
    }
    free(top_data_f);
    free(bot_data_f);

    if (sfile_rho != NULL) free(sfile_rho);
    if (sfile_p   != NULL) free(sfile_p);
    if (sfile_s   != NULL) free(sfile_s);
    if (sfile_qp  != NULL) free(sfile_qp);
    if (sfile_qs  != NULL) free(sfile_qs);

    return 0;
}
