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
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
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
#ifndef _SACHDR_H
#define _SACHDR_H
/* this is the include file for the
        C language SAC codes
*/

#define True 1
#define False 0

#ifdef MSDOS
#define SW4_TYPE long
#else
#define SW4_TYPE sw4_type
#endif

/* data structures */

struct sachdr_ {
  float rhdr[70];
  SW4_TYPE ihdr[40];
  char chdr[24][8];
};

/* function prototypes */
void scmxmn(float *x, sw4_type npts, float *depmax, float *depmin, float *depmen);

void brsac(sw4_type npts, char *name, float **data, sw4_type *nerr);
void arsac(sw4_type npts, char *name, float **data, sw4_type *nerr);
void getfhv(char *strcmd, float *fval, sw4_type *nerr);
void getnhv(char *strcmd, sw4_type *ival, sw4_type *nerr);
void getkhv(char *strcmd, char *cval, sw4_type *nerr);
void getlhv(char *strcmd, sw4_type *lval, sw4_type *nerr);
void bwsac(sw4_type npts, const char *name, float *data);
void awsac(sw4_type npts, const char *name, float *data);
void setfhv(const char *strcmd, float fval, sw4_type *nerr);
void setnhv(const char *strcmd, sw4_type ival, sw4_type *nerr);
void setkhv(const char *strcmd, char *cval, sw4_type *nerr);
void setlhv(const char *strcmd, sw4_type lval, sw4_type *nerr);
void newhdr(void);
void inihdr(void);
void getihv(char *strcmd, char *strval, sw4_type *nerr);
void setihv(const char *strcmd, const char *strval, sw4_type *nerr);
sw4_type streql(const char *str1, const char *str2);

#endif
