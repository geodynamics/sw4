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
#ifndef _SACHDR_H
#define _SACHDR_H
/* this is the include file for the
	C language SAC codes
*/

#define True 1
#define False 0

#ifdef MSDOS
#define INT long
#else
#define INT int
#endif


/* data structures */

struct sachdr_  {
	float rhdr[70];
	INT ihdr[40];
	char chdr[24][8];
	}  ;


/* function prototypes */
void scmxmn(float *x, int npts, float *depmax, float *depmin, float *depmen);

void brsac (int npts,char *name,float **data,int *nerr);
void arsac (int npts,char *name,float **data,int *nerr);
void getfhv(char *strcmd,float *fval,int *nerr);
void getnhv(char *strcmd,int *ival,int *nerr);
void getkhv(char *strcmd,char *cval,int *nerr);
void getlhv(char *strcmd,int *lval,int *nerr);
void bwsac (int npts, const char *name,float *data);
void awsac (int npts, const char *name,float *data);
void setfhv(const char *strcmd,float  fval,int *nerr);
void setnhv(const char *strcmd,int  ival,int *nerr);
void setkhv(const char *strcmd,char *cval,int *nerr);
void setlhv(const char *strcmd,int lval,int *nerr);
void newhdr(void);
void inihdr(void);
void getihv(char *strcmd,char *strval,int *nerr);
void setihv(const char *strcmd, const char *strval,int *nerr);
int streql(const char *str1, const char *str2);

#endif
