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
#ifndef SACSUB_C
#define SACSUB_C

#include	<stdio.h>
#include	<stdlib.h>
#include 	<string.h>
#include	<string.h>
#include	<errno.h>
#include "sacsubc.h"
#include <ctype.h>
struct sachdr_ sachdr;
#ifdef MSDOS
#include <fcntl.h>
#endif
/* string patterns */

static int  NRSTR=70;
const char *rstr[] = {
"DELTA   ", "DEPMIN  ", "DEPMAX  ", "SCALE   ", "ODELTA  ", 
"B       ", "E       ", "O       ", "A       ", "FMT     ", 
"T0      ", "T1      ", "T2      ", "T3      ", "T4      ", 
"T5      ", "T6      ", "T7      ", "T8      ", "T9      ", 
"F       ", "RESP0   ", "RESP1   ", "RESP2   ", "RESP3   ", 
"RESP4   ", "RESP5   ", "RESP6   ", "RESP7   ", "RESP8   ", 
"RESP9   ", "STLA    ", "STLO    ", "STEL    ", "STDP    ", 
"EVLA    ", "EVLO    ", "EVEL    ", "EVDP    ", "FHDR40  ", 
"USER0   ", "USER1   ", "USER2   ", "USER3   ", "USER4   ",
"USER5   ", "USER6   ", "USER7   ", "USER8   ", "USER9   ", 
"DIST    ", "AZ      ", "BAZ     ", "GCARC   ", "SB      ", 
"SDELTA  ", "DEPMEN  ", "CMPAZ   ", "CMPINC  ", "XMINIMUM", 
"XMAXIMUM", "YMINIMUM", "YMAXIMUM", "ADJTM   ", "FHDR65  ", 
"FHDR66  ", "FHDR67  ", "FHDR68  ", "FHDR69  ", "FHDR70  " 
	};
static int NISTR=40;
const char *istr[] = {
"NZYEAR  ", "NZJDAY  ", "NZHOUR  ", "NZMIN   ", "NZSEC   ", 
"NZMSEC  ", "NVHDR   ", "NINF    ", "NHST    ", "NPTS    ", 
"NSNPTS  ", "NSN     ", "NXSIZE  ", "NYSIZE  ", "NHDR15  ", 
"IFTYPE  ", "IDEP    ", "IZTYPE  ", "IHDR4   ", "IINST   ", 
"ISTREG  ", "IEVREG  ", "IEVTYP  ", "IQUAL   ", "ISYNTH  ", 
"IDHR11  ", "IDHR12  ", "IDHR13  ", "IDHR14  ", "IDHR15  ", 
"IDHR16  ", "IDHR17  ", "IDHR18  ", "IDHR19  ", "IDHR20  ", 
"LEVEN   ", "LPSPOL  ", "LOVROK  ", "LCALDA  ", "LHDR5   " 
	};
static int NCSTR=24;
const char *cstr[] = {
"KSTNM   ", "KEVNM   ", "KEVNMC  ", "KHOLE   ", 
"KO      ", "KA      ", "KT0     ", "KT1     ", 
"KT2     ", "KT3     ", "KT4     ", "KT5     ", 
"KT6     ", "KT7     ", "KT8     ", "KT9     ", 
"KF      ", "KUSER0  ", "KUSER1  ", "KUSER2  ", 
"KCMPNM  ", "KNETWK  ", "KDATRD  ", "KINST   "
	};

static int NESTR=50;
const char *estr[] = {
	"ITIME   ", "IRLIM   ", "IAMPH   ", "IXY     ", "IUNKN   ", 
	"IDISP   ", "IVEL    ", "IACC    ", "IB      ", "IDAY    ", 
	"IO      ", "IA      ", "IT0     ", "IT1     ", "IT2     ", 
	"IT3     ", "IT4     ", "IT5     ", "IT6     ", "IT7     ", 
	"IT8     ", "IT9     ", "IRADNV  ", "ITANNV  ", "IRADEV  ", 
	"ITANEV  ", "INORTH  ", "IEAST   ", "IHORZA  ", "IDOWN   ", 
	"IUP     ", "ILLLBB  ", "IWWSN1  ", "IWWSN2  ", "IHGLP   ", 
	"ISRO    ", "INUCL   ", "IPREN   ", "IPOSTN  ", "IQUAKE  ", 
	"IPREQ   ", "IPOSTQ  ", "ICHEM   ", "IOTHER  ", "IGOOD   ", 
	"IGLCH   ", "IDROP   ", "ILOWSN  ", "IRLDTA  ", "IVOLTS  "
	} ;
static int NIVAL=50;
const char *Istr[] = {
	"IFTYPE  ", "IFTYPE  ", "IFTYPE  ", "IFTYPE  ", "IDEP    ",
	"IDEP    ", "IDEP    ", "IDEP    ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IRADNV  ", "ITANNV  ", "IRADEV  ",
	"ITANEV  ", "INORTH  ", "IEAST   ", "IHORZA  ", "IDOWN   ",
	"IUP     ", "ILLLBB  ", "IWWSN1  ", "IWWSN2  ", "IHGLP   ",
	"ISRO    ", "IEVTYP  ", "IEVTYP  ", "IEVTYP  ", "IEVTYP  ",
	"IEVTYP  ", "IEVTYP  ", "IEVTYP  ", "IQUAL   ", "IQUAL   ",
	"IQUAL   ", "IQUAL   ", "IQUAL   ", "ISYNTH  ", "IDEP    "    
};

float  fhdr_default = -12345.;
int    ihdr_default = -12345 ;
const char  *chdr_default = "-12345  ";

void brsac (int LN,char *name,float **data,int *nerr)
/*
	LN	I*4	length of data array
	name	C*	Name of file to be opened
	data	R*4	Array of trace values
	nerr	I*4	-1 file does not exist
			-2 data points in file exceed dimension
			-3 file is not a SAC binary file

	NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
	LN POINTS ARE USED

  This routine reads waveform data written in SAC binary format.
  Adapted from FORTRAN version written by by Hafidh A. A. Ghalib, 1988.
*/
{
	int maxpts, nread;
	int i;
	int has12345;
	FILE *fptr;
	/*
		determine if the file exists for reading
	*/
	if((fptr=fopen(name,"rb")) == NULL){
		sachdr.ihdr[ 9] = 0;
		*nerr = -1;
perror("fopen error in brsac:");
		return;
	}
#ifdef MSDOS
	setmode(fileno(fptr), O_BINARY);
#endif
	/*
		file exists get header information only
	*/
	
	*nerr = 0;
	if(fread(&sachdr,sizeof( struct sachdr_),1,fptr)!=1){
		*nerr=-3;
		fclose(fptr);
		return;
	}
	/* ensure that we have enough space */
	if(sachdr.ihdr[ 9] > LN){
			maxpts = LN ;
			sachdr.ihdr[ 9] = LN ;
			*nerr = -2 ;
	} else {
			maxpts = sachdr.ihdr[ 9] ;
			*nerr = 0 ;
	}
	/* there must be at least one -12345 for this to be a sac file!! */
	has12345 = 0;
	for (i=0;i < 70 ; i++){
		if(sachdr.rhdr[i] == -12345.0)has12345++;
	}
	for (i=0;i < 40 ; i++){
		if(sachdr.ihdr[i] == -12345)has12345++;
	}
	for (i=0;i < 24 ; i++){
		if(strncmp(sachdr.chdr[i],  "-12345  ",8)==0)has12345++;
	}
	if(has12345 == 0){
		*nerr = -3;
		fclose(fptr);
		return;
	}
	
	/* allocate the data values */
	if((*data = (float *)calloc(maxpts,sizeof(float))) == NULL){
		*nerr = -2;
		fclose(fptr);
		return;
	}
	nread = fread(*data,sizeof(float),maxpts,fptr);
	fclose (fptr) ;
	/* safety */
	if(sachdr.ihdr[9] > LN){
			maxpts = LN;
			sachdr.ihdr[ 9] = LN;
	} else {
		maxpts = sachdr.ihdr[ 9];
	}
}

void bwsac (int LN, const char *name,float *data)
/*
	LN	I*4	length of data array
	name	C*	Name of file to be opened
	data	R*4	Array of trace values


  This routine writes waveform data written in SAC binary format.
  Adapted from FORTRAN version written by by Hafidh A. A. Ghalib, 1988.
*/
{
	int nwrite;
	FILE *fptr;
	if((fptr=fopen(name,"wb")) == NULL){
perror("fopen error in bwsac:");
		return;
	}
#ifdef MSDOS
	setmode(fileno(fptr), O_BINARY);
#endif
	fwrite(&sachdr,sizeof( struct sachdr_),1,fptr);
	nwrite = fwrite(data,sizeof(float),LN,fptr);
	fclose (fptr) ;
}

void arsac (int LN, char *name,float **data,int *nerr)
/*
	IRU	I*4	FILE handle for IO
	LN	I*4	length of data array
	name	C*	Name of file to be opened
	data	R*4	Array of trace values
	nerr	I*4	-1 file does not exist
			-2 data points in file exceed dimension
			-3 file is not a SAC ascii file

	NOTE IF THE ASCII FILE HAS MORE THAN LN POINTS, THEN ONLY
	LN POINTS ARE USED

  This routine reads waveform data written in SAC ascii format.
  Adapted from FORTRAN version written by by Hafidh A. A. Ghalib, 1988.
*/
{
	int maxpts, nread;
	int i, j, j1, j2, k, l;
	int has12345;
	FILE *fptr;
	float rval;
	int ival;
	char cval[9];
	char cval0[9], cval1[9], cval2[9];
	char cstr[81];
	char fstr[16];
	char *endptr;
	int ls;
	/*
		determine if the file exists for reading
	*/
	if((fptr=fopen(name,"r")) == NULL){
		sachdr.ihdr[ 9] = 0;
		*nerr = -1;
perror("fopen error in arsac:");
		return;
	}
	/*
		file exists get header information only
	*/
	
	*nerr = 0;
	/* this is ugly - SAC ascii was originally a FORTRAN specification			scanf likes field separators so problems with
		reading a FORTRAN 4g15.7
	*/
	/* read in   float header values 5 per line */
	for(i=0;i< NRSTR; i+=5){
		fgets(cstr,80,fptr);
		for(k=0,j=0;k<5;k++,j+=15){
			strncpy(fstr,&cstr[j],15);
			fstr[15] = '\0';
			sachdr.rhdr[i+k] = (float)strtod(fstr,&endptr);
		}
	}
	/* read in integer header values 5 per line */
	for(i=0;i< NISTR; i+=5){
		fgets(cstr,80,fptr);
		for(k=0,j=0;k<5;k++,j+=10){
			strncpy(fstr,&cstr[j],10);
			fstr[10] = '\0';
			sachdr.ihdr[i+k] = (int)strtol(fstr,&endptr,10);
		}
	}
	/* read in strings - 3 eight character strings per line */
	for(i=0;i< NCSTR; i+=3){
		/* use fgets to get a string - scanf ignores a
			valid blank field */
		fgets(cstr, 80, fptr);
		strncpy(cval0,&cstr[ 0],8);
		strncpy(cval1,&cstr[ 8],8);
		strncpy(cval2,&cstr[16],8);
		cval0[8] = '\0';
		cval1[8] = '\0';
		cval2[8] = '\0';
		for(k=0 ; k < 3 ; k++){
			switch(k){
				case 0: 
					strncpy(&cval[0],cval0,9);	
					break;
				case 1: 
					strncpy(&cval[0],cval1,9);
					break;
				case 2: 
					strncpy(&cval[0],cval2,9);
					break;
			}
			ls = strlen(cval);
			/* check for printing characters */
			for(l=0; l< 8;l++){
				if(isprint(cval[l]) == 0 || l >= ls)
					cval[l] = ' ';
			}
			cval[8] = '\0';
			strncpy(sachdr.chdr[i+k], cval, 8);
		}
	}
	/* ensure that we have enough space */
	if(sachdr.ihdr[ 9] > LN){
			maxpts = LN ;
			sachdr.ihdr[ 9] = LN ;
			*nerr = -2 ;
	} else {
			maxpts = sachdr.ihdr[ 9] ;
			*nerr = 0 ;
	}
	/* there must be at least one -12345 for this to be a sac file!! */
	has12345 = 0;
	for (i=0;i < 70 ; i++){
		if(sachdr.rhdr[i] == -12345.0)has12345++;
	}
	for (i=0;i < 40 ; i++){
		if(sachdr.ihdr[i] == -12345)has12345++;
	}
	for (i=0;i < 24 ; i++){
		if(strncmp(sachdr.chdr[i],  "-12345  ",8)==0)has12345++;
	}
	if(has12345 == 0){
		*nerr = -3;
		fclose(fptr);
		return;
	}
	/* allocate the data values */
	if((*data = (float *)calloc(maxpts,sizeof(float))) == NULL){
		*nerr = -2;
		fclose(fptr);
		return;
	}
	for(j=0; j< maxpts; j++){
		fscanf(fptr,"%g",&rval);
		 (*data)[j] = rval;
	}
		
	fclose (fptr) ;
	/* safety */
	if(sachdr.ihdr[9] > LN){
			maxpts = LN;
			sachdr.ihdr[ 9] = LN;
	} else {
		maxpts = sachdr.ihdr[ 9];
	}
}

void awsac (int LN, const char *name,float *data)
/*
	LN	I*4	length of data array
	name	C*	Name of file to be opened
	data	R*4	Array of trace values


  This routine writes waveform data written in SAC binary format.
  Adapted from FORTRAN version written by by Hafidh A. A. Ghalib, 1988.
*/
{
	int nwrite;
	FILE *fptr;
	int i, j, j1, j2, k;
	char cval[9];
	float zero;
	if((fptr=fopen(name,"w")) == NULL){
perror("fopen error in awsac:");
		return;
	}
	/* write real header block */
	for(j=0 ; j < NRSTR; ){
		fprintf(fptr,"%#15.7g",sachdr.rhdr[j]);
		j++;
		if( (j%5) == 0)
		fprintf(fptr,"\n");
	}
	/* write integer header block */
	j1 = 0;
	j2 = 5;
	for(i=0 ; i < 8 ; i++){
		for(j=j1 ; j < j2 ; j++)
			fprintf(fptr,"%10d",sachdr.ihdr[j]);
		fprintf(fptr,"\n");
		j1 += 5;
		j2 += 5;
	}
	/* write character header block */
	j1 = 0;
	j2 = 3;
	for(i=0 ; i < 8 ; i++){
		for(j=j1 ; j < j2 ; j++){
			for(k=0;k<8;k++){
				if(isprint(sachdr.chdr[j][k]))
					cval[k] = sachdr.chdr[j][k];
				else
					cval[k] = ' ';
			}
			cval[8] = '\0';
			fprintf(fptr,"%8s",cval);
		}
		fprintf(fptr,"\n");
		j1 += 3;
		j2 += 3;
	}
	j1 = 0;
	j2 = 4;
	zero = 0.0;
	for(i=0 ; i < LN; i+= 5){
		for(j=i ; j < i+5; j++){
			if(j < LN )
				fprintf(fptr,"%#15.7g",data[j]);
			else
				fprintf(fptr,"%#15.7g",zero);
		}
		fprintf(fptr,"\n");
	}
	fclose (fptr) ;
}

void getfhv(char *strcmd,float *fval,int *nerr)
/*
	Get float header value

	strcmd	C*8	String to key on
	fval	R*4	Real value returned
	nerr	I*4	Error condition
				0 no error
				1336 Value not defined
				1337 Header variable does not exist
*/
{
	int i;
	*nerr = -1;
	*fval = -12345. ;
	for(i=0;i< NRSTR ;i++){
		if(streql(strcmd,rstr[i])) {
			*fval = sachdr.rhdr[i];
			*nerr = 0;
			break ;
		}
	}
}

void getnhv(char *strcmd,int *ival,int *nerr)
/*
	Get integer header value

	strcmd	C*8	String to key on
	ival	R*4	Real value returned
	nerr	I*4	Error condition
				0 no error
				1336 Value not defined
				1337 Header variable does not exist
*/
{
	int i;
	*ival = -12345 ;
	*nerr = -1;
	for(i=0;i< NISTR ;i++){
		if(streql(strcmd,istr[i])) {
			*ival = sachdr.ihdr[i];
			*nerr = 0;
			break ;
		}
	}
}

void getkhv(char *strcmd,char *cval,int *nerr)
/*
	Get character header value

	strcmd	C*8	String to key on
	cval	C	Real value returned
	nerr	I*4	Error condition
				0 no error
				1336 Value not defined
				1337 Header variable does not exist
*/
{
	int i;
	*nerr = -1;
	strncpy(cval,"-12345  ",8) ;
	for(i=0;i< NCSTR ;i++){
		if(streql(strcmd,cstr[i])) {
			strncpy(cval,sachdr.chdr[i],8);
			/* force proper C string */
			cval[8] = '\0';
			*nerr = 0;
			break ;
		}
	}
}

void getlhv(char *strcmd,int *lval,int *nerr)
/*
	Get logical header value

	strcmd	C*8	String to key on
	lval	I*4	Real value returned
	nerr	I*4	Error condition
				0 no error
				1336 Value not defined
				1337 Header variable does not exist
*/
{
	int ival;
	*nerr = -1;
	getnhv(strcmd,&ival,nerr);
	if(ival == 0)
		*lval = False;
	else
		*lval = True;
}
	
void getihv(char *strcmd,char *strval,int *nerr)
/*
	Get enumerated header value

	strcmd	C*8	String to key on
	strval	C*8	real value set
	nerr	I*4	Error condition
			0  no error
			1336 Header variable undefined
			1337 Header variable does not exist
*/
{
	int nval;
	getnhv(strcmd,&nval,nerr);
	strcpy(strval,"        ");
	if(*nerr == 0){
		if(nval >= 1 &&  nval <= NESTR){
			strncpy(&strval[0],estr[nval -1],8);
			strval[9] = '\0';
		}
	}
}

void setfhv(const char *strcmd, float fval, int *nerr)
/*
	Set float header value

        strcmd  C*8     String to key on
        fval    R*4     real value set
        nerr    I*4     Error condition
                                0  no error
                                1337 Header variable does not exist
*/
{
	int i;
	for(i = 0 ; i < NRSTR; i++)
		if(streql(strcmd,rstr[i])) {
			sachdr.rhdr[i] = fval;
			break;
	}
}

void setnhv(const char *strcmd, int ival, int *nerr)
/*
	Set integer header value

        strcmd  C*8     String to key on
        ival    I*4     integer value set
        nerr    I*4     Error condition
                                0  no error
                                1337 Header variable does not exist
*/
{
	int i;
	for(i = 0 ; i < NISTR; i++)
		if(streql(strcmd,istr[i])) {
			sachdr.ihdr[i] = ival;
			break;
	}
}

void setkhv(const char *strcmd, char *cval, int *nerr)
/*
	Set integer header value

        strcmd  C*8     String to key on
        cval    I*4     String value set
        nerr    I*4     Error condition
                                0  no error
                                1337 Header variable does not exist
*/
{
	int i, ls, j;
	for(i = 0 ; i < NCSTR; i++)
		if(streql(strcmd,cstr[i])) {
			ls = strlen(cval);
			for(j=0;j<8;j++){
				if(j < ls)
					sachdr.chdr[i][j] = cval[j];
				else
					sachdr.chdr[i][j] = ' ';
			}
			break;
		}
}

void setlhv(const char *strcmd, int lval, int *nerr)
/*
	Set integer header value

        strcmd  C*8     String to key on
        cval    I*4     String value set
        nerr    I*4     Error condition
                                0  no error
                                1337 Header variable does not exist
*/
{
	if(lval == 0)
		setnhv(strcmd,0,nerr);
	else
		setnhv(strcmd,1,nerr);
}

void setihv(const char *strcmd, const char *strval, int *nerr)
/*
        Set enumerated header value
 
        strcmd  C*8     String to key on
        strval  C*8     real value set
        nerr    I*4     Error condition
                        0  no error
                        1336 Header variable undefined
                        1337 Header variable does not exist
*/
{
	int i;
	if(streql(strcmd,"IDEP    ") && streql(strval,"IUNKN   "))
		setnhv("IDEP    ",5,nerr);
	else if(streql(strcmd,"IZTYPE  ") && streql(strval,"IUNKN   "))
		setnhv("IZTYPE  ",5,nerr);
	else if(streql(strcmd,"IEVTYP  ") && streql(strval,"IOTHER  "))
		setnhv("IEVTYP  ",44,nerr);
	else if(streql(strcmd,"IQUAL   ") && streql(strval,"IOTHER  "))
		setnhv("IWUAL   ",44,nerr);
	else {
		*nerr = 1336;
		/* IFTYPE */
		for(i=0 ; i < NESTR; i++){
			if(streql(strcmd,Istr[i])&&streql(strval,estr[i])){
				setnhv(strcmd,i+1,nerr);
				break;
			}
		}
	}
}

void newhdr()
{
	inihdr();
	sachdr.ihdr[15]  = 1;	/* ITIME */
	sachdr.ihdr[35] = 1;	/* LEVEN = TRUE */
	sachdr.ihdr[37] = 1;	/* LOVROK = TRUE */
	sachdr.ihdr[38] = 1;	/* LCALDA = TRUE */
	
}

void inihdr()
/* initialize sac header */
{
	int i,j;
	for(i=0 ; i < NRSTR; i++)
		sachdr.rhdr[i] = fhdr_default;
	for(i=0 ; i < NISTR; i++)
		sachdr.ihdr[i] = ihdr_default;
	sachdr.ihdr[6] = 6;
	sachdr.ihdr[7] = 0;
	sachdr.ihdr[8] = 0;
	for(i=35 ; i < 40; i++)
		sachdr.ihdr[i] = 0;
	for(i=0 ; i < NCSTR; i++)
		strncpy(sachdr.chdr[i],chdr_default,8);
}

void wsac1(char *ofile,float *y,int npts,float btime,float dt,int *nerr)
{
/*
c----- 
c	PURPOSE: WRITE AN EVENLY SPACED SAC FILE
c
c	write a binary sac file with evenly sampled data
c	ofile	Char	name of file
c	y	R	array of values
c	npts	I	number of points in data
c	btime	R	start time
c	dt	R	sample interval
c	maxpts	I	maximum number of points to read
c	nerr	I	error return
c-----
*/
	float e;
	float depmax, depmin, depmen;
	scmxmn(y,npts,&depmax,&depmin,&depmen);
	setfhv("DEPMAX", depmax, nerr);
	setfhv("DEPMIN", depmin, nerr);
	setfhv("DEPMEN", depmen, nerr);
	setnhv("NPTS    ",npts,nerr);
	setfhv("DELTA   ",dt  ,nerr);
	setfhv("B       ",btime  ,nerr);
	setihv("IFTYPE  ","ITIME   ",nerr);
	e = btime + (npts -1 )*dt;
	setfhv("E       ",e     ,nerr);
	setlhv("LEVEN   ",1,nerr);
	setlhv("LOVROK  ",1,nerr);
	setlhv("LCALDA  ",1,nerr);
	bwsac(1,ofile,y);
	*nerr = 0;
}



int  streql(const char *str1, const char *str2)
{
	int l1, l2,i;
	int retval;
	char str_1[9], str_2[9];
	l1 = strlen(str1);
	if(l1<0)l1=0;
	if(l1>8)l1=8;
	for(i=0;i<8;i++){
		if(i < l1)
			str_1[i] = str1[i];
		else
			str_1[i] = ' ';
	}
	str_1[8] = '\0';
	l2 = strlen(str2);
	if(l2<0)l2=0;
	if(l2>8)l2=8;
	for(i=0;i<8;i++){
		if(i < l2)
			str_2[i] = str2[i];
		else
			str_2[i] = ' ';
	}
	str_2[8] = '\0';
	if(strncmp(str_1, str_2, 8) == 0)
			return(True);
	return (False);
}


void scmxmn(float *x, int npts, float *depmax, float *depmin, float *depmen)
{
	int i;
	double sum;
	*depmax = -1.0e+38 ;
	*depmin =  1.0e+38 ;
	sum = 0.0e+00 ;
	for( i=0 ; i < npts ; i++){
		if( x[i] > *depmax) *depmax = x[i] ;
		if( x[i] < *depmin) *depmin = x[i] ;
		sum +=  x[i] ;
 	}
	if(npts > 0)
		*depmen = sum / npts ;
	else
		*depmen = -12345. ;
}
	


#endif /* sacsubc.h */
