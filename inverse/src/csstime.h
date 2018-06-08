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
//  WPP LICENSE
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// # terms and conditions of the GNU General Public License for more details.
// # 
// # You should have received a copy of the GNU General Public License along with
// # this program; if not, write to the Free Software Foundation, Inc., 59 Temple
// # Place, Suite 330, Boston MA 02111-1307 USA.
// # ----------------------------------------------------------------------
//  WPP LICENSE
// # 
// # You should have received a copy of the GNU General Public License along with
// # this program; if not, write to the Free Software Foundation, Inc., 59 Temple
// # Place, Suite 330, Boston MA 02111-1307 USA.
// # ----------------------------------------------------------------------
//  WPP LICENSE
// # 
// # You should have received a copy of the GNU General Public License along with
// # this program; if not, write to the Free Software Foundation, Inc., 59 Temple
// # Place, Suite 330, Boston MA 02111-1307 USA.
// # ----------------------------------------------------------------------
//  WPP LICENSE
// You should have received a copy of the GNU General Public License along with
// this program; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston MA 02111-1307 USA.
// ----------------------------------------------------------------------
//  WPP LICENSE
#ifndef _CSSTIME_H
#define _CSSTIME_H

#define ISLEAP(yr)	((!(yr % 4) && yr % 100) || !(yr % 400))

struct date_time{
	double epoch;
	long date;
	int year;
	int month;
	char mname[4];
	int day;
	long doy;
	int hour;
	int minute;
	float second;
};

void htoe(struct date_time *dt) ;	/* convert from human to epoch	*/
double dtoepoch(long date) ;		/* convert julian date to epoch	*/
void month_day(struct date_time *dt) ;	/* from epoch fill in monty/day	*/
void etoh(struct date_time *dt) ;	/* epoch to human		*/
void mdtodate(struct date_time *dt);	/* from epoch to YYYY DOY */
void timestr(struct date_time *dt, char *str) ;
					/* 1999 12 31 23:59:59.999	*/
void timeprintstr(struct date_time *dt,char *str) ;
					/* epoch jday mon 12,1999 23:59:59.999*/

#endif
