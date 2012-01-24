#  WPP LICENSE
# #  WPP LICENSE
# # ----------------------------------------------------------------------
# # WPP - Wave propagation Program
# # ----------------------------------------------------------------------
# # Copyright (C) 2008, Lawrence Livermore National Security, LLC.  
# # Produced at the Lawrence Livermore National Laboratory
# # 
# # Written by:
# # 
# # Bjorn Sjogreen   (sjogreen2@llnl.gov)
# # Anders Petersson  (andersp@llnl.gov)
# # 
# # Alums:
# # Stefan Nilsson      
# # Daniel Appelo
# # Kathleen McCandless (mccandless2@llnl.gov)
# # Caroline Bono
# # 
# # CODE-227123 All rights reserved.
# # 
# # This file is part of WPP, v2.0
# # 
# # Please also read docs/GPLLICENSE.txt which contains 
# # "Our Notice and GNU General Public License"
# # 
# # This program is free software; you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation; version 2, dated June 1991.
# # 
# # This program is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # terms and conditions of the GNU General Public License for more details.
# # 
# # You should have received a copy of the GNU General Public License along with
# # this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# # Place, Suite 330, Boston MA 02111-1307 USA.
# # ----------------------------------------------------------------------
#  WPP LICENSE
# # 
# # You should have received a copy of the GNU General Public License along with
# # this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# # Place, Suite 330, Boston MA 02111-1307 USA.
# # ----------------------------------------------------------------------
#  WPP LICENSE
# # 
# # You should have received a copy of the GNU General Public License along with
# # this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# # Place, Suite 330, Boston MA 02111-1307 USA.
# # ----------------------------------------------------------------------
#  WPP LICENSE
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston MA 02111-1307 USA.
# ----------------------------------------------------------------------
import sys, os, math

def usage():
    print
    print "Usage: python sumSources.py wpp.in"
    print
    sys.exit(0)



if __name__=='__main__':

    import sys
    if len(sys.argv) != 2:
        usage()

    wppfile = sys.argv[1]

    if not os.access(wppfile, os.R_OK):
        print "Cannot access file: " , wppfile

    f = open(wppfile, 'r')

    lines = f.readlines()

    numsources = 0
    momentsum = 0

    for l in lines:

        if l[0:6] == 'source':
            # Found source
            numsources += 1

            sourceTokens = l.split(' ')

            for token in sourceTokens:
                if token[0:2] == 'm0':
                    # Found moment term
                    [name, value] = token.split('=')
                    momentsum += eval(value)

            
    print "------------------------------------------------------"
    print " Processed ", numsources, " sources from: ", wppfile
    if numsources > 0:
        print " Total Seismic moment: ", momentsum, " Nm "
    if momentsum > 0:
        print " Moment magnitude: ", 2./3. * (math.log10(momentsum) - 9.1) 
    print "------------------------------------------------------"







