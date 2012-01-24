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
#
# Run me as: python makempeg.py
#

def EncodeMPEG(moviename, xres, yres, wildcard, firstfile, lastfile, delta):
    import os
    paramFile = "%s.params" % moviename
    f = open(paramFile, "w")
    f.write("MPEG-2 Test Sequence, 30 frames/sec\n");
    f.write(wildcard + " /* name of source files */\n");
    f.write("-        /* name of reconstructed images ('-': don't store) */\n")
    f.write("-        /* name of intra quant matrix file     ('-': default matrix) */\n") 
    f.write("-         /* name of non intra quant matrix file ('-': default matrix) */\n")
    f.write("stat.out  /* name of statistics file ('-': stdout ) */\n")
    f.write("3         /* input picture file format: 0=*.Y,*.U,*.V, 1=*.yuv, 2=*.ppm 3=*.tif*/\n")
    f.write("356     /* number of frames %d*/\n"%((lastfile-firstfile)/delta))
    f.write("0         /* number of first frame */\n")#%(firstfile))
    f.write("00:00:00:00 /* timecode of first frame */\n")
    f.write("15        /* N (# of frames in GOP) */\n")
    f.write("3         /* M (I/P frame distance) */\n")
    f.write("0         /* ISO/IEC 11172-2 stream */\n")
    f.write("0         /* 0:frame pictures, 1:field pictures */\n")
    f.write("%d       /* horizontal_size */\n"%(xres))
    f.write("%d       /* vertical_size */\n"%(yres))
    f.write("1        /* aspect_ratio_information 1=square pel, 2=4:3, 3=16:9, 4=2.11:1 */\n")
    f.write("4         /* frame_rate_code 1=23.976, 2=24, 3=25, 4=29.97, 5=30 frames/sec. */\n")
    f.write("5000000.0 /* bit_rate (bits/s) */\n")
    f.write("112   /* vbv_buffer_size (in multiples of 16 kbit) */\n")
    f.write("0 /* low_delay  */\n")
    f.write("0 /* constrained_parameters_flag */\n")
    f.write("1 /* Profile ID: Simple = 5, Main = 4, SNR = 3, Spatial = 2, High = 1 */\n")
    f.write("4 /* Level ID:   Low = 10, Main = 8, High 1440 = 6, High = 4  */\n")
    f.write("0 /* progressive_sequence */\n")
    f.write("1 /* chroma_format: 1=4:2:0, 2=4:2:2, 3=4:4:4 */\n")
    f.write("0 /* video_format: 0=comp., 1=PAL, 2=NTSC, 3=SECAM, 4=MAC, 5=unspec. */\n")
    f.write("5 /* color_primaries */\n")
    f.write("5 /* transfer_characteristics */\n")
    f.write("5 /* matrix_coefficients */\n")
    f.write("%d       /* display_horizontal_size */\n"%(xres))
    f.write("%d       /* display_vertical_size */\n"%(yres))
    f.write("0     /* intra_dc_precision (0: 8 bit, 1: 9 bit, 2: 10 bit, 3: 11 bit */\n")
    f.write("1 /* top_field_first */\n")
    f.write("0 0 0 /* frame_pred_frame_dct (I P B) */\n")
    f.write("0 0 0 /* concealment_motion_vectors (I P B) */\n")
    f.write("1 1 1 /* q_scale_type  (I P B) */\n")
    f.write("1 0 0 /* intra_vlc_format (I P B)*/\n")
    f.write("0 0 0 /* alternate_scan (I P B) */\n")
    f.write("0 /* repeat_first_field */\n")
    f.write("0 /* progressive_frame */\n")
    f.write("0 /* P distance between complete intra slice refresh */\n")
    f.write("0 /* rate control: r (reaction parameter) */\n")
    f.write("0 /* rate control: avg_act (initial average activity) */\n")
    f.write("0 /* rate control: Xi (initial I frame global complexity measure) */\n")
    f.write("0 /* rate control: Xp (initial P frame global complexity measure) */\n")
    f.write("0 /* rate control: Xb (initial B frame global complexity measure) */\n")
    f.write("0 /* rate control: d0i (initial I frame virtual buffer fullness) */\n")
    f.write("0 /* rate control: d0p (initial P frame virtual buffer fullness) */\n")
    f.write("0 /* rate control: d0b (initial B frame virtual buffer fullness) */\n")
    f.write("2 2 11 11 /* P:  forw_hor_f_code forw_vert_f_code search_width/height */\n")
    f.write("1 1 3  3  /* B1: forw_hor_f_code forw_vert_f_code search_width/height */\n")
    f.write("1 1 7  7  /* B1: back_hor_f_code back_vert_f_code search_width/height */\n")
    f.write("1 1 7  7  /* B2: forw_hor_f_code forw_vert_f_code search_width/height */\n")
    f.write("1 1 3  3  /* B2: back_hor_f_code back_vert_f_code search_width/height */\n")

    f.close();
    cmd = "./mpeg2encode " + paramFile + ' ' + moviename
    print cmd
    os.system(cmd)
    #os.unlink(paramFile)

# Call it
EncodeMPEG("seychelles.mpeg", 1088, 528, "Seawave%05d", 2000, 10000, 100)
