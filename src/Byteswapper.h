//  WPP LICENSE
// ----------------------------------------------------------------------
// WPP - Wave propagation Program
// ----------------------------------------------------------------------
// Copyright (C) 2007, Lawrence Livermore National Security, LLC.  
// Produced at the Lawrence Livermore National Laboratory
// 
// Written by:
// 
// Daniel   Appelo     (appelo2@llnl.gov) 
// Kathleen McCandless (mccandless2@llnl.gov)
// Anders   Petersson  (andersp@llnl.gov)
// Bjorn    Sjogreen   (sjogreen2@llnl.gov)
// 
// Alums:
// Stefan Nilsson      
// 
// CODE-227123 All rights reserved.
// 
// Thie file is part of WPP, v1.1
// 
// Please also read docs/GPLLICENSE.txt which contains 
// "Our Notice and GNU General Public License"
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2, dated June 1991.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// terms and conditions of the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along with
// this program; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston MA 02111-1307 USA.
// ----------------------------------------------------------------------
//  WPP LICENSE
#ifndef WPP_BYTESWAPPER_H
#define WPP_BYTESWAPPER_H

class Byteswapper
{
   inline void swap_bytes (volatile void *ptr, unsigned int i, unsigned int j)
   {
      volatile char *t = (volatile char *)(ptr);
  
     char tmp = t[i];
     t[i] = t[j];
     t[j] = tmp;
   }
 
   inline void swap_bytes8 (volatile void *ptr)
   {
      swap_bytes (ptr, 0, 7);
      swap_bytes (ptr, 1, 6);
      swap_bytes (ptr, 2, 5);
      swap_bytes (ptr, 3, 4);
   }
 
 
   inline void swap_bytes4 (volatile void *ptr )
   {
      swap_bytes (ptr, 0, 3);
      swap_bytes (ptr, 1, 2);
   }
 
 
public:

   bool bigend()
   {
     union
     {
        long l;
        char c[sizeof (long)];
     } u;
     u.l = 1;
     return (u.c[sizeof (long) - 1] == 1);
   }

   void byte_rev( void* data, int sz, const char* typ )
   {
      volatile char *t = (volatile char *)(data);
      if( strcmp(typ,"double") == 0 )
      {
         for( int i=0 ; i<sz ; i++ )
         {
            swap_bytes8( t );
            t += 8;
         }
      }
      else if( strcmp(typ,"float") == 0 || strcmp(typ,"int") == 0 )
      {
         for( int i=0 ; i<sz ; i++ )
         {
            swap_bytes4( t );
            t += 4;
         }
      }
   }
};
 
#endif
