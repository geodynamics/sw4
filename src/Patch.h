#ifndef SW4_PATCH_H
#define SW4_PATCH_H

class Patch
{
public:   
   int m_procid; // ID of processor to send patch to, or to receive patch from. 
                 // Sending or receiving depending on context.
   int m_ib, m_ie, m_jb, m_je, m_kb, m_ke; // Size of array patch in my processor

   Patch( int dims[6], int procid );

   template<class T>
      void pack( T* array, AllDims& dims, T* array_patch );

   template<class T>
      void unpack( T* array, AllDims& dims, T* array_patch );

   template<class T>
      void selfcopy( AllDims& src, T* src_array, AllDims& dest, T* dest_array );
      
   size_t size();
};

#endif
