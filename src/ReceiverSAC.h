#ifndef SW4_RECEIVER_SAC_H
#define SW4_RECEIVER_SAC_H
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace sw4 {
// Read the whole trace before exposing metadata or samples to the receiver.
struct SACTrace {
   std::array<float,70> real;
   std::array<int,40> integer;
   std::array<char,192> text;
   std::vector<float> samples;
   std::string field(int index) const {
      std::string value(text.data()+8*index,8);
      const auto end=value.find_last_not_of(" \0",std::string::npos,2);
      return end==std::string::npos ? "" : value.substr(0,end+1);
   }
   bool read(const std::string& path) {
      static_assert(sizeof(float)==4 && sizeof(int)==4,"SAC uses 32-bit fields");
      FILE* file=std::fopen(path.c_str(),"rb");
      if(!file) return false;
      bool ok=std::fread(real.data(),4,70,file)==70 &&
              std::fread(integer.data(),4,40,file)==40 &&
              std::fread(text.data(),1,192,file)==192;
      if(!ok) { std::fclose(file); return false; }
      bool swap=false;
      auto reverse=[](void* word) { auto p=static_cast<unsigned char*>(word); std::reverse(p,p+4); };
      if(ok && integer[6]!=6 && integer[6]!=7) {
         swap=true;
         for(auto& value:integer) reverse(&value);
         for(auto& value:real) reverse(&value);
      }
      const int n=integer[9];
      if(!ok || (integer[6]!=6 && integer[6]!=7) || n<1 ||
         !std::isfinite(real[0]) || real[0]<=0 || !std::isfinite(real[5]) ||
         integer[35]!=1 || integer[15]!=1) { std::fclose(file); return false; }
      // Bound allocation by the payload's actual size, including SAC v7 footer.
      const long start=std::ftell(file);
      ok=std::fseek(file,0,SEEK_END)==0;
      const long length=std::ftell(file);
      ok=ok && length>=start && (length-start)/4>=n && std::fseek(file,start,SEEK_SET)==0;
      if(ok) {
         samples.resize(n);
         ok=std::fread(samples.data(),4,n,file)==static_cast<size_t>(n);
         if(swap) for(auto& value:samples) reverse(&value);
         for(auto value:samples) if(!std::isfinite(value)) ok=false;
      }
      std::fclose(file);
      return ok;
   }
};
}
#endif
