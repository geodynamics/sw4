#ifndef __STATMACHINE__
#define __STATMACHINE__
#include <unordered_map>

template<typename T1, typename T2>
  class StatMachine{
 public:
  std::unordered_map<T1,std::vector<T2> > map;
  void insert(T1 arg1, T2 arg2){
    auto got = map.find(arg1);
    if (got==map.end()){
      map.emplace(arg1,std::vector<T2>());
    }
    map[arg1].push_back(arg2);
  }
  void print(ofstream &ofile){
    ofile<<"#Key Mean Median Min Max Count\n";
    
    for ( auto it : map){
      ofile<<it.first<<" ";
      std::sort(it.second.begin(),it.second.end());
      T2 sum(0);
      for (auto v : it.second) sum+= v;
      sum=sum/it.second.size();
      ofile<<sum<<" "<<it.second[it.second.size()/2]<<" "<<it.second[0]<<" "<<it.second.back()<<" "<<it.second.size()<<"\n";
    }
  }
  template<typename Func1, typename Func2>
    void print(ofstream &ofile, Func1 &&f1, Func2 &&f2){
    ofile<<"#Key Mean Median Min Max Psum Count\n";
    T2 grand_total(0);
    for ( auto it : map){
      std::sort(it.second.begin(),it.second.end());
      ofile<<f1(it.first)<<" ";
      T2 sum(0);
      for (auto v : it.second) { sum+= v; grand_total+=v;}
      T2 psum=sum;
      sum=sum/it.second.size();
      ofile<<f2(it.first,sum)<<" "
	   <<f2(it.first,it.second[it.second.size()/2])<<" "
	   <<f2(it.first,it.second[0])<<" "
	   <<f2(it.first,it.second.back())<<" "<<psum<<" "<<it.second.size()<<"\n";
    }
    ofile<<"# Grand total "<<grand_total<<"\n";
  }
  void printhistory(ofstream &ofile){
    int c=0;
    for ( auto it : map){
      ofile<<"#"<<it.first<<"\n";
      for (auto v : it.second) ofile<<c++<<" "<<v<<"\n";
      ofile<<"&\n";
    }
  }
  
};
#endif
