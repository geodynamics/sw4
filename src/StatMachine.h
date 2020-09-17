#ifndef __STATMACHINE__
#define __STATMACHINE__
#include <unordered_map>

class StatMachineBase{
 public:
  static bool ProfilerOn;
};
template <typename T1, typename T2>
  class StatMachine:public StatMachineBase {
 public:
  std::unordered_map<T1, std::vector<T2> > map;
  void insert(T1 arg1, T2 arg2) {
    if (ProfilerOn){
      auto got = map.find(arg1);
      if (got == map.end()) {
	map.emplace(arg1, std::vector<T2>());
      }
      map[arg1].push_back(arg2);
    }
  }
  void print(ofstream &ofile) {
    ofile << "#Key Mean Median Min Max Count\n";

    for (auto it : map) {
      ofile << it.first << " ";
      std::vector<T2> copy(it.second);
      std::sort(copy.begin(), copy.end());
      T2 sum(0);
      for (auto v : copy) sum += v;
      sum = sum / copy.size();
      ofile << sum << " " << copy[copy.size() / 2] << " "
            << copy[0] << " " << copy.back() << " "
            << copy.size() << "\n";
    }
  }
  template <typename Func1, typename Func2>
    void print(ofstream &ofile, Func1 &&f1, Func2 &&f2,std::string &&str) {
    ofile << "#Key Mean Median Min Max Psum Count "<<str<<"\n";
    T2 grand_total(0);
    for (auto it : map) {
      std::vector<T2> copy(it.second);
      std::sort(copy.begin(), copy.end());
      ofile << f1(it.first) << " ";
      T2 sum(0);
      for (auto v : copy) {
        sum += v;
        grand_total += v;
      }
      T2 psum = sum;
      sum = sum / copy.size();
      ofile << f2(it.first, sum) << " "
            << f2(it.first, copy[it.second.size() / 2]) << " "
            << f2(it.first, copy[0]) << " "
            << f2(it.first, copy.back()) << " " << psum << " "
            << copy.size() << "\n";
    }
    ofile << "# Grand total " << grand_total << " "<<str<<" \n";
  }
  void printhistory(ofstream &ofile) {
    int c = 0;
    for (auto it : map) {
      ofile << "#" << it.first << "\n";
      for (auto v : it.second) ofile << c++ << " " << v << "\n";
      ofile << "&\n";
    }
  }
};
#endif
