#!/usr/bin/env python3
"""Check the production window setters and misfit weights without MPI or a solver.

The harness extracts the native methods and window-weight block; it tests
single-window selection, two-window union, overlap, and switching back to one.
"""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile


def body_end(source, start):
    opening = source.index("{", start)
    level = 0
    # These selected methods contain no braces in string literals or comments.
    for index in range(opening, len(source)):
        if source[index] == "{":
            level += 1
        elif source[index] == "}":
            level -= 1
            if level == 0:
                return index + 1
    raise ValueError("unterminated native function")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cxx", default="c++")
    parser.add_argument("--raja", action="store_true", help="Select unified GPU preprocessor paths")
    args = parser.parse_args()
    source = (Path(__file__).resolve().parents[1] / "src/TimeSeries.C").read_text()
    text = re.sub(r"^\s*#include[^\n]*", "", source, flags=re.M)
    command = [args.cxx, "-E", "-P", "-x", "c++"]
    if args.raja:
        command.append("-DSW4_USE_RAJA")
    source = subprocess.check_output(command + ["-"], input=text, text=True)
    setters = []
    for match in re.finditer(r"void TimeSeries::set_window\s*\([^)]*\)\s*\{", source):
        setters.append(source[match.start():body_end(source, match.start())])
    if len(setters) == 1:
        header = (Path(__file__).resolve().parents[1] / "src/TimeSeries.h").read_text()
        header = re.sub(r"^\s*#include[^\n]*", "", header, flags=re.M)
        header = subprocess.check_output(command + ["-"], input=header, text=True)
        for match in re.finditer(r"void set_window\s*\([^)]*\)\s*\{", header):
            method = header[match.start():body_end(header, match.start())]
            setters.append(method.replace("void set_window", "void TimeSeries::set_window", 1))
    assert len(setters) == 2, "two- and four-value window setters are required"
    defaults = []
    for field in ("m_winL", "m_winR", "m_winL2", "m_winR2"):
        match = re.search(r"\b" + field + r"\(([^)]+)\)", source)
        assert match, "window state must be initialized: " + field
        defaults.append(field + "=" + match.group(1))
    start = source.index("TimeSeries::misfit(")
    match = re.search(r"float_sw4\s+wghx\s*,\s*wghy\s*,\s*wghz\s*;", source[start:])
    assert match
    start += match.start()
    match = re.search(r"if\s*\(\s*m_use_win\s*\)", source[start:])
    assert match
    start += match.start()
    weight = source[start:body_end(source, start)]
    helper = ""
    if "window_weight(t)" in weight:
        match = re.search(r"float_sw4 TimeSeries::window_weight\([^)]*\) const\s*\{", source)
        assert match
        helper = source[match.start():body_end(source, match.start())]
    if "receiver_window_weight(" in weight:
        match = re.search(r"static float_sw4 receiver_window_weight\([^)]*\)\s*\{", source)
        assert match
        helper = source[match.start():body_end(source, match.start())]
    product = ""
    match = re.search(r"^float_sw4 TimeSeries::product_wgh\(", source, re.M)
    if match:
        product = source[match.start():body_end(source, match.start())]
    stub = r'''
#include <cmath>
#include <iostream>
#include <stdexcept>
#define float_sw4 double
using namespace std;
static double receiver_window_weight(double,double,double,double,double,double);
class TimeSeries {
 public:
  bool m_use_win=false;
  double m_dt=0.01, m_t0=3.9;
  bool m_use_x=true, m_use_y=true, m_use_z=true;
  int mLastTimeStep=5;
  double data[3][6], *mRecordedSol[3];
  TimeSeries() { for(int q=0;q<3;++q) { mRecordedSol[q]=data[q]; for(int i=0;i<6;++i) data[q][i]=1; } }
  double product_wgh(TimeSeries&) const;
  double DEFAULTS;
  void set_window(double, double);
  void set_window(double, double, double, double);
  double window_weight(double) const;
  double weight(double t) const {
    double itau=1/(5*m_dt), wghv=0.42, wghx=wghv, wghy=wghv, wghz=wghv;
    WEIGHT
    return wghx;
  }
};
'''.replace("DEFAULTS", ", ".join(defaults)).replace("WEIGHT", weight)
    main = r'''
static void close(double actual, double expected) {
  if (std::abs(actual-expected)>1e-12)
    throw std::runtime_error("window weight differs: actual="+std::to_string(actual)+
                             " expected="+std::to_string(expected));
}
static double interval(double t, double left, double right) {
  return 0.5*(std::tanh((t-left)/0.05)-std::tanh((t-right)/0.05));
}
int main() {
  TimeSeries single;
  close(single.weight(1.5),0.42);
  single.set_window(1,2);
  for(double t : {0.0,1.0,1.5,2.0,3.0,5.0}) close(single.weight(t),interval(t,1,2));
  single.set_window(1,2,4,5);
  for(double t : {0.0,1.0,1.5,2.0,3.0,4.0,4.5,5.0,6.0})
    close(single.weight(t),std::fmax(interval(t,1,2),interval(t,4,5)));
  single.set_window(1,3,2,4);
  for(double t : {0.0,1.0,2.0,2.5,3.0,4.0,5.0})
    close(single.weight(t),std::fmax(interval(t,1,3),interval(t,2,4)));
  single.set_window(1,2);
  close(single.weight(3.5),interval(3.5,1,2));
#ifdef HAVE_PRODUCT
  single.set_window(1,2,4,5);
  double expected=0;
  for(int i=0;i<=single.mLastTimeStep;++i) {
    double t=single.m_t0+i*single.m_dt;
    expected+=std::fmax(interval(t,1,2),interval(t,4,5));
  }
  close(single.product_wgh(single),3*expected);
  single.m_use_x=false;
  close(single.product_wgh(single),2*expected);
#endif
  std::cout<<"PASS: native single/two-window weights, overlap, and mode reset\n";
}
'''
    with tempfile.TemporaryDirectory(prefix="sw4-window-unit-") as name:
        root = Path(name)
        unit = root / "windows.cpp"
        unit.write_text(("#define HAVE_PRODUCT\n" if product else "") + stub + "\n".join(setters) + helper + product + main)
        subprocess.run([args.cxx, "-std=c++11", str(unit), "-o", str(root/"windows")], check=True)
        subprocess.run([str(root/"windows")], check=True)


if __name__ == "__main__":
    main()
