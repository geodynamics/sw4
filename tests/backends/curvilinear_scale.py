#!/usr/bin/env python3
"""On-demand sanitized regression for large curvilinear interface acceptance.

Compile the actual production scale expression and stopping/failure checks in a
small host probe. No huge grid, MPI allocation, pytest or CTest registration is
needed. Retain the extracted probe, compiler logs and JSON evidence.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess


def extract(text):
    scale = re.search(r'(?:const\s+)?float_sw4\s+scalef\s*=([^;]+);', text)
    if not scale:
        raise ValueError('Missing production scale expression')
    # Preserve the guard body, including the collective abort. In the probe an
    # MPI stub throws, so both invalid-scale and acceptance behavior are checked.
    guard = re.search(r'if\(\s*!std::isfinite\(scalef\).*?\n\s*\}', text, re.S)
    if not guard:
        raise ValueError('Missing pre-iteration invalid-scale guard')
    iteration = re.search(r'while\(\s*(std::isfinite\(maxres\).*?)\s*\)\s*\{', text, re.S)
    failure = re.search(r'// Never copy an unconverged interface[^\n]*\n\s*if\(\s*(.*?)\s*\)\s*\{', text, re.S)
    if not iteration or not failure:
        raise ValueError('Missing production convergence checks')
    return scale[1], guard[0], iteration[1], failure[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[2]/'src/CurvilinearInterface2.C')
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--cxx', default='g++')
    args = parser.parse_args()
    source = args.source.resolve(strict=True)
    scale, guard, iteration, failure = extract(source.read_text())
    root = args.work_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    report = {'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'scale_expression': scale, 'checks': {}}
    # Includes the largest safe square, the first overflowing square, a
    # rectangular sign wrap, a zero wrap, and the maximum positive int axes.
    dimensions = [(2,2), (1025,1025), (8193,8193), (46341,46341),
                  (46342,46342), (50001,50001), (65537,32769),
                  (65537,65537), (2147483647,2147483647)]
    for precision in ('double', 'float'):
        probe = root/f'probe-{precision}.C'
        probe.write_text('''#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
using float_sw4 = '''+precision+''';
struct EW { int m_global_nx[1], m_global_ny[1];
 int m_cartesian_communicator=0; int getRank() { return 0; } };
void MPI_Abort(int, int) { throw std::runtime_error("collective abort"); }
bool valid_scale(float_sw4 scalef) {
 EW ew; EW* m_ew=&ew; int m_gc=0;
 try {
'''+guard+'''
 } catch (const std::runtime_error&) { return false; }
 return true;
}
int main(int argc, char** argv) {
 if(argc!=3) return 2;
 EW ew; ew.m_global_nx[0]=std::atoi(argv[1]); ew.m_global_ny[0]=std::atoi(argv[2]);
 EW* m_ew=&ew; int m_gc=0;
 const float_sw4 scalef='''+scale+''';
 const float_sw4 maxres0=1, m_reltol=1e-6, m_abstol=1e-6;
 const int iter=0, m_maxit=500;
 if(!valid_scale(scalef)) return 3;
 for(float_sw4 invalid: {float_sw4(0), float_sw4(-1),
     std::numeric_limits<float_sw4>::infinity(),
     std::numeric_limits<float_sw4>::quiet_NaN()})
   if(valid_scale(invalid)) return 4;
 float_sw4 maxres=1;
 if(!('''+iteration+''') || !('''+failure+''')) return 5;
 maxres=1e-30;
 if(('''+iteration+''') || ('''+failure+''')) return 6;
 for(float_sw4 invalid: {std::numeric_limits<float_sw4>::infinity(),
                       std::numeric_limits<float_sw4>::quiet_NaN()}) {
   maxres=invalid;
   if(('''+iteration+''') || !('''+failure+''')) return 7;
 }
 std::cout << std::setprecision(17) << scalef << "\\n";
}
''')
        binary = root/f'probe-{precision}'
        command = [args.cxx, '-std=c++11', '-O2', '-fsanitize=undefined',
                   '-fno-sanitize-recover=all', str(probe), '-o', str(binary)]
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        (root/f'compile-{precision}.log').write_text(result.stdout+result.stderr)
        records = []
        for nx, ny in dimensions:
            run = subprocess.run([str(binary),str(nx),str(ny)], text=True, capture_output=True)
            expected = (nx-1)*(ny-1)
            if run.returncode:
                raise ValueError(f'{precision} {nx}x{ny}: {run.stderr}')
            actual = float(run.stdout.strip())
            if not math.isfinite(actual) or actual <= 0:
                raise ValueError('Nonpositive or nonfinite scale')
            if precision=='double' and actual != float(expected):
                raise ValueError(f'Incorrect double scale: {actual} != {expected}')
            if precision=='float' and not math.isclose(actual,expected,rel_tol=2e-7):
                raise ValueError(f'Incorrect float scale: {actual} != {expected}')
            records.append({'nx':nx,'ny':ny,'integer_area':expected,'scale':actual,
                            'unconverged_rejected':True,'small_residual_accepted':True})
        report['checks'][precision] = {'compile_command':command,'cases':records,
                                      'invalid_scales_rejected':4,'nonfinite_residuals_rejected':2}
    (root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: sanitized production scaling and convergence checks, double/float, 18 dimensions')


if __name__ == '__main__':
    main()
