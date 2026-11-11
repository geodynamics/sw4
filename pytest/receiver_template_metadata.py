#!/usr/bin/env python3
"""Check the SW4 receiver reader with real HDF5.

Run replay locally for metadata and read failures. Run groups inside an
authorized allocation for real MPI communicator isolation.
"""
import argparse
import os
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


STUBS = r'''
#include <mpi.h>
#include <hdf5.h>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
#define float_sw4 double
struct aborted {};
static int rank_id, opens, visits;
static const MPI_Comm station_group = (MPI_Comm)17;
static size_t cursor;
struct message { int count; MPI_Datatype type; vector<char> bytes; };
static vector<message> tape;
static int test_bcast(void *data, int count, MPI_Datatype type, int, MPI_Comm comm) {
  if (comm != station_group) throw runtime_error("metadata escaped the simulation communicator");
  size_t bytes = count * (type == MPI_DOUBLE ? sizeof(double) :
                         type == MPI_INT ? sizeof(int) : sizeof(char));
  if (rank_id == 0) {
    message value = {count, type, vector<char>(bytes)};
    if (bytes) memcpy(value.bytes.data(), data, bytes);
    tape.push_back(value);
  } else {
    if (cursor == tape.size()) throw runtime_error("missing MPI message");
    const message &value = tape[cursor++];
    if (value.count != count || value.type != type)
      throw runtime_error("inconsistent MPI stream");
    if (bytes) memcpy(data, value.bytes.data(), bytes);
  }
  return 0;
}
static int test_abort(MPI_Comm, int) { throw aborted(); }
static double test_time() {
  return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count();
}
static hid_t test_open(const char *name, unsigned flags, hid_t fapl) {
  ++opens;
  return H5Fopen(name, flags, fapl);
}
static herr_t test_iterate(hid_t file, H5_index_t index, H5_iter_order_t order,
                          hsize_t *position, H5L_iterate_t callback, void *data) {
  ++visits;
  return H5Literate(file, index, order, position, callback, data);
}
// The old all-rank reader is replayed independently using serial file access.
static herr_t test_mpio(hid_t, MPI_Comm, MPI_Info) { return 0; }
#define MPI_Bcast test_bcast
#define MPI_Abort test_abort
#define MPI_Wtime test_time
#define H5Fopen test_open
#undef H5Literate
#define H5Literate test_iterate
#define H5Pset_fapl_mpio test_mpio
#define ASSERT(x) do { if (!(x)) throw runtime_error(#x); } while (0)
#define CHECK_INPUT(x, msg) ASSERT(x)
class EW {
 public:
  MPI_Comm m_1d_communicator = station_group;
  bool topo;
  int mNumberOfCartesianGrids = 1;
  vector<double> m_zmin = {0};
  int getRank() { return rank_id; }
  int getVerbosity() { return 0; }
  bool topographyExists() { return topo; }
  void computeCartesianCoord(double &x, double &y, double lon, double lat) {
    x = (lon + 120) * 10; y = (lat - 30) * 10;
  }
};
struct fid { int id; };
class TimeSeries {
 public:
  enum receiverMode { Displacement, Velocity };
  string name, output, input;
  double x, y, z, wl = 0, wr = 0, shift = 0, scale = 1;
  double w1=0, w2=0, w3=0, w4=0;
  int mode, every, down, event;
  bool topo, xyz, observation = false, window = false;
  bool usex = true, usey = true, usez = true;
  fid *handle = NULL;
  TimeSeries *ts0 = NULL;
  static int next_id;
  TimeSeries(EW *, string, string n, receiverMode m, bool, bool, bool,
             string o, double a, double b, double c, bool t, int e, int d,
             bool component, int ev) : name(n), output(o), x(a), y(b), z(c),
             mode(m), every(e), down(d), event(ev), topo(t), xyz(component) {}
  bool getUseHDF5() { return true; }
  string getPath() { return "."; }
  string gethdf5FileName() { return output; }
  TimeSeries *getTS0Ptr() { return ts0; }
  void setTS0Ptr(TimeSeries *p) { ts0 = p; }
  void allocFid() { handle = new fid{++next_id}; }
  fid *getFidPtr() { return handle; }
  void setFidPtr(fid *p) { handle = p; }
  bool myPoint() { return (int(x) % 3) == rank_id; }
  bool readSACHDF5(EW *, string path, bool) { observation = true; input = path; return true; }
  void set_utc_to_simulation_utc() {}
  void set_window(double a, double b) { window = true; wl = a; wr = b; }
  void set_window(double a, double b, double c, double d) { w1=a; w2=b; w3=c; w4=d; }
  void exclude_component(bool a, bool b, bool c) { usex = a; usey = b; usez = c; }
  void add_shift(double value) { shift += value; }
  void set_scalefactor(double value) { scale = value; }
};
int TimeSeries::next_id = 0;
'''

REPLAY_MAIN = r'''
int main(int argc, char **argv) {
  if (argc != 5) return 2;
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
  bool failure = string(argv[3]) == "failure";
  int expected_reads = atoi(argv[4]);
  for (rank_id = 0; rank_id < 3; ++rank_id) {
    cursor = 0; opens = visits = 0; TimeSeries::next_id = 0;
    EW ew; ew.topo = true;
    vector<vector<TimeSeries *> > stations(1);
    bool failed = false;
    try {
      // Revisit the first output file to check sharing across receiver sets.
      for (int set = 0; set < (failure ? 1 : 3); ++set)
        readStationHDF5(&ew, argv[1], set == 1 ? "second.h5" : "first.h5",
                       137 + set, 3 + set, TimeSeries::Velocity, 0, &stations,
                       100, 100, true, set == 1, set == 1, 1, 12,
                       true, false, true, 1.25, true, 2.5);
    } catch (const aborted &) { failed = true; }
    if (failed != failure) throw runtime_error("inconsistent read failure");
    int wanted = rank_id == 0 || expected_reads != 0 ? (failure ? 1 : 3) : 0;
    if (opens != wanted || (!failure && visits != wanted))
      throw runtime_error("unexpected template open/traversal count");
    if (rank_id != 0 && cursor != tape.size())
      throw runtime_error("unconsumed MPI messages");
    ofstream output(string(argv[2]) + "." + to_string(rank_id));
    output << setprecision(17);
    for (TimeSeries *s : stations[0]) {
      output << s->name << " " << s->output << " " << s->x << " " << s->y
             << " " << s->z << " " << s->topo << " " << s->xyz << " "
             << s->mode << " " << s->every << " " << s->down << " "
             << s->event << " " << s->handle->id << " " << s->observation
             << " " << s->input << " " << s->window << " " << s->wl
             << " " << s->wr << " " << s->usex << " " << s->usey << " "
             << s->usez << " " << s->shift << " " << s->scale << " "
             << s->w1 << " " << s->w2 << " " << s->w3 << " " << s->w4 << "\n";
    }
    // Avoid deleting a shared handle until all receivers have been captured.
    for (TimeSeries *s : stations[0]) if (s->ts0 == s) delete s->handle;
    for (TimeSeries *s : stations[0]) delete s;
  }
}
'''


def station(stream, name, coordinates, *, geographic=False, orientation=None,
            left=None, right=None):
    group = stream.create_group(name)
    group.create_dataset("STLA,STLO,STDP" if geographic else "STX,STY,STZ",
                         data=np.asarray(coordinates, dtype=np.float64))
    if orientation is not None:
        group.create_dataset("ISNSEW", data=np.asarray([orientation], dtype=np.int32))
    if left is not None:
        group.create_dataset("WindowL", data=np.asarray([left], dtype=np.float64))
    if right is not None:
        group.create_dataset("WindowR", data=np.asarray([right], dtype=np.float64))
    return group


def extract(source):
    # Compile the actual changed station reader, replacing its solver dependencies.
    text = source[:source.index("static herr_t traverse_func2")]
    for header in ("EW.h", "Require.h", "TimeSeries.h", "readhdf5.h"):
        text = text.replace('#include "' + header + '"', "")
    return text + "\n#endif\n#endif\n"


GROUP_MAIN = r'''
int main(int argc,char **argv) {
  MPI_Init(&argc,&argv);
  try {
    int world_rank,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    if(world_size!=4 || argc!=3) throw runtime_error("requires four MPI ranks and two templates");
    int group=world_rank/2;
    MPI_Comm_split(MPI_COMM_WORLD,group,world_rank,&station_group);
    MPI_Comm_rank(station_group,&rank_id);
    EW ew; ew.topo=false;
    vector<vector<TimeSeries *> > stations(1);
    readStationHDF5(&ew,argv[group+1],"output.h5",137,3,
                   TimeSeries::Velocity,0,&stations,100,100,false,
                   false,false,0,0,true,true,true,0,false,1);
    int count=group==0 ? 8 : 2053;
    if((int)stations[0].size()!=count) throw runtime_error("incorrect group receiver count");
    for(int k=0;k<count;++k) {
      ostringstream name;
      name << (group==0 ? 'a' : 'b') << setfill('0') << setw(6) << k;
      TimeSeries *s=stations[0][k];
      if(s->name!=name.str() || s->x!=k%99 || s->y!=group+1 || s->z!=0)
        throw runtime_error("receiver metadata crossed simulation groups");
    }
    int total_opens,total_visits;
    MPI_Allreduce(&opens,&total_opens,1,MPI_INT,MPI_SUM,station_group);
    MPI_Allreduce(&visits,&total_visits,1,MPI_INT,MPI_SUM,station_group);
    if(total_opens!=1 || total_visits!=1) throw runtime_error("template not read exactly once per group");
    for(TimeSeries *s:stations[0]) if(s->ts0==s) delete s->handle;
    for(TimeSeries *s:stations[0]) delete s;
    MPI_Comm_free(&station_group);
    MPI_Barrier(MPI_COMM_WORLD);
    if(world_rank==0) cout << "PASS: two independent MPI groups retain their own receiver metadata; one read per group\n";
  } catch(const exception &error) {
    cerr << error.what() << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  MPI_Finalize();
}
'''


def run_replay(args):
    root = Path(__file__).resolve().parents[1]
    sources = {"candidate": (root / "src/readhdf5.C").read_text(),
               "baseline": subprocess.check_output(
                   ["git", "show", args.baseline_ref + ":src/readhdf5.C"],
                   cwd=root, text=True)}
    with tempfile.TemporaryDirectory(prefix="sw4-station-unit-") as temporary:
        work = Path(temporary)
        programs = {}
        for label, source in sources.items():
            unit = work / (label + ".cpp")
            stubs = STUBS
            if args.single_precision and label == "candidate":
                stubs = stubs.replace("#define float_sw4 double", "#define float_sw4 float")
            unit.write_text(stubs + extract(source) + REPLAY_MAIN)
            programs[label] = work / label
            subprocess.run([args.cxx, "-std=c++11", "-O0", "-g", "-DUSE_HDF5",
                            "-I" + str(args.hdf5_prefix / "include"),
                            "-I" + str(args.mpi_include), str(unit),
                            "-L" + str(args.hdf5_prefix / "lib"),
                            "-Wl,-rpath," + str(args.hdf5_prefix / "lib"),
                            "-lhdf5", "-o", str(programs[label])], check=True)
        fixtures = []
        mixed = work / "mixed.h5"
        with h5py.File(mixed, "w") as stream:
            stream.create_dataset("header", data=[0])
            station(stream, "a_cartesian", [1, 2, 3], orientation=0, left=2)
            station(stream, "b_default_orientation", [2, 3, 4])
            station(stream, "c_geographic", [31, -117, 2], geographic=True,
                    orientation=0, right=6)
            both = station(stream, "d_both_coordinates", [3, 4, 5], orientation=1)
            both.create_dataset("STLA,STLO,STDP", data=[32., -118., 1.])
            station(stream, "e_outside", [1000, 3, 4])
            station(stream, "f_above_cartesian_grid", [4, 5, -2])
            station(stream, "g_long_" + "s" * 250, [5, 6, 7], left=4, right=9)
            group = stream.create_group("h_nonstation")
            group.create_dataset("WindowL", data=[5.])
            station(stream, "i_after_nonstation", [6, 7, 8])
            depth = station(stream, "j_absolute_depth", [7, 8, 9])
            depth.create_dataset("USEZVALUE", data=np.array([1], dtype=np.int32))
            depth.create_dataset("WINDOWS", data=[1., 2., 3., 4.])
            station(stream, "k_after_windows", [8, 9, 10])
            scalar = station(stream, "l_scalar_fields", [9, 10, 11])
            scalar["ISNSEW"] = np.int32(0)
            scalar["USEZVALUE"] = np.int32(1)
            scalar["WindowL"] = np.float64(1)
            scalar["WindowR"] = np.float64(4)
        fixtures.append(mixed)
        for count in (0, 1024, 2053):
            path = work / ("count-" + str(count) + ".h5")
            with h5py.File(path, "w") as stream:
                for index in range(count):
                    station(stream, "s{:06d}".format(index),
                            [index % 99, (index * 3) % 99, 0],
                            orientation=index % 2,
                            left=2.0 if index == 1023 else None,
                            right=9.0 if index == 1024 else None)
            fixtures.append(path)
        # The baseline carries scalar bounds between stations. Check these
        # fields against the fixture contract instead of inheriting that result.
        # Each pair gives the result with no command bounds and with [1, 12].
        window_cases = {
            "a_cartesian": ((1, 2, 1e38), (1, 2, 12)),
            "c_geographic": ((1, -1, 6), (1, 1, 6)),
            "g_long_" + "s" * 250: ((1, 4, 9), (1, 4, 9)),
            "l_scalar_fields": ((1, 1, 4), (1, 1, 4)),
            "s001023": ((1, 2, 1e38), (1, 2, 12)),
            "s001024": ((1, -1, 9), (1, 1, 9)),
        }
        defaults = ((0, 0, 0), (1, 1, 12))
        for path in fixtures:
            for label, program in programs.items():
                subprocess.run([str(program), str(path), str(work / label),
                                "success", "1" if label == "baseline" else "0"],
                               check=True, stdout=subprocess.PIPE, text=True)
            for rank in range(3):
                suffix = "." + str(rank)
                left = (work / ("baseline" + suffix)).read_text().splitlines()
                right = (work / ("candidate" + suffix)).read_text().splitlines()
                assert len(left) == len(right), (path, rank)
                for a, b in zip(left, right):
                    a, b = a.split(" "), b.split(" ")
                    expected = (0, 0, 0)
                    if b[12] == "1":  # Observation metadata is set on the owner.
                        expected = window_cases.get(b[0], defaults)[b[1] == "second.h5"]
                    np.testing.assert_allclose(
                        [float(value) for value in b[14:17]], expected,
                        rtol=1e-6 if args.single_precision else 0, atol=0,
                        err_msg=str((path, rank, b[0], b[1])))
                    a, b = a[:14] + a[17:], b[:14] + b[17:]
                    assert len(a) == len(b), (path, rank)
                    for x, y in zip(a, b):
                        if x != y:
                            assert args.single_precision, (path, rank, a, b)
                            np.testing.assert_allclose(float(x), float(y),
                                                       rtol=1e-6, atol=1e-7,
                                                       err_msg=str((path, rank)))
            print("PASS:", path.name, "station-local scalar windows and baseline metadata",
                  "on three replayed ranks; only candidate root opens/traverses the template")
        malformed = work / "bad.h5"
        with h5py.File(malformed, "w") as stream:
            group = stream.create_group("station")
            group.create_dataset("STX,STY,STZ", data=[b"a", b"b", b"c"])
        partial = work / "bad-after-full-batch.h5"
        with h5py.File(partial, "w") as stream:
            for index in range(1024):
                station(stream, "a{:04d}".format(index), [1, 2, 3])
            group = stream.create_group("z_bad")
            group.create_dataset("STX,STY,STZ", data=[b"a", b"b", b"c"])
        invalid_sizes = []
        fields = (("ISNSEW", 1), ("USEZVALUE", 1), ("WindowL", 1),
                  ("WindowR", 1), ("STX,STY,STZ", 3),
                  ("STLA,STLO,STDP", 3), ("WINDOWS", 4))
        for field, expected in fields:
            for count in (expected - 1, expected + 1):
                path = work / (field.replace(",", "-") + "-size-" + str(count) + ".h5")
                with h5py.File(path, "w") as stream:
                    group = stream.create_group("station")
                    if not field.startswith("ST"):
                        group["STX,STY,STZ"] = [1., 2., 3.]
                    group[field] = np.ones(count, dtype=np.float64)
                invalid_sizes.append(path)
        partial_size = work / "bad-size-after-full-batch.h5"
        with h5py.File(partial_size, "w") as stream:
            for index in range(1024):
                station(stream, "a{:04d}".format(index), [1, 2, 3])
            group = station(stream, "z_bad", [1, 2, 3])
            group["WINDOWS"] = np.ones(5)
        invalid_sizes.append(partial_size)
        for path in [work / "missing.h5", malformed, partial] + invalid_sizes:
            subprocess.run([str(programs["candidate"]), str(path),
                            str(work / "failure"), "failure", "0"],
                           check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print("PASS:", path.name, "propagates read failure to every replayed rank")
    print("PASS: replayed MPI metadata and failure cases")


def run_groups(args):
    root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=False)
    stubs=STUBS.replace('static const MPI_Comm station_group = (MPI_Comm)17;',
                        'static MPI_Comm station_group;')
    for macro in ['#define MPI_Bcast test_bcast','#define MPI_Abort test_abort',
                  '#define MPI_Wtime test_time']:
        stubs=stubs.replace(macro,'')
    old='mode(m), every(e), down(d), event(ev), topo(t), xyz(component) {}'
    new='''mode(m), every(e), down(d), event(ev), topo(t), xyz(component) {
      int one=1,total;
      MPI_Allreduce(&one,&total,1,MPI_INT,MPI_SUM,station_group);
      if(total!=2) throw runtime_error("receiver construction used the wrong group");
    }'''
    assert old in stubs
    stubs=stubs.replace(old,new)
    unit=root/'groups.cpp'
    unit.write_text(stubs+extract((args.source_root/'src/readhdf5.C').read_text())+GROUP_MAIN)
    program=root/'groups'
    command=[args.cxx or str(args.mpi_prefix/'bin/mpicxx'),'-std=c++11','-DUSE_HDF5',
             '-I'+str(args.hdf5_prefix/'include'),'-I'+str(args.mpi_prefix/'include'),
             str(unit),'-L'+str(args.hdf5_prefix/'lib'),
             '-Wl,-rpath,'+str(args.hdf5_prefix/'lib'),'-lhdf5',
             '-Wl,-rpath,'+str(args.mpi_prefix/'lib'),'-o',str(program)]
    subprocess.run(command,check=True)
    files=[]
    for group,count in enumerate((8,2053)):
        filename=root/('group-'+str(group)+'.h5');files.append(filename)
        with h5py.File(filename,'w') as stream:
            for k in range(count):
                station=stream.create_group(('a' if group==0 else 'b')+'{:06d}'.format(k))
                station['STX,STY,STZ']=np.asarray([k%99,group+1,0],dtype=np.float64)
    launch=['srun','--overlap','--exact','--nodes=1','--ntasks=4',
            '--ntasks-per-node=4','--cpus-per-task=2']
    if args.job_id: launch.append('--jobid='+args.job_id)
    launch+=[str(program),*[str(path) for path in files]]
    result=subprocess.run(launch,cwd=root,env=dict(os.environ,OMP_NUM_THREADS='1'),
                          text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                          timeout=120)
    (root/'mpi.log').write_text(result.stdout)
    assert result.returncode==0,result.stdout[-8000:]
    print(result.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    hdf5 = argparse.ArgumentParser(add_help=False)
    hdf5.add_argument("--hdf5-prefix", type=Path, required=True,
                      help="HDF5 installation with include/ and lib/")
    modes = parser.add_subparsers(dest="mode", required=True)

    replay = modes.add_parser("replay", parents=[hdf5], help="Run local reader and failure checks with replayed MPI")
    replay.add_argument("--mpi-include", type=Path, required=True)
    replay.add_argument("--cxx", default="c++")
    replay.add_argument("--baseline-ref", required=True,
                        help="CPU baseline revision for metadata parity")
    replay.add_argument("--single-precision", action="store_true",
                        help="Compare float metadata state with the double-precision reference")
    replay.set_defaults(check=run_replay)

    groups = modes.add_parser("groups", parents=[hdf5], help="Test two MPI groups on four ranks in an allocation")
    groups.add_argument("--mpi-prefix", type=Path, required=True)
    groups.add_argument("--work-dir", type=Path, required=True,
                        help="New shared directory for the executable and results")
    groups.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    groups.add_argument("--job-id")
    groups.add_argument("--cxx", help="MPI C++ wrapper (default: MPI prefix/bin/mpicxx)")
    groups.set_defaults(check=run_groups)

    args = parser.parse_args()
    args.check(args)


if __name__ == "__main__":
    main()
