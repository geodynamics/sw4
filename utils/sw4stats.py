#!/usr/bin/env python3
# Reads stdout from SW4 run and outputs stats of interest
# along with a data for plotting variability
# Ramesh Pankajakshan May 3 2023
import time
import sys
import statistics
import math

def stats(times):

    print("#Mean =",statistics.mean(times))
    print("#Median =",statistics.median(times))
    print("#STDEV =",statistics.pstdev(times))
    return statistics.mean(times),statistics.median(times)

def mstats(times):
    # Drop outliers and look at stats
    print("\n\n #Stats after dropping outliers \n\n")
    start=int(len(times)/3)
    times.sort()
    middle=times[start:start+start]
    #print("Middle =",len(middle)," Full = ",len(times))
    stats(middle)

def main():
    sample="Jun 23 15:59:03 2022"
    s2="Thu Jun 23 16:07:42 2022"
    pattern = "%m %d %H:%M:%S %Y"
    
    i1=time.strptime(s2)
    epoch = int(time.mktime(time.strptime(s2)))
    count=0
    timestep=[]
    mem_usage=0
    step=0
    laststep=0
    gpssources=0.0
    source=0
    factor=1000000.0
    first_step=0
    mpi_init=0
    with open(sys.argv[1]) as infile:
        for line in infile:
            line=line.rstrip()
            if "After MPI_Init" in line:
                data=line.split()
                date=" ".join(data[2:])
                mpi_init=int(time.mktime(time.strptime(date)))
            if "number of grid point  sources" in line:
                gpsources = int(int(line.split()[6])/factor)
            elif "number of unique" in line:
                ugpsources = int(int(line.split()[6])/factor)
            elif "Number of time step" in line:
                total_steps = int(line.split()[5])
                dt= float(line.split()[7])
            elif "Number of point sources" in line:
                sources=int(int(line.split()[7])/factor)
                #print("Source = ",sources)
            elif "Min & Max" in line:
                mmax=line.split()[6]
                mem_usage+=float(mmax.split(",")[1])
            elif "Time step" in line:
                data=line.split()
                if len(data)<6 :
                    continue
                date=" ".join(data[6:])
                epoch=int(time.mktime(time.strptime(date)))
                t=0
                if count==0:
                    first_step=epoch
                if count>0:
                    t=epoch-last
                    print(data[2], t)
                    timestep.append(int(t))
                    if count==1:
                        step=int(data[2])-laststep
                    
                last=epoch
                laststep=int(data[2])
                count=count+1
    print("\n\n Timing stats \n\n")
    mean,median=stats(timestep)
    print("#min = ",min(timestep),"max = ",max(timestep))
    print("\n\n")
    print("#Number of point sources ",sources," M")
    print("# number of g p sources= ",gpsources," M")
    print("# number of unique g p sources= ",ugpsources," M")
    print("#Steps between timestamps",step)
    print("#Mem Usage ",mem_usage," GB")
    norm_time=mean/step/mem_usage*1000
    print("#Total steps to solution ",total_steps)
    print("#dt = ",dt)
    print("#Adjusted time per step",norm_time," s/per 1K steps/per GB")
    #print("#MPI_Init=",mpi_init," First step",first_step)
    startup_time = first_step-mpi_init
    print("#Startup time is ",startup_time," s")
    print("#Solve Time  ",total_steps/step*mean/3600," hours", total_steps/step*median/3600,"hours")
    TTS=(total_steps/step*mean+startup_time)/3600
    TTS_MEDIAN = (total_steps/step*median+startup_time)/3600
    print("#Total time to solution",TTS,"  hours",TTS_MEDIAN," hours")
    mstats(timestep)

if __name__=="__main__":
        main()
