# SW4 VisIt plugin - BETA instructions

This README file contains vague instructions on building a [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) plugin
that can read [SW4](https://computation.llnl.gov/projects/serpentine-wave-propagation) files, 2D `\*.sw4img` and 3D or "volimage" `\*.3Dimg` files.

**NOTE**: This plugin is not supported by the VisIt or SW4 team in any substantial
way, but it's being provided as a potentially useful tool for visualization
of SW4 data.

## Building against installed version of VisIt on MacOS
The following steps worked on MacOS 10.13.6. There are some hack-ish 
steps that we are trying to figure out, but for now this seems to work.

1. Change to this directory, `<sw4 home>/tools/sw4img`

2. Locate installed version of VisIt, and set the environment variables for VisIt top-level directory and plugin directory:

    `export VISITARCHHOME=/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64`  
    `export VISITPLUGININSTPRI=/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/plugins`

3. Run `xml2cmake` (located in `$VISITARCHHOME/bin`):

    `xml2cmake -clobber sw4.xml`  
    `cmake .`  
    `make clean`  
    `make`  
\
If you download [CMake](https://cmake.org), the executable will be in `/Applications/CMake.app/Contents/bin`

4. This may work, or it may produce build error with an incorrect, long path:

    `Scanning dependencies of target Isw4imgDatabase`   
    `[  5%] Building CXX object CMakeFiles/Isw4imgDatabase.dir/sw4imgPluginInfo.C.o`  
    ``make[2]: *** No rule to make target `/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/Projects/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/VisIt/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/Thirdparty/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/2.13RC/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/visit/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/zlib/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/1.2.7/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/i386-apple-darwin16_clang/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/lib/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/lib/libz.dylib', needed by `/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86_64/plugins/databases/libIsw4imgDatabase.dylib'.  Stop.``  
    `make[1]: *** [CMakeFiles/Isw4imgDatabase.dir/all] Error 2`  

5. The solution is to *edit* the cmake files in each of these CMakeFiles subdirectories:  

    `Esw4imgDatabase_par.dir`  
    `Esw4imgDatabase_ser.dir`  
    `Isw4imgDatabase.dir`  
    `Msw4imgDatabase.dir`  
\
In each one, there’s are two files:  

    `build.cmake`  
    `link.txt`  
\
In each of those, there are 1 or more references to the crazy long path, in front of `libz.dylib`, `libmpi.dylib`, and `libpmpi.dylib`. Fix each of those paths to point to the correct directory (the path given by `$VISITARCHHOME/lib`)

6. Then rebuild:  

    `> make clean`  
    `> make`  
    `...`  
    `[100%] Built target Msw4imgDatabase`  
\
And it should install the 4 new libraries in `$VISITPLUGININSTPRI/databases`:  

    `-rwxr-xr-x  1 user  admin  161408 Nov  6 16:18
libEsw4imgDatabase_par.dylib`  
    `-rwxr-xr-x  1 user  admin  157296 Nov  6 16:18
libEsw4imgDatabase_ser.dylib`  
    `-rwxr-xr-x  1 user  admin   28500 Nov  6 16:18 libIsw4imgDatabase.dylib`  
    `-rwxr-xr-x  1 user  admin  157336 Nov  6 16:18 libMsw4imgDatabase.dylib`  
    
## Building for LLNL LC version of VisIt 

Adding the plugin on LC has a similar problem as described in the MacOS section, namely the repeating of the library path. For example, after downloading the sw4img files to `~/.visit/linux-x86_64/plugins/databases/SW4` and performing step 3 listed above, the `build.make` file in `./CMakeFiles/Esw4imgDatabase_par.dir` has a line that looks like this:

`/g/g14/ford17/.visit/2.13.2/linux-x86_64-toss3/plugins/databases/libEsw4imgDatabase_par.so: /usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/usr/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/workspace/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/wsa/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/visit/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/visit/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/thirdparty_shared/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/2.13.0/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/toss3/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/zlib/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/1.2.7/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/linux-x86_64_gcc-4.9/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/lib/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/libz.so`  
\
It is your job to trim it (and the other build.make and link.txt files) to this:

`/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/lib/usr/gapps/visit/2.13.2/linux-x86_64-toss3/lib/libz.so`  
\
You can then run make in the SW4 directory and it should compile the plugin for use with LC hosted 3Dimg files.