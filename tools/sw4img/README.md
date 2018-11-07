# SW4 Visit plugin - BETA instructions

This README file contains vague instructions on building a Visit plugin
that can read SW4 files, 2D `\*.sw4img` and 3D or "volimage" `\*.3Dimg` files.

**NOTE**: This plugin is not supported by the Visit or SW4 team in any substantial
way, but it's being provided as a potentially useful tool for visualization
of SW4 data.

## Building against installed version of Visit on MacOS
The following steps worked on MacOS 10.13.6. There are some hack-ish 
steps that we are trying to figure out, but for now this seems to work.

1. Change to this directory, `<sw4 home>/tools/sw4img`

2. Locate installed version of Visit, and set the environment variables for Visit top-level directory and plugin directory:

    `export VISITARCHHOME=/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86\_64`  
    `export VISITPLUGININSTPRI=/Applications/VisIt.app/Contents/Resources/2.13.2/darwin-x86\_64/plugins`

3. Run `xml2cmake` (located in `$VISITARCHHOME/bin`):

    `xml2cmake -clobber sw4.xml`  
    `cmake .`  
    `make clean`  
    `make`

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
In each one, there’s are two files:  

    `build.cmake`  
    `link.txt`  
In each of those, there are 1 or more references to the crazy long path, in front of `libz.dylib`, `libmpi.dylib`, and `libpmpi.dylib`. Fix each of those paths to point to the correct directory (the path given by `$VISITARCHHOME/lib`)

6. Then rebuild:  

    `> make clean`  
    `> make`  
    `...`  
    `[100%] Built target Msw4imgDatabase`  
And it should install the 4 new libraries in `$VISITPLUGININSTPRI/databases`:  

    `-rwxr-xr-x  1 user  admin  161408 Nov  6 16:18
libEsw4imgDatabase_par.dylib`  
    `-rwxr-xr-x  1 user  admin  157296 Nov  6 16:18
libEsw4imgDatabase_ser.dylib`  
    `-rwxr-xr-x  1 user  admin   28500 Nov  6 16:18 libIsw4imgDatabase.dylib`  
    `-rwxr-xr-x  1 user  admin  157336 Nov  6 16:18 libMsw4imgDatabase.dylib`  
