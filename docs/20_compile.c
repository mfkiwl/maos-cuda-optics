/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
   \page page20_compile  Compile the Code

   \tableofcontents

   This step can be skipped if the user choose to use the precompiled binary files,
   which are available for all the popular plateforms: Windows, Mac, and
   Linux. Please notice that the job scheduling and monitoring utilities are not
   included in the released binaries (the binaries may not be up to data).

   \section sect-requirement Requirements

   - C99 compliant compiler: GCC, Intel ICC, or clang.

   - (Optional) GTK+. For drawdaemon and monitor.

   - FFTW version 3. Will download from MAOS site if not available in the system.
    
   - Optimized blas and lapack. Blas and lapack can usually be installed through
   the linux distribution official repository. Will download from MAOS site if
   not available in the system.
    
   - (Optional) libwebsocket. For web based job monitoring. It requires cmake to compile.
    
   \section sect-prep Preparing the folders and compiling
    
   We recommend using three different folders to 1) store the source tree, 2)
   compile the code, and 3) run the simulation.

   The source tree is where the source is extracted. Obtain the software in the
   form of the compressed tar ball: maos_version.tar.gz. Extract it into some
   directory and name the directory src_dir. Example:

   \verbatim
   cd ~/work/programming
   tar zxvf maos_version.tar.gz 
   cd maos_version
   export src_dir=$(pwd)
   \endverbatim

   If you choose to use always up to date git version instead, an extra step of running autogen.sh is needed
   \verbatim
   cd ~/work/programming
   git clone https://github.com/lianqiw/maos.git maos_git
   cd maos_git
   ./autogen.sh
   export src_dir=$(pwd)
   \endverbatim

   Next we create another folder, where we are going to compile the code. example:

   \verbatim
   cd ~/work/maos
   mkdir comp_optim && cd comp_optim
   $src_dir/configure
   make -j4 #this will compile using 4 threads.
   \endverbatim

   For clusters that use job management tools, disable the built-in scheduler by using
   \verbatim
   $src_dir/configure --disable-scheduler
   \endverbatim

   The compiled executable is maos in the sub-folder “bin” of the compiling
   folder. You do not have to do "make install" to run the simulations.

   \subsection sect-compiler Compiler options 

   The default compiler for the system is used with optimization flags. There
   are a few optional switches or additional components that can be enabled when
   calling configure.
    
   To use a different compiler:
   \verbatim
   $src_dir/configure CC=icc    #ICC compiler
   $src_dir/configure CC=clang  #clang compiler
   $src_dir/configure CC=gcc4.5 #use gcc4.5
   $src_dir/configure --enable-debug #enable debugging
   \endverbatim

   \subsection sect-cuda GPU acceleration

   To enable GPU computing, first install nvidia graphics driver and cuda
   software development toolkit (SDK). Issue the following command to ensure
   graphics configuration is correct: \c nvidia-smi. You should see all your
   nvidia graphics gpus are listed.
    
   The configure can automatically detect cuda if nvcc is in the path, or if
   cuda is installed in /usr/local/cuda. For other locations, you
   need to pass the cuda directory to configure:

   \verbatim $src_dir/configure --with-cuda=cuda_dir \endverbatim
    
   When GPU compute is compiled in maos and Nvidia GPUs are available, maos
   will make use of all the GPUs to run simulation. If that is not desired, the
   GPU compute can be disabled at run time by passing -g-1 to maos. Different
   set of GPU can also be selected by passing one or multiple -gn: 
    
   \verbatim 
   $maos -g-1    #disable GPU compute 
   $maos -g0 -g2 #Only use gpu 0 and 2
   $maos -G2     #Automatically choose 2 gpus.
   \endverbatim
    
   \subsection sect-mex Matlab Mex Routines

   To compile mex routines for matlab, pass the matlab directory to configure if it is not in \c $PATH:
    
   \verbatim
   $src_dir/configure --with-matlab=matlab_dir
   \endverbatim

   The compiled mex routines are in \c mex folder in the compilation
   directory. The \c read and \c write are for handling .bin and .fits files
   used by maos. The \c aolib exports many useful routines for use in
   matlab, including maos itself.
   
    
    
 

   \subsection sect-mac-gtk Installing GTK+ in MAC OS and Compile Monitor, Drawdaemon

   It is simpliest to install gtk using macport or homebrew. For manual install,
   follow the instructions on page gtk-osx.sourceforge.net/ to install gtk+ for
   osx. The environmental variables are set by jhbuild, here we take important
   ones and put them in our .bash_profile to avoid running jhbuild each
   time. Put the following in your .bash_profile. Adjust the path accordingly if
   you changed gtk-osx installing options.

   First, if xcode is not installed in full, make sure the following packages are installed

   \verbatim
   DeveloperToolsCLI.pkg
   gcc4.x.pkg
   MACOSX10.x.pkg
   DevSDK.pkg
   OpenGLSDK.pkg
   DeveloperToolsSystemSupport
   \endverbatim
   Then download and install pkg-config from pkgconfig.freedesktop.org/releases/     


   \verbatim
   export PKG_CONFIG_PATH="/usr/lib/pkgconfig:$HOME/gtk/inst/lib/pkgconfig:$HOME/gtk/inst/share/pkgconfig"
   export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/gtk/inst/lib"
   export PATH="$PATH:$HOME/gtk/inst/bin"
   export C_INCLUDE_PATH="$HOME/gtk/inst/include"
   \endverbatim
   Now rerun MAOS configure and make. Monitor and drawdaemon should appear in bin folder now.

   To make the app bundles. 
   First, install quartz engine and ige-mac-buildler. 
   sourceforge.net/apps/trac/gtk-osx/wiki/Bundle 
   jhbuild build gtk-quartz-engine
    
   Second, symbolic link bin/drawdaemon and bin/monitor to $HOME/gtk/inst/bin
   And cd to scripts/bundle folder, run build.sh. .app files should appear in your desktop.
    
   \subsection sect-monitor Using the job monitor
    
   When GTK+ libraries are present in the system, additional executives will
   be compiled, \c drawdaemon, \c drawres, \c drawbin and \c monitor in the \c bin folder. 

   - The plotting utility \c drawdaemon launched automatically when plotting
   commands are issued in the software.

   - The monitoring utility \c monitor can be used to monitor jobs in this and linked machines.

   - The \c drawres can be used to automatically plot the results folder. 

   - The \c drawbin can be used to automatically plot the OPDs or maps.


   In order to use monitor to monitor jobs in different machines, follow these steps:

   Create a folder .aos in your home folder, create two files in there. Put
   these files in all the machines you want to monitor. 

   - \c ~/.aos/hosts Each line contains a machine's hostname you want to monitor. It also supports the
   full format \c prettyname=full_hostname:port 

   - \c ~/.aos/port contains the port you want the \c scheduler to bind to so that the
   monitor can connect to it.  Any non-previledged port numbers are fine (like
   10000).  Different users should have different port number.  The port hust to
   be the same on all the machines. Make sure the firewall allows them.

*/
