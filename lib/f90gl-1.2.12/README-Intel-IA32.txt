As a courtesy to users of Intel Visual Fortran, we provide a copy of the f90GL sources, samples and build scripts, along with a precompiled version of f90GL and the GLUT library it depends on.

They have been built using Intel(R) Visual Fortran 9.1 and Intel(R) C++ Compiler 9.1 for IA-32.

To use this with Intel Visual Fortran:

1.  Copy lib-IA32\glut32.dll to your Windows System32 folder or some folder on your PATH
2.  Modify the INCLUDE/USE path to include the lib-IA32 folder from this package.  From
    the IDE, this is Tools..Options..Intel(R) Fortran..Project Directories..Include.
    If using the command line, use /include or modify the INCLUDE environment variable
3.  If using the IDE, add the .LIB files from the lib-IA32 folder to your project
    (Project..Add Existing Item).  If using the command line, you will need to
    name these libraries on your link command.

If you wish to use the build scripts in this folder, copy the contents of the lib-IA32 folder
to the lib folder.

Folder "examples" contains example programs.  Batch file mf8nio.bat will build for Intel
Visual Fortran.  To build examples, start a Fortran Command Prompt window (Start..Programs..
Intel(R) Software Development Tools..Intel(R) Fortran Compiler x.x..Build Environment for
Fortran IA-32 Applications) and type the command:

	mf8nio all

You can build individual examples by naming the example on the command line, for example:

	mf8nio blender

These files are provided as a courtesy of Intel and are not supported by Intel. Please
do not contact Intel Fortran support for assistance with f90GL or GLUT.  f90gl was produced by an agency of the U.S. Government, and is not subject to copyright in the United States. For information on 
f90GL, please see http://math.nist.gov/f90gl/ 

The GLUT library is copyrighted by Mark J. Kilgard. GLUT is not in the public domain, but it is freely distributable without licensing fees. GLUT is provided without gurantee or warrantee expressed or implied. For more information about GLUT, please see http://www.opengl.org/resources/libraries/glut/

Intel, the Intel logo, Itanium, Xeon, Pentium and PCA are trademarks or registered 
trademarks of Intel Corporation or its subsidiaries in the United States and other countries. 
*Other brands and names are the property of their respective owners. 

This text Copyright © 2007 Intel Corporation. All Rights Reserved. 