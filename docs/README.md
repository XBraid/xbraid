## Using Doxygen
<!--
 - Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 - Produced at the Lawrence Livermore National Laboratory.
 - This file is part of XBraid.  See file COPYRIGHT for details.
 -
 - XBraid is free software; you can redistribute it and/or modify it under the
 - terms of the GNU Lesser General Public License (as published by the Free
 - Software Foundation) version 2.1 dated February 1999.
 -->

To build the documentation, doxygen must be version 1.8 or greater.
\f$\chi\f$Braid documentation uses a 
[markdown](http://www.stack.nl/~dimitri/doxygen/manual/markdown.html) syntax
both in source file comments and in \*.md files.  

To make the documentation,

      $ make user_manual 
      $ acroread user_manual.pdf

or to make a more extensive reference manual for developers, 

      $ make developer_manual 
      $ acroread developer_manual.pdf

Developers can run doxygen from a precompiled binary, 
which may or may not work for your machine, 

      /usr/casc/hypre/braid/share/doxygen/bin/doxygen

or build doxygen from

      /usr/casc/hypre/braid/share/doxygen.tgz

- Compiling doxygen requires a number of dependencies
  like Bison, GraphViz and Flex.  Configure will tell 
  you what you're missing
- Unpack doxygen.tgz, then from the doxygen directory

      ./configure --prefix some_dir_in_your_path
      make
      make install

### Documentation Strategy
-  The doxygen comments are to be placed in the header files.
-  A sample function declaration using the documenation approach
   using markdown (including typesetting equations) is in braid.h
   for the function braid_Init()
-  A sample structure documentation is in _braid.h for _braid_Core_struct
-  Descriptors for files can also be added, as at the top of braid.h
-  The Doxygen manual is at 
   http://www.stack.nl/~dimitri/doxygen/manual/index.html

### \f$\chi\f$Braid Doxygen details

The user and developer manuals are ultimately produced by Latex.  The formatting 
of the manuals is configured according to the following.
-  docs/local_doxygen.sty           
  + Latex style file used
-  docs/user_manual_header.tex      
  + User manual title page and header info
-  docs/developer_manual_header.tex
  + Developer manual title page and header info
-  *.md                      
  + Any file ending in .md is extra documentation in markdown format, 
    like Introduction.md or the various Readme.md files in each directory.  
    This material can be read in plain-text or when it's compiled by Doxygen and Latex.
-  docs/user_manual.conf             
  + Doxygen configure file for the user manual
  + The FILE_NAMES tag is a filter to only include the user interface routines in braid.h
  + The INPUT tag orders the processing of the files and hence the section ordering
-  docs/reference_manual.conf       
  + Same as user_manual.conf, but the FILE_NAMES tag does not exclude any 
    file from processing.
- docs/img                         
  + Contains the images

-  To regenerate generic doxygen latex files, type
  
         $ doxygen -w latex header.tex footer.tex doxygen.sty doxy.conf

   If this is done, then the .conf file must be changed to use the new header file
   and to copy the local_doxygen.sty file to the latex directory.


