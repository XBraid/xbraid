## Using Doxygen

To build the documentation, doxygen must be version 1.8 or greater.
Warp documentation uses a 
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

      /usr/casc/hypre/warp/share/doxygen/bin/doxygen

or build doxygen from

      /usr/casc/hypre/warp/share/doxygen.tgz

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
   using markdown (including typesetting equations) is in warp.h
   for the function warp_Init()
-  A sample structure documentation is in _warp.h for _warp_Core_struct
-  Descriptors for files can also be added, as at the top of warp.h
-  The Doxygen manual is at 
   http://www.stack.nl/~dimitri/doxygen/manual/index.html

### Warp Doxygen details
The latex manuals are built according to 
-  docs/local_doxygen.sty           
  + Latex style file used
-  docs/user_manual_header.tex      
  + User manual title page and header info
-  docs/developer_manual_header.tex
  + Developer manual title page and header info
-  filename.md                      
  + Extra material in markdown format, like Abstract.md and Introduction.md
-  docs/user_manual.conf             
  + Doxygen configure file for the user manual
  + The FILE_NAMES tag filters to only include the user interface routines in warp.h
  + The INPUT tag orders the processing of the files and hence the section ordering
-  docs/reference_manual.conf       
  + Same as user_manual.conf, out the FILE_NAMES tag doesn't exclude anything
- docs/img                         
  + Contains the images

-  To regenerate generic doxygen latex files, type
  
         $ doxygen -w latex header.tex footer.tex doxygen.sty doxy.conf

   The .conf file must then be changed to use the new header file
   and to copy the local_doxygen.sty file to the latex directory.


