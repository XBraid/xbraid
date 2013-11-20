## Using Doxygen

To make the documentation:

      $ doxygen user_manual.conf
      $ cd user_manual
      $ make
      $ acroread refman.pdf

or to make a more extensive reference manual for developers, 

      $ doxygen reference_manual.conf
      $ cd reference_manual
      $ make
      $ acroread refman.pdf

Here is some background on doxygen
-  http://www.stack.nl/~dimitri/doxygen/manual/index.html
-  Developers should run doxygen from /usr/casc/hypre/warp/share/doxygen/bin/doxygen
-  The doxygen comments are to be placed in the header files.

-  A sample function declaration using the markdown approach
   to typesetting with equations is warp_Init() in warp.h
-  A sample structure is  _warp_Core_struct in _warp.h
-  Descriptors for files can also be added, as at the top of warp.h

-  The latex manuals are built according to 

   -  docs/local_doxygen.sty           :: Latex style file
   -  docs/user_manual_header.tex      :: User manual title page and header info
   -  docs/developer_manual_header.tex :: Developer manual title page and header info
   -  docs/Introduction.md             :: Extra material that goes at the front of the PDF
   -  docs/user_manual.conf            :: Only includes the user interface routines in warp.h
   -  docs/reference_manual.conf       :: Includes everything, and the kitchen sink
   -  docs/img                         :: Contains the images

-  To regenerate generic doxygen latex files, type
  
         $ doxygen -w latex header.tex footer.tex doxygen.sty doxy.conf

   The .conf file must then be changed to use the new header file
   and to copy the local_doxygen.sty file to the latex directory.


