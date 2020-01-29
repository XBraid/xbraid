<!--
  - Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
  - Produced at the Lawrence Livermore National Laboratory.
  - 
  - This file is part of XBraid. For support, post issues to the XBraid Github page.
  - 
  - This program is free software; you can redistribute it and/or modify it under
  - the terms of the GNU General Public License (as published by the Free Software
  - Foundation) version 2.1 dated February 1999.
  - 
  - This program is distributed in the hope that it will be useful, but WITHOUT ANY
  - WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
  - PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
  - License for more details.
  - 
  - You should have received a copy of the GNU Lesser General Public License along
  - with this program; if not, write to the Free Software Foundation, Inc., 59
  - Temple Place, Suite 330, Boston, MA 02111-1307 USA
 -->

# File naming conventions

User interface routines in braid begin with `braid_` and all other internal
non-user routines begin with `_braid_`.  This helps to prevent name clashes when
working with other libraries and helps to clearly distinguish user routines that
are supported and maintained.

To keep things somewhat organized, all user header files and implementation
files should have names that begin with `braid`, for example, `braid.h`,
`braid.c`, `braid_status.c`, ...  There should be no user interface prototypes
or implementations that appear elsewhere.

Note that it is okay to include internal prototypes and implementations in these
user interface files when it makes sense (say, as supporting routines), but this
should generally be avoided.

