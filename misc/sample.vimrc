"" Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
"" Produced at the Lawrence Livermore National Laboratory.
"" This file is part of WARP.  See file COPYRIGHT for details.
""
"" WARP is free software; you can redistribute it and/or modify it under the
"" terms of the GNU Lesser General Public License (as published by the Free
"" Software Foundation) version 2.1 dated February 1999.
 

" Add these lines to your .vimrc file to mimic the ellemtel C style 
" given by emacs

" Syntax highlighting
:syn on         
" Show line numbers
:set nu
" Autoindent
:set ai
" Tab and shift settings
:set tabstop=3
:set shiftwidth=3
:set expandtab
" Choose your color scheme
:colorscheme desert
" Automatically do language dependent indenting
filetype plugin indent on
" Set 'textwidth' to 78 characters.
autocmd FileType text setlocal textwidth=78
