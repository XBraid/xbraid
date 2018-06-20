"" Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
"" Produced at the Lawrence Livermore National Laboratory. Written by 
"" Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
"" Dobrev, et al. LLNL-CODE-660355. All rights reserved.
"" 
"" This file is part of XBraid. For support, post issues to the XBraid Github page.
"" 
"" This program is free software; you can redistribute it and/or modify it under
"" the terms of the GNU General Public License (as published by the Free Software
"" Foundation) version 2.1 dated February 1999.
"" 
"" This program is distributed in the hope that it will be useful, but WITHOUT ANY
"" WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
"" PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
"" License for more details.
"" 
"" You should have received a copy of the GNU Lesser General Public License along
"" with this program; if not, write to the Free Software Foundation, Inc., 59
"" Temple Place, Suite 330, Boston, MA 02111-1307 USA
""
 

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
