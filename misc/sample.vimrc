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
