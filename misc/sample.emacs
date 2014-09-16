;;; 
;;; Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
;;; Produced at the Lawrence Livermore National Laboratory. Written by 
;;; Jacob Schroder schroder2@llnl.gov, Rob Falgout falgout2@llnl.gov,
;;; Tzanio Kolev kolev1@llnl.gov, Ulrike Yang yang11@llnl.gov, 
;;; Veselin Dobrev dobrev1@llnl.gov, et al. 
;;; LLNL-CODE-660355. All rights reserved.
;;; 
;;; This file is part of XBraid. Email schroder2@llnl.gov on how to download. 
;;; 
;;; This program is free software; you can redistribute it and/or modify it under
;;; the terms of the GNU General Public License (as published by the Free Software
;;; Foundation) version 2.1 dated February 1999.
;;; 
;;; This program is distributed in the hope that it will be useful, but WITHOUT ANY
;;; WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
;;; PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
;;; License for more details.
;;; 
;;; You should have received a copy of the GNU Lesser General Public License along
;;; with this program; if not, write to the Free Software Foundation, Inc., 59
;;; Temple Place, Suite 330, Boston, MA 02111-1307 USA
;;;

;;;--------------------------------------------------------------------------
;;; 
;;; casc customizations for .emacs file. can be used with emacs or xemacs.
;;;
;;; to use this file just by itself, do
;;;     emacs  -l sample.emacs file_to_edit
;;;
;;; to get emacs to automatically use this file on startup, add this to
;;; your .emacs file.
;;;
;;; author: chandrika kamath
;;;         lawrence livermore national laboratory
;;;         livermore, ca 94551
;;;
;;; date: october 24, 1997
;;;
;;; modified:  1) Added the font lock for c++ mode.
;;;               Added indent-tabs-mode to force indentation with spaces
;;;               chandrika kamath, december 9, 1997.
;;;--------------------------------------------------------------------------

;;;
;;; C/C++ mode
;;;
;;; c-auto-newline automatically inserts a newline before and after a
;;; brace, and after a ";". To toggle this on/off, use C-x C-a. To 
;;; turn it off permanently, comment out the line that sets c-auto-newline.
;;; 
;;; the ellemtel style sets the opening brace on a new line. to set the
;;; brace at the end of the same line, use the stroustrup style.
;;;

(defun my-c-mode-common-hook ()
  (c-set-style "ellemtel")
  (setq c-tab-always-indent t)
  (setq c-basic-offset 3)
  (c-set-offset 'knr-argdecl-intro 0)
  (setq c-auto-newline t)
  (setq c-continued-statement-offset 3)
  (setq c-brace-offset 0)
  (c-set-offset 'inclass '+)
  (setq line-number-mode t)
  (setq indent-tabs-mode nil)
)

(add-hook 'c-mode-common-hook 'my-c-mode-common-hook)

;;;
;;; to turn on color in cc-mode (to set background color, see the manpages for
;;; emacs or xemacs)
;;;


