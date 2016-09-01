"""
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of fortebiopkg.                                     !
   !                                                                          !
   !    Version 0.1.0                                                         !
   !                                                                          !
   !    Copyright (C) 2016                                                    !
   !    Boston AutoLytics, LLC                                                           !
   !                                                                          !
   !                                                                          !
   !--------------------------------------------------------------------------!
   ! AUTHORSHIP INFO:                                                         !
   !--------------------------------------------------------------------------!
   !                                                                          !
   ! MAIN AUTHOR:   R. Paul Nobrega                                           !
   !                                                                          !
   !--------------------------------------------------------------------------!
   File Description:
   ================
   This is the main package import file for fortebiopkg.
"""

import fileio as fileio
import datacorrection as datacorrection
import datafitting as datafitting


__all__ = ['fileio', 'datacorrection', 'datafitting']



