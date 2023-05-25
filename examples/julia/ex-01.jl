#=BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory.
 * 
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER=#

# main ex-01
# include("XBraid.jl")
include("../../braid/braid.jl/XBraid.jl")
using .XBraid
using MPI

MPI.Init()
comm = MPI.COMM_WORLD

#=
using 1-element vector here since standard Float64
values are immutable and don't have stable addresses.
=#
function my_init(app, t)
    t == 0.0 && return [1.0]
    return [0.456]
end

# This must mutate u in place
function my_step!(app, status, u, ustop, tstart, tstop)
    # backward Euler
    u[] = 1.0 / (1.0 + tstop - tstart) * u[]
end

#=
called by the solver to give access to 
solution values.
=#
function my_access(app, status, u)
    XBraid.status_GetWrapperTest(status) && return
    t = XBraid.status_GetT(status)
    ti = XBraid.status_GetTIndex(status)
    print("t: $(t[]),\tu: $(u[])\n")
end

#=
all other usual wrapper functions are either not required, 
or are set to reasonable defaults for operating on julia Arrays. 😎
=#

test = false
if test
    test_app = XBraid.BraidApp(nothing, comm, my_step!, my_init, my_access)

    open("ex_01.test.out", "w") do file
        cfile = Libc.FILE(file)
        XBraid.testInitAccess(test_app, 0.0, cfile)
        XBraid.testClone(test_app, 0.0, cfile)
        XBraid.testSum(test_app, 0.0, cfile)
        XBraid.testSpatialNorm(test_app, 0.0, cfile)
        XBraid.testBuf(test_app, 0.0, cfile)
    end
    MPI.Barrier(comm)
end

ntime = 10
tstart = 0.0
tstop = tstart + ntime / 2.0;
core = XBraid.Init(comm, tstart, tstop, ntime, my_step!, my_init, my_access)

XBraid.SetPrintLevel(core, 2)
XBraid.SetMaxLevels(core, 2)
XBraid.SetAbsTol(core, 1.e-6)
XBraid.SetCFactor(core, -1, 2)

XBraid.Drive(core)

# no need for braid_Destroy(core), julia will take care of it :)
MPI.Finalize()