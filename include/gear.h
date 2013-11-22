// -------------------------------------------------------------------------------
// Copyright (c) 2012, Junggon Kim
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met: 
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution. 
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -------------------------------------------------------------------------------

//================================================================================
//         GEAR: Geometric Engine for Articulated Rigid-body system
// 
//                                                               junggon@gmail.com
//================================================================================

#include "greal.h"
#include "rmatrix3j.h"
#include "liegroup.h"
#include "liegroup_rmatrix3_ext.h"
#include "gcoordinate.h"
#include "gelement.h"
#include "gbody.h"
#include "gjoint.h"
#include "gjoint_fixed.h"
#include "gjoint_prismatic.h"
#include "gjoint_planar.h"
#include "gjoint_translational.h"
#include "gjoint_revolute.h"
#include "gjoint_universal.h"
#include "gjoint_spherical.h"
#include "gjoint_free.h"
#include "gjoint_composite.h"
#include "gforce.h"
#include "gspringdamper.h"
#include "gconstraint.h"
#include "gconstraint_jointloop.h"
#include "gsystem.h"
#include "gsystem_constrained.h"
#include "gsystem_ik.h"


