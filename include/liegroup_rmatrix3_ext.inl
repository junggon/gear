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

inline RMatrix Ad(const SE3 &T, const RMatrix &J)
{
	se3 S;
	RMatrix re(J.RowSize(), J.ColSize());
	for (int i=0; i<J.ColSize(); i++)
	{
		S = se3(J[6*i], J[6*i+1], J[6*i+2], J[6*i+3], J[6*i+4], J[6*i+5]);
		//put_se3_to_matrix(re, Ad(T,S), 6*i);
		matSet(&re[6*i], Ad(T,S).GetArray(), 6);
	}
	return re;
}

inline RMatrix dAd(const SE3 &T, const RMatrix &J)
{
	dse3 S;
	RMatrix re(J.RowSize(), J.ColSize());
	for (int i=0; i<J.ColSize(); i++)
	{
		S = dse3(J[6*i], J[6*i+1], J[6*i+2], J[6*i+3], J[6*i+4], J[6*i+5]);
		//put_dse3_to_matrix(re, dAd(T,S), 6*i);
		matSet(&re[6*i], dAd(T,S).GetArray(), 6);
	}
	return re;
}

inline RMatrix ad(const se3 &S, const RMatrix &J)
{
	se3 S2;
	RMatrix re(J.RowSize(), J.ColSize());
	for (int i=0; i<J.ColSize(); i++)
	{
		S2 = se3(J[6*i], J[6*i+1], J[6*i+2], J[6*i+3], J[6*i+4], J[6*i+5]);
		//put_se3_to_matrix(re, ad(S,S2), 6*i);
		matSet(&re[6*i], ad(S,S2).GetArray(), 6);
	}
	return re;
}

inline RMatrix ad(const RMatrix &S, const RMatrix &J)
{
	se3 SS(S[0], S[1], S[2], S[3], S[4], S[5]);
	return ad(SS, J);
}

inline RMatrix dad(const se3 &S, const RMatrix &J)
{
	dse3 S2;
	RMatrix JJ(J.RowSize(), J.ColSize());
	for (int i=0; i<J.ColSize(); i++)
	{
		S2 = dse3(J[6*i], J[6*i+1], J[6*i+2], J[6*i+3], J[6*i+4], J[6*i+5]);
		//put_dse3_to_matrix(JJ, dad(S,S2), 6*i);
		matSet(&JJ[6*i], dad(S,S2).GetArray(), 6);
	}
	return JJ;
}

inline RMatrix Cross(const Vec3 &a, const RMatrix &B)
{
	Vec3 b;
	RMatrix re(B.RowSize(), B.ColSize());
	for (int i=0; i<B.ColSize(); i++) {
		b = Vec3(B[3*i], B[3*i+1], B[3*i+2]);
		//put_Vec3_to_matrix(re, Cross(a,b), 3*i);
		matSet(&re[3*i], Cross(a,b).GetArray(), 3);
	}
	return re;
}

