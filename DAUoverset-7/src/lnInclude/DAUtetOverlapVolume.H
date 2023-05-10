/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2014 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::DAUtetOverlapVolume

Description
    Calculates the overlap volume of two cells using tetrahedral decomposition

SourceFiles
    DAUtetOverlapVolume.C
    DAUtetOverlapVolumeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef DAUtetOverlapVolume_H
#define DAUtetOverlapVolume_H

#include "FixedList.H"
#include "labelList.H"
#include "DAUtreeBoundBox.H"
#include "Tuple2.H"
#include "tetrahedron.H"
#include "tetPointRef.H" // by SBLEE
#include "tetPoints.H" // by SBLEE

namespace Foam
{

class primitiveMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class DAUtetOverlapVolume Declaration
\*---------------------------------------------------------------------------*/

class DAUtetOverlapVolume
{
    // Private classes

        typedef FixedList<tetPoints, 200> tetIntersectionList;

        //- tetPoints handling : sum resulting volumes
        class sumMomentOp
        {
        public:
            Tuple2<scalar, point> vol_;

            inline sumMomentOp()
            :
                vol_(0.0, Zero)
            {}

            inline void operator()(const tetPoints& tet)
            {
                const tetPointRef t(tet.tet());
                scalar tetVol = t.mag();
                vol_.first() += tetVol;
                vol_.second() += (tetVol*t.centre());
            }
        };

        class dummyOp
        {
        public:
            inline void operator()(const tetPoints& tet)
            {}
        };

        class sumVolOp
        {
        public:
            scalar vol_;
            inline sumVolOp()
            :
                vol_(0.0)
            {}
            inline void operator()(const tetPoints& tet)
            { vol_+=tet.tet().mag(); }
        };

        class storeOp
        {
            tetIntersectionList& tets_;
            label& nTets_;

        public:
            inline storeOp(tetIntersectionList& tets, label& nTets)
            :
                tets_(tets),
                nTets_(nTets)
            {}
            inline void operator()(const tetPoints& tet)
            { tets_[nTets_++]=tet; }
        };

        //- tetPoints combining : check for overlap
        class hasOverlapOp
        {
        public:

            const scalar threshold_;
            sumVolOp iop_;
            bool ok_;

            inline hasOverlapOp(const scalar threshold)
            :
                threshold_(threshold),
                iop_(),
                ok_(false)
            {}

            //- Overlap two tets
            inline bool operator()(const tetPoints& A, const tetPoints& B)
            {
                tetTetOverlap<sumVolOp>(A, B, iop_);
                ok_ = (iop_.vol_ > threshold_);
                return ok_;
            }
        };

        //- tetPoints combining : sum overlap volume
        class sumOverlapOp
        {
        public:

            sumVolOp iop_;

            inline sumOverlapOp()
            :
                iop_()
            {}

            //- Overlap two tets
            inline bool operator()(const tetPoints& A, const tetPoints& B)
            {
                tetTetOverlap<sumVolOp>(A, B, iop_);
                return false;
            }
        };

        //- tetPoints combining : sum overlap volume
        class sumOverlapMomentOp
        {
        public:

            sumMomentOp iop_;

            inline sumOverlapMomentOp()
            :
                iop_()
            {}

            //- Overlap two tets
            inline bool operator()(const tetPoints& A, const tetPoints& B)
            {
                tetTetOverlap<sumMomentOp>(A, B, iop_);
                return false;
            }
        };


    // Private member functions

        //- Tet overlap calculation
        template<class tetPointsOp>
        static void tetTetOverlap
        (
            const tetPoints& tetA,
            const tetPoints& tetB,
            tetPointsOp& insideOp
        );

        //- Cell overlap calculation
        template<class tetsOp>
        static void cellCellOverlapMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,
            const primitiveMesh& meshB,
            const label cellBI,
            const DAUtreeBoundBox& cellBbB,
            tetsOp& combineTetsOp
        );

        //- Return a const treeBoundBox
        static DAUtreeBoundBox pyrBb
        (
            const pointField& points,
            const face& f,
            const point& fc
        );


public:

    //- Runtime type information
    ClassName("DAUtetOverlapVolume");


    // Constructors

        //- Null constructor
        DAUtetOverlapVolume();


    // Public members

        //- Return a list of cells in meshA which overlaps with cellBI in
        // meshB
        labelList overlappingCells
        (
            const polyMesh& meshA,
            const polyMesh& meshB,
            const label cellBI
        ) const;

        //- Return true if overlap volume is greater than threshold
        bool cellCellOverlapMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,
            const primitiveMesh& meshB,
            const label cellBI,
            const DAUtreeBoundBox& cellBbB,
            const scalar threshold = 0.0
        ) const;

        //- Calculates the overlap volume
        scalar cellCellOverlapVolumeMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,

            const primitiveMesh& meshB,
            const label cellBI,
            const DAUtreeBoundBox& cellBbB
        ) const;

        //- Calculates the overlap volume and moment
        Tuple2<scalar, point> cellCellOverlapMomentMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,

            const primitiveMesh& meshB,
            const label cellBI,
            const DAUtreeBoundBox& cellBbB
        ) const;

        inline static point planeIntersection
        (
            const FixedList<scalar,4>& d,
            const tetPoints& t,
            const label negI,
            const label posI
        )
        {
            return (d[posI]*t[negI]-d[negI]*t[posI])/(-d[negI]+d[posI]);
        }

        template<class TetOp>
        inline static void decomposePrism
        (
            const FixedList<point, 6>& points,
            TetOp& op
        )
        {
            op(tetPoints(points[1], points[3], points[2], points[0]));
            op(tetPoints(points[1], points[2], points[3], points[4]));
            op(tetPoints(points[4], points[2], points[3], points[5]));
        }
        template<class AboveTetOp, class BelowTetOp>
        inline static void tetSliceWithPlane
        (
            const plane& pl,
            const tetPoints& tet,
            AboveTetOp& aboveOp,
            BelowTetOp& belowOp
        )
        {
            FixedList<scalar, 4> d;
            label nPos = 0;
            forAll(tet, i)
            {
                d[i] = ((tet[i]-pl.refPoint()) & pl.normal());
                if (d[i] > 0) nPos++;
            }
            if (nPos == 4) aboveOp(tet);
            else if (nPos == 3)
            {
                label i0 = -1;
                forAll(d, i)
                {
                    if (d[i]<=0)
                    {
                        i0 = i;
                        break;
                    }
                }
                label i1 = d.fcIndex(i0);
                label i2 = d.fcIndex(i1);
                label i3 = d.fcIndex(i2);

                point p01 = planeIntersection(d, tet, i0, i1);
                point p02 = planeIntersection(d, tet, i0, i2);
                point p03 = planeIntersection(d, tet, i0, i3);
                if (i0 == 0 || i0 == 2)
                {
                    tetPoints t(tet[i0], p01, p02, p03);
                    belowOp(t);
                    FixedList<point,6> p;
                    p[0]=tet[i1];
                    p[1]=tet[i3];
                    p[2]=tet[i2];
                    p[3] = p01;
                    p[4] = p02;
                    p[5] = p03;
                    decomposePrism(p, aboveOp);
                }
                else
                {
                    tetPoints t(p01, p02, p03, tet[i0]);
                    belowOp(t);
                    FixedList<point,6> p;
                    p[0]=tet[i3];
                    p[1]=tet[i1];
                    p[2]=tet[i2];
                    p[3] = p03;
                    p[4] = p01;
                    p[5] = p02;
                    decomposePrism(p, aboveOp);
                }
            }
            else if (nPos == 2)
            {
                label pos0=-1;
                label pos1=-1;
                forAll(d, i)
                {
                    if (d[i] > 0)
                    {
                        if (pos0 == -1) pos0 = i;
                        else pos1 = i;
                    }
                }
                const edge posEdge(pos0, pos1);
                if (posEdge == edge(0, 1))
                {
                    point p02 = planeIntersection(d, tet, 0, 2);
                    point p03 = planeIntersection(d, tet, 0, 3);
                    point p12 = planeIntersection(d, tet, 1, 2);
                    point p13 = planeIntersection(d, tet, 1, 3);
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[0];
                        p[1] = p02;
                        p[2] = p03;
                        p[3] = tet[1];
                        p[4] = p12;
                        p[5] = p13;
                        decomposePrism(p, aboveOp);
                    }
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[2];
                        p[1] = p02;
                        p[2] = p12;
                        p[3] = tet[3];
                        p[4] = p03;
                        p[5] = p13;
                        decomposePrism(p, belowOp);
                    }
                }
                else if (posEdge == edge(1, 2))
                {
                    point p01 = planeIntersection(d, tet, 0, 1);
                    point p13 = planeIntersection(d, tet, 1, 3);
                    point p02 = planeIntersection(d, tet, 0, 2);
                    point p23 = planeIntersection(d, tet, 2, 3);
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[1];
                        p[1] = p01;
                        p[2] = p13;
                        p[3] = tet[2];
                        p[4] = p02;
                        p[5] = p23;
                        decomposePrism(p, aboveOp);
                    }
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[3];
                        p[1] = p23;
                        p[2] = p13;
                        p[3] = tet[0];
                        p[4] = p02;
                        p[5] = p01;
                        decomposePrism(p, belowOp);
                    }
                }
                else if (posEdge == edge(2, 0))
                {
                    point p01 = planeIntersection(d, tet, 0, 1);
                    point p03 = planeIntersection(d, tet, 0, 3);
                    point p12 = planeIntersection(d, tet, 1, 2);
                    point p23 = planeIntersection(d, tet, 2, 3);
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[2];
                        p[1] = p12;
                        p[2] = p23;
                        p[3] = tet[0];
                        p[4] = p01;
                        p[5] = p03;
                        decomposePrism(p, aboveOp);
                    }
                   {
                        FixedList<point, 6> p;
                        p[0] = tet[1];
                        p[1] = p12;
                        p[2] = p01;
                        p[3] = tet[3];
                        p[4] = p23;
                        p[5] = p03;
                        decomposePrism(p, belowOp);
                    }
                }
                else if (posEdge == edge(0, 3))
                {
                    point p01 = planeIntersection(d, tet, 0, 1);
                    point p02 = planeIntersection(d, tet, 0, 2);
                    point p13 = planeIntersection(d, tet, 1, 3);
                    point p23 = planeIntersection(d, tet, 2, 3);
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[3];
                        p[1] = p23;
                        p[2] = p13;
                        p[3] = tet[0];
                        p[4] = p02;
                        p[5] = p01;
                        decomposePrism(p, aboveOp);
                    }
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[2];
                        p[1] = p23;
                        p[2] = p02;
                        p[3] = tet[1];
                        p[4] = p13;
                        p[5] = p01;
                        decomposePrism(p, belowOp);
                    }
                }
                else if (posEdge == edge(1, 3))
                {
                    point p01 = planeIntersection(d, tet, 0, 1);
                    point p12 = planeIntersection(d, tet, 1, 2);
                    point p03 = planeIntersection(d, tet, 0, 3);
                    point p23 = planeIntersection(d, tet, 2, 3);
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[1];
                        p[1] = p12;
                        p[2] = p01;
                        p[3] = tet[3];
                        p[4] = p23;
                        p[5] = p03;
                        decomposePrism(p, aboveOp);
                    }
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[2];
                        p[1] = p12;
                        p[2] = p23;
                        p[3] = tet[0];
                        p[4] = p01;
                        p[5] = p03;
                        decomposePrism(p, belowOp);
                    }
                }
                else if (posEdge == edge(2, 3))
                {
                    point p02 = planeIntersection(d, tet, 0, 2);
                    point p12 = planeIntersection(d, tet, 1, 2);
                    point p03 = planeIntersection(d, tet, 0, 3);
                    point p13 = planeIntersection(d, tet, 1, 3);
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[2];
                        p[1] = p02;
                        p[2] = p12;
                        p[3] = tet[3];
                        p[4] = p03;
                        p[5] = p13;
                        decomposePrism(p, aboveOp);
                    }
                    {
                        FixedList<point, 6> p;
                        p[0] = tet[0];
                        p[1] = p02;
                        p[2] = p03;
                        p[3] = tet[1];
                        p[4] = p12;
                        p[5] = p13;
                        decomposePrism(p, belowOp);
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << "Missed edge:" << posEdge
                        << abort(FatalError);
                }

            }
            else if (nPos == 1)
            {
                label i0 = -1;
                forAll(d, i)
                {
                    if (d[i] > 0)
                    {
                        i0 = i;
                        break;
                    }
                }
                label i1 = d.fcIndex(i0);
                label i2 = d.fcIndex(i1);
                label i3 = d.fcIndex(i2);

                point p01 = planeIntersection(d, tet, i0, i1);
                point p02 = planeIntersection(d, tet, i0, i2);
                point p03 = planeIntersection(d, tet, i0, i3);
                if (i0 == 0 || i0 == 2)
                {
                    tetPoints t(tet[i0], p01, p02, p03);
                    aboveOp(t);
                    FixedList<point, 6> p;
                    p[0] = tet[i1];
                    p[1] = tet[i3];
                    p[2] = tet[i2];
                    p[3] = p01;
                    p[4] = p03;
                    p[5] = p02;
                    decomposePrism(p, belowOp);
                }
                else
                {
                    tetPoints t(p01, p02, p03, tet[i0]);
                    aboveOp(t);
                    FixedList<point, 6> p;
                    p[0] = tet[i3];
                    p[1] = tet[i1];
                    p[2] = tet[i2];
                    p[3] = p03;
                    p[4] = p01;
                    p[5] = p02;
                    decomposePrism(p, belowOp);
                }
            }
            else    // nPos == 0
            {
                belowOp(tet);
            }
        }

        template<class AboveTetOp, class BelowTetOp>
        inline static void sliceWithPlane
        (
            tetPointRef tPR,
            const plane& pl,
            AboveTetOp& aboveOp,
            BelowTetOp& belowOp
        )
        {
            tetSliceWithPlane(pl, tetPoints(tPR.a(), tPR.b(), tPR.c(), tPR.d()), aboveOp, belowOp);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DAUtetOverlapVolumeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
