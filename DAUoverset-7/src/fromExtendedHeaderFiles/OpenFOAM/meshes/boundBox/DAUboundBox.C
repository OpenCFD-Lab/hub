/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "DAUboundBox.H"
#include "PstreamReduceOps.H"
#include "tmp.H"
#include "plane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::DAUboundBox Foam::DAUboundBox::greatBox
(
    point::uniform(-ROOTVGREAT),
    point::uniform(ROOTVGREAT)
);

const Foam::DAUboundBox Foam::DAUboundBox::invertedBox
(
    point::uniform(ROOTVGREAT),
    point::uniform(-ROOTVGREAT)
);

/*
const Foam::faceList Foam::DAUboundBox::faces
({
    // Point and face order as per hex cellmodel
    face({0, 4, 7, 3}),  // 0: x-min, left
    face({1, 2, 6, 5}),  // 1: x-max, right
    face({0, 1, 5, 4}),  // 2: y-min, bottom
    face({3, 7, 6, 2}),  // 3: y-max, top
    face({0, 3, 2, 1}),  // 4: z-min, back
    face({4, 5, 6, 7})   // 5: z-max, front
});
*/
// START - by SBLEE
const Foam::faceList Foam::DAUboundBox::faces()
{
    faceList faces(6);
    forAll(faces, fI) faces[fI].setSize(4);
    faces[0][0] = 0;
    faces[0][1] = 4;
    faces[0][2] = 7;
    faces[0][3] = 3;

    faces[1][0] = 1;
    faces[1][1] = 2;
    faces[1][2] = 6;
    faces[1][3] = 5;

    faces[2][0] = 0;
    faces[2][1] = 1;
    faces[2][2] = 5;
    faces[2][3] = 4;

    faces[3][0] = 3;
    faces[3][1] = 7;
    faces[3][2] = 6;
    faces[3][3] = 2;

    faces[4][0] = 0;
    faces[4][1] = 3;
    faces[4][2] = 2;
    faces[4][3] = 1;

    faces[5][0] = 4;
    faces[5][1] = 5;
    faces[5][2] = 6;
    faces[5][3] = 7;
    
    return faces;
}
// END - by SBLEE

/*
const Foam::FixedList<Foam::vector, 6> Foam::DAUboundBox::faceNormals
({
    vector(-1,  0,  0), // 0: x-min, left
    vector( 1,  0,  0), // 1: x-max, right
    vector( 0, -1,  0), // 2: y-min, bottom
    vector( 0,  1,  0), // 3: y-max, top
    vector( 0,  0, -1), // 4: z-min, back
    vector( 0,  0,  1)  // 5: z-max, front
});
*/
// START - by SBLEE
const Foam::FixedList<Foam::vector, 6> Foam::DAUboundBox::faceNormals()
{
    FixedList<Foam::vector, 6> fN;
    fN[0]=vector(-1,  0,  0); // 0: x-min, left
    fN[1]=vector( 1,  0,  0); // 1: x-max, right
    fN[2]=vector( 0, -1,  0); // 2: y-min, bottom
    fN[3]=vector( 0,  1,  0); // 3: y-max, top
    fN[4]=vector( 0,  0, -1); // 4: z-min, back
    fN[5]=vector( 0,  0,  1); // 5: z-max, front
    return fN;
}
// END - by SBLEE

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DAUboundBox::DAUboundBox(const UList<point>& points, bool doReduce)
:
    DAUboundBox()
{
    add(points);

    if (doReduce)
    {
        reduce();
    }
}


Foam::DAUboundBox::DAUboundBox(const tmp<pointField>& tpoints, bool doReduce)
:
    DAUboundBox()
{
    add(tpoints);

    if (doReduce)
    {
        reduce();
    }
}


Foam::DAUboundBox::DAUboundBox
(
    const UList<point>& points,
    const labelUList& indices,
    bool doReduce
)
:
    DAUboundBox()
{
    add(points, indices);

    if (doReduce)
    {
        reduce();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::DAUboundBox::points() const
{
    //auto tpt = tmp<pointField>::New(8);
    auto tpt = tmp<pointField>(new pointField(8)); // by SBLEE
    auto& pt = tpt.ref();

    pt[0] = min_;                                   // min-x, min-y, min-z
    pt[1] = point(max_.x(), min_.y(), min_.z());    // max-x, min-y, min-z
    pt[2] = point(max_.x(), max_.y(), min_.z());    // max-x, max-y, min-z
    pt[3] = point(min_.x(), max_.y(), min_.z());    // min-x, max-y, min-z
    pt[4] = point(min_.x(), min_.y(), max_.z());    // min-x, min-y, max-z
    pt[5] = point(max_.x(), min_.y(), max_.z());    // max-x, min-y, max-z
    pt[6] = max_;                                   // max-x, max-y, max-z
    pt[7] = point(min_.x(), max_.y(), max_.z());    // min-x, max-y, max-z

    return tpt;
}


Foam::tmp<Foam::pointField> Foam::DAUboundBox::faceCentres() const
{
    //auto tpts = tmp<pointField>::New(6);
    auto tpts = tmp<pointField>(new pointField(6)); // by SBLEE
    auto& pts = tpts.ref();

    forAll(pts, facei)
    {
        pts[facei] = faceCentre(facei);
    }

    return tpts;
}


Foam::point Foam::DAUboundBox::faceCentre(const direction facei) const
{
    point pt = DAUboundBox::centre();

    if (facei > 5)
    {
        FatalErrorInFunction
            << "face should be [0..5]"
            << abort(FatalError);
    }

    switch (facei)
    {
        case 0: pt.x() = min().x(); break;  // 0: x-min, left
        case 1: pt.x() = max().x(); break;  // 1: x-max, right
        case 2: pt.y() = min().y(); break;  // 2: y-min, bottom
        case 3: pt.y() = max().y(); break;  // 3: y-max, top
        case 4: pt.z() = min().z(); break;  // 4: z-min, back
        case 5: pt.z() = max().z(); break;  // 5: z-max, front
    }

    return pt;
}


void Foam::DAUboundBox::inflate(const scalar s)
{
    const vector ext = vector::one*s*mag();

    min_ -= ext;
    max_ += ext;
}


void Foam::DAUboundBox::reduce()
{
    Foam::reduce(min_, minOp<point>());
    Foam::reduce(max_, maxOp<point>());
}


bool Foam::DAUboundBox::intersect(const DAUboundBox& bb)
{
    min_ = ::Foam::max(min_, bb.min_);
    max_ = ::Foam::min(max_, bb.max_);

    return valid();
}


bool Foam::DAUboundBox::intersects(const plane& pln) const
{
    // Require a full 3D box
    if (nDim() != 3)
    {
        return false;
    }

    bool above = false;
    bool below = false;

    tmp<pointField> tpts(points());
    const auto& pts = tpts();

    for (const point& p : pts)
    {
        //if (pln.sideOfPlane(p) == plane::FRONT)
        if (pln.sideOfPlane(p) == plane::NORMAL) // by SBLEE
        {
            above = true;
        }
        else
        {
            below = true;
        }
    }

    return (above && below);
}


bool Foam::DAUboundBox::contains(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (!contains(p))
        {
            return false;
        }
    }

    return true;
}


bool Foam::DAUboundBox::containsAny(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (contains(p))
        {
            return true;
        }
    }

    return false;
}


Foam::point Foam::DAUboundBox::nearest(const point& pt) const
{
    // Clip the point to the range of the bounding box
    const scalar surfPtx = Foam::max(Foam::min(pt.x(), max_.x()), min_.x());
    const scalar surfPty = Foam::max(Foam::min(pt.y(), max_.y()), min_.y());
    const scalar surfPtz = Foam::max(Foam::min(pt.z(), max_.z()), min_.z());

    return point(surfPtx, surfPty, surfPtz);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const DAUboundBox& bb)
{
    if (os.format() == IOstream::ASCII)
    {
        os << bb.min_ << token::SPACE << bb.max_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&bb.min_),
            sizeof(DAUboundBox)
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, DAUboundBox& bb)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> bb.min_ >> bb.max_;
    }
    else
    {
        /*
        Detail::readContiguous<DAUboundBox>
        (
            is,
            reinterpret_cast<char*>(&bb.min_),
            sizeof(DAUboundBox)
        );
        */
        is.read(reinterpret_cast<char*>(&bb.min_), sizeof(DAUboundBox)); // by SBLEE
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
