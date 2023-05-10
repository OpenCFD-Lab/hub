/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "dynamicOversetFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "cellCellStencilObject.H"
#include "zeroGradientFvPatchFields.H"
#include "lduPrimitiveProcessorInterface.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicOversetFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicOversetFvMesh, IOobject);
    //addToRunTimeSelectionTable(dynamicFvMesh, dynamicOversetFvMesh, doInit); // by SBLEE
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicOversetFvMesh::updateAddressing() const
{
    const cellCellStencilObject& overlap = Stencil::New(*this);

    // The (processor-local part of the) stencil determines the local
    // faces to add to the matrix. tbd: parallel
    const labelListList& stencil = overlap.cellStencil();

    // Get the base addressing
    //const lduAddressing& baseAddr = dynamicMotionSolverFvMesh::lduAddr();
    const lduAddressing& baseAddr = dynamicFvMesh::lduAddr();

    // Add to the base addressing
    labelList lowerAddr;
    labelList upperAddr;
    label nExtraFaces;

    const globalIndex globalNumbering(baseAddr.size());
    labelListList localFaceCells;
    labelListList remoteFaceCells;

    labelList globalCellIDs(overlap.cellInterpolationMap().constructSize());
    forAll(baseAddr, cellI)
    {
        globalCellIDs[cellI] = globalNumbering.toGlobal(cellI);
    }
    overlap.cellInterpolationMap().distribute(globalCellIDs);


    reverseFaceMap_ = fvMeshPrimitiveLduAddressing::addAddressing
    (
        baseAddr,
        stencil,
        nExtraFaces,
        lowerAddr,
        upperAddr,
        stencilFaces_,
        globalNumbering,
        globalCellIDs,
        localFaceCells,
        remoteFaceCells
    );

    // Now for the tricky bits. We want to hand out processor faces according
    // to the localFaceCells/remoteFaceCells. Ultimately we need
    // per entry in stencil:
    // - the patch (or -1 for internal faces)
    // - the face (is either an internal face index or a patch face index)

    stencilPatches_.setSize(stencilFaces_.size());

    // Per processor to owner (local)/neighbour (remote)
    List<DynamicList<label>> procOwner(Pstream::nProcs());
    List<DynamicList<label>> dynProcNeighbour(Pstream::nProcs());
    forAll(stencil, celli)
    {
        const labelList& nbrs = stencil[celli];
        stencilPatches_[celli].setSize(nbrs.size());
        stencilPatches_[celli] = -1;

        forAll(nbrs, nbri)
        {
            if (stencilFaces_[celli][nbri] == -1)
            {
                const label nbrCelli = nbrs[nbri];
                label globalNbr = globalCellIDs[nbrCelli];
                label proci = globalNumbering.whichProcID(globalNbr);
                label remoteCelli = globalNumbering.toLocal(proci, globalNbr);

                // Overwrite the face to be a patch face
                stencilFaces_[celli][nbri] = procOwner[proci].size();
                stencilPatches_[celli][nbri] = proci;
                procOwner[proci].append(celli);
                dynProcNeighbour[proci].append(remoteCelli);

                //Pout<< "From neighbour proc:" << proci
                //    << " allocating patchFace:" << stencilFaces_[celli][nbri]
                //    << " to get remote cell " << remoteCelli
                //    << " onto local cell " << celli << endl;
            }
        }
    }
    labelListList procNeighbour(dynProcNeighbour.size());
    forAll(procNeighbour, i)
    {
        procNeighbour[i] = std::move(dynProcNeighbour[i]);
    }
    labelListList mySendCells;
    Pstream::exchange<labelList, label>(procNeighbour, mySendCells);

    label nbri = 0;
    forAll(procOwner, proci)
    {
        if (procOwner[proci].size())
        {
            nbri++;
        }
        if (mySendCells[proci].size())
        {
            nbri++;
        }
    }
    remoteStencilInterfaces_.setSize(nbri);
    nbri = 0;

    // E.g. if proc1 needs some data from proc2 and proc2 needs some data from
    //      proc1. We first want the pair : proc1 receive and proc2 send
    //      and then the pair : proc1 send, proc2 receive


    labelList procToInterface(Pstream::nProcs(), -1);

    forAll(procOwner, proci)
    {
        if (proci < Pstream::myProcNo() && procOwner[proci].size())
        {
            procToInterface[proci] = nbri;
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    procOwner[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+2
                )
            );
        }
        else if (proci > Pstream::myProcNo() && mySendCells[proci].size())
        {
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    mySendCells[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+2
                )
            );
        }
    }
    forAll(procOwner, proci)
    {
        if (proci > Pstream::myProcNo() && procOwner[proci].size())
        {
            procToInterface[proci] = nbri;
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    procOwner[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+3
                )
            );
        }
        else if (proci < Pstream::myProcNo() && mySendCells[proci].size())
        {
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    mySendCells[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+3
                )
            );
        }
    }


    // Rewrite stencilPatches now we know the actual interface (procToInterface)
    for (auto& patches : stencilPatches_)
    {
        for (auto& interface : patches)
        {
            if (interface != -1)
            {
                interface = procToInterface[interface]+boundary().size();
            }
        }
    }


    // Get addressing and interfaces of all interfaces


    UPtrList<const labelUList> patchAddr;
    {
        const fvBoundaryMesh& fvp = boundary();

        patchAddr.setSize(fvp.size() + remoteStencilInterfaces_.size());

        //allInterfaces_ = dynamicMotionSolverFvMesh::interfaces();
        allInterfaces_ = dynamicFvMesh::interfaces();
        allInterfaces_.setSize(patchAddr.size());

        forAll(fvp, patchi)
        {
            patchAddr.set(patchi, &fvp[patchi].faceCells());
        }
        forAll(remoteStencilInterfaces_, i)
        {
            const label patchi = fvp.size()+i;
            const lduPrimitiveProcessorInterface& pp =
                remoteStencilInterfaces_[i];

            //Pout<< "at patch:" << patchi
            //    << " have procPatch:" << pp.type()
            //    << " from:" << pp.myProcNo()
            //    << " to:" << pp.neighbProcNo()
            //    << " with fc:" << pp.faceCells().size() << endl;

            patchAddr.set(patchi, &pp.faceCells());
            allInterfaces_.set(patchi, &pp);
        }
    }
    const lduSchedule ps
    (
        lduPrimitiveMesh::nonBlockingSchedule<processorLduInterface>
        (
            allInterfaces_
        )
    );

    lduPtr_.reset
    (
        new fvMeshPrimitiveLduAddressing
        (
            nCells(),
            std::move(lowerAddr),
            std::move(upperAddr),
            patchAddr,
            ps
        )
   );


    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicOversetFvMesh::dynamicOversetFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicMotionSolverListFvMesh(io, doInit)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicOversetFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        dynamicMotionSolverListFvMesh::init(doInit);
    }

    active_ = false;

    // Load stencil (but do not update)
    (void)Stencil::New(*this, false);

    // Assume something changed
    return true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicOversetFvMesh::~dynamicOversetFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduAddressing& Foam::dynamicOversetFvMesh::lduAddr() const
{
    if (!active_)
    {
        //return dynamicMotionSolverFvMesh::lduAddr();
        return dynamicFvMesh::lduAddr();
    }
    //if (!lduPtr_)
    if (!(lduPtr_.valid())) // by SBLEE
    {
        // Build extended addressing
        updateAddressing();
    }
    //return *lduPtr_;
    return lduPtr_();  // by SBLEE
}

bool Foam::dynamicOversetFvMesh::update()
{
    //if (dynamicMotionSolverFvMesh::update())
    if (dynamicMotionSolverListFvMesh::update())
    {
        // Calculate the local extra faces for the interpolation. Note: could
        // let demand-driven lduAddr() trigger it but just to make sure.
        updateAddressing();

        // Addressing and/or weights have changed. Make interpolated cells
        // up to date with donors
        interpolateFields();

        return true;
    }

    return false;
}


Foam::word Foam::dynamicOversetFvMesh::baseName(const word& name)
{
    //if (!name.empty() && name.back()=="_0") // by SBLEE
    //if (name.ends_with("_0"))
    //{
    //    return baseName(name.substr(0, name.size()-2));
    //}

    return name;
}


bool Foam::dynamicOversetFvMesh::interpolateFields()
{
    // Add the stencil suppression list
    wordHashSet suppressed(Stencil::New(*this).nonInterpolatedFields());

    // Use whatever the solver has set up as suppression list
    const dictionary* dictPtr
    (
        //this->schemesDict().findDict("oversetInterpolationSuppressed")
        this->schemesDict().subDictPtr("oversetInterpolationSuppressed") // by SBLEE
    );
    if (dictPtr)
    {
        suppressed.insert(dictPtr->toc());
    }

    interpolate<volScalarField>(suppressed);
    interpolate<volVectorField>(suppressed);
    interpolate<volSphericalTensorField>(suppressed);
    interpolate<volSymmTensorField>(suppressed);
    interpolate<volTensorField>(suppressed);

    return true;
}



bool Foam::dynamicOversetFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    //bool ok = dynamicMotionSolverFvMesh::writeObject(streamOpt, valid);
    //bool ok = dynamicMotionSolverListFvMesh::writeObject(streamOpt, valid); // by SBLEE

    // For postprocessing : write cellTypes and zoneID
    {
        const cellCellStencilObject& overlap = Stencil::New(*this);

        const labelUList& cellTypes = overlap.cellTypes();

        volScalarField volTypes
        (
            IOobject
            (
                "cellTypes",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            //dimensionedScalar(dimless, Zero),
            dimensionedScalar("",dimless, Zero), // by SBLEE
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(volTypes.internalField(), cellI)
        {
            volTypes[cellI] = cellTypes[cellI];
        }
        volTypes.correctBoundaryConditions();
        //volTypes.writeObject(streamOpt, valid);
        volTypes.write(); // by SBLEE
    }
    {
        volScalarField volZoneID
        (
            IOobject
            (
                "zoneID",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            //dimensionedScalar(dimless, Zero),
            dimensionedScalar("",dimless, Zero),  // by SBLEE
            zeroGradientFvPatchScalarField::typeName
        );

        const cellCellStencilObject& overlap = Stencil::New(*this);
        const labelIOList& zoneID = overlap.zoneID();

        forAll(zoneID, cellI)
        {
            volZoneID[cellI] = zoneID[cellI];
        }
        volZoneID.correctBoundaryConditions();
        //volZoneID.writeObject(streamOpt, valid);
        volZoneID.write(); // by SBLEE
    }

    //return ok;
    return true; // by SBLEE
}


// ************************************************************************* //
