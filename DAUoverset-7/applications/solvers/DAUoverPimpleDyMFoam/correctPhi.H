if (mesh.changing())
{
    volVectorField::Boundary& bfld = U.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (bfld[patchi].fixesValue())
        {
            bfld[patchi].initEvaluate();
        }
    }

    surfaceScalarField::Boundary& phiBfld = phi.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (bfld[patchi].fixesValue())
        {
            bfld[patchi].evaluate();

            phiBfld[patchi] = bfld[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
}
    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        //dimensionedScalar(p.dimensions(), Zero),
        dimensionedScalar("",p.dimensions(), Zero),
        pcorrTypes
    );

{
    dimensionedScalar rAUf("rAUf", dimTime, 1.0);

    const cellCellStencilObject& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();
    const labelIOList& zoneIDs = overlap.zoneID();

    while (pimple.correctNonOrthogonal())
    {
        label nZones = gMax(zoneIDs)+1;

        //label refCellI2 = -1;
        labelList refCells(nZones, -1);
        labelList refZones(nZones, -1);

        forAll(zoneIDs, cellI)
        {
            label zoneId = zoneIDs[cellI];
            if
            (
                refCells[zoneId] == -1
             && cellTypes[cellI] == cellCellStencil::CALCULATED
             && refZones[zoneId] == -1
            )
            {
                refCells[zoneId] = cellI;
                refZones[zoneId] = zoneId;
            }
        }
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );

        // Only set reference for cells that are CALCULATED
        {
            DynamicList<label> validCells(refCells.size());
            forAll(refCells, zoneId)
            {
                if (refCells[zoneId] != -1)
                {
                    validCells.append(refCells[zoneId]);
                }
            }
            //pcorrEqn.setReferences(validCells, scalar(0), true);
            // START - by SBLEE
            //if(pcorr.needReference())  // by CWKIM
            //{
                forAll(validCells, i)
                {
                    label cellI=validCells[i]; 
                    pcorrEqn.source()[cellI]+=pcorrEqn.diag()[cellI]*scalar(0); // dummy
                    pcorrEqn.diag()[cellI] += pcorrEqn.diag()[cellI];
                }
            //}
            // END - by SBLEE
        }
        const dictionary& d = mesh.solverDict // by SBLEE
        (
            pcorr.select
            (
                pimple.finalInnerIter()
            )
        );
        //Bypass virtual layer
//        mesh.fvMesh::solve(pcorrEqn, d);
        pcorrEqn.solve(d); // by CWKIM

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }

    }
}
