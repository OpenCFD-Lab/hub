/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::calculatedProcessorGAMGInterfaceField

Description
    GAMG agglomerated processor interface field.

SourceFiles
    calculatedProcessorGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef calculatedProcessorGAMGInterfaceField_H
#define calculatedProcessorGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "calculatedProcessorGAMGInterface.H"
#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class calculatedProcessorGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class calculatedProcessorGAMGInterfaceField
:
    public GAMGInterfaceField,
    public processorLduInterfaceField
{
    // Private data
        typedef scalarField solveScalarField; // by SBLEE

        //- Local reference cast into the processor interface
        const calculatedProcessorGAMGInterface& procInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;


        // Sending and receiving

            //- Outstanding request
            mutable label outstandingSendRequest_;

            //- Outstanding request
            mutable label outstandingRecvRequest_;

            //- Scalar send buffer
            mutable solveScalarField scalarSendBuf_;

            //- Scalar receive buffer
            mutable solveScalarField scalarReceiveBuf_;



    // Private Member Functions

        //- No copy construct
        calculatedProcessorGAMGInterfaceField
        (
            const calculatedProcessorGAMGInterfaceField&
        ) = delete;

        //- No copy assignment
        void operator=(const calculatedProcessorGAMGInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("calculatedProcessor");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        calculatedProcessorGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterface
        );

        //- Construct from GAMG interface and fine level interface field
        calculatedProcessorGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const bool doTransform,
            const int rank
        );


    //- Destructor
    virtual ~calculatedProcessorGAMGInterfaceField() = default;


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return procInterface_.size();
            }


        // Interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            virtual void updateInterfaceMatrix
            (
                scalarField&,
                const scalarField&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const
            {
                NotImplemented;
            }


        //- Processor interface functions

            //- Return communicator used for comms
            virtual label comm() const
            {
                return procInterface_.comm();
            }

            //- Return processor number
            virtual int myProcNo() const
            {
                return procInterface_.myProcNo();
            }

            //- Return neighbour processor number
            virtual int neighbProcNo() const
            {
                return procInterface_.neighbProcNo();
            }

            //- Does the interface field perform the transformation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return procInterface_.forwardT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
            }

            // START - by SBLEE
            template<class Type>
            void addToInternalField
            (
                Field<Type>& result,
                const bool add,
                const labelUList& faceCells,
                const scalarField& coeffs,
                const Field<Type>& vals
            ) const
            {
                if (add)
                {
                    forAll(faceCells, elemi)
                    {
                        result[faceCells[elemi]] += coeffs[elemi]*vals[elemi];
                    }
                }
                else
                {
                    forAll(faceCells, elemi)
                    {
                        result[faceCells[elemi]] -= coeffs[elemi]*vals[elemi];
                    }
                }
            }
            // END - by SBLEE
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
