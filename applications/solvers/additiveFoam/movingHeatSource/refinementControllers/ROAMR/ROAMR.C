/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
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

#include "ROAMR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementControllers
{
    defineTypeNameAndDebug(ROAMR, 0);
    addToRunTimeSelectionTable(refinementController, ROAMR, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementControllers::ROAMR::ROAMR
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    uniformIntervals(sources, dict, mesh, true),
    
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    cellsPerProc_(coeffs_.lookupOrDefault<int>("cellsPerProc", 20000))
{
    //- Take first guess at interval size
    Info << "Calculating initial AMR interval size..." << endl;
    
    //- Get average cell volume and cross-sectional area
    scalar vAvg = gSum(mesh_.V()) / mesh_.nCells();
    
    //- Get approximate cross-sectional area of each beam
    scalar scanArea = 0.0;
    scalar maxLen = 0.0;
    scalar scanEnd = 0.0;
    forAll(sources_, i)
    {
        scanArea += boundingBox_ / 2.0 
                    * max
                      (
                        sources_[i].dimensions().x(),
                        sources_[i].dimensions().y()
                      )                    
                    * sources_[i].dimensions().z();
                    
        maxLen = max(sources_[i].beam().length(), maxLen);
        scanEnd = max(sources_[i].beam().endTime(), scanEnd);
    }
    
    //- Calculate number of intervals to optimize cells per processor
    scalar targetCells = Pstream::nProcs() * cellsPerProc_;
    
    if (targetCells > mesh_.nCells())
    {
        intervals_ = maxLen * scanArea / vAvg / (targetCells - mesh_.nCells())
                     * (Foam::pow(2.0, 3.0 * nLevels_) - 1.0);
    }
    
    // Enforce nInvervals >= 1
    intervals_ = max(intervals_, 1.0);
    
    Info << "Settiing initial number of intervals to: " << intervals_ << endl;
    
    //- Set number of intervals in uniformIntervals class    
    intervalSize_ = scanEnd / intervals_;
    
    Info << "Setting initial interval size to: " << intervalSize_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementControllers::ROAMR::update()
{
    //- Scale interval size based on current cells/proc
    scalar currCellsPerProc = mesh_.nCells() / Pstream::nProcs();
    
    Info << "Current cells per processor: " << currCellsPerProc << endl;
    
    if (currCellsPerProc < cellsPerProc_)
    {
        intervalSize_ =
            min
            (
                intervalSize_ * cellsPerProc_ / currCellsPerProc,
                1.1 * intervalSize_
            );
    }
    else
    {
        intervalSize_ =
            max
            (
                intervalSize_ * currCellsPerProc / cellsPerProc_,
                0.9 * intervalSize_
            );
    }
    
    //- Update refinement field using uniform intervals functions
    if (uniformIntervals::update())
    {
        Info << "Recalculating AMR interval size to target " << cellsPerProc_
         << " cells per processor." << endl;
        
        return true;
    }
    
    else
    {
        return false;
    }
}


bool Foam::refinementControllers::ROAMR::read()
{
    if (uniformIntervals::read())
    {
        refinementDict_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        refinementDict_.lookup("cellsPerProc") >> cellsPerProc_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
