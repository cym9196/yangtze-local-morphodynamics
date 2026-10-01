// 中文说明：求解 Exner 床面守恒方程，并限制单步床面变化以保证毫米网格稳定。
/*---------------------------------------------------------------------------*\
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of ScourFOAM.

    ScourFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ScourFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ScourFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

const faMesh& aMesh = bed.aMesh();

areaScalarField& dH = dHptr.ref();
areaScalarField& dVsed = dVsedPtr.ref();

// faMesh horizontal projection
const projectedFaMesh& meshProj = bed.aProjMesh();
// area of projected mesh faces
const scalarField& SfProj = meshProj.S();

#include "balanceBedload.H"

//#include "bedFormVelocity.H"

if (switchSuspension=="on")
{
    #include "erosionDeposition.H"
}

scalar deltaT = runTime.time().deltaT().value();  // current time step

// solve exner equation explicitly with Adams-Bashforth scheme
if (ABorder==0)
{
    // Adams-Bashforth order 0, euler explicit scheme
    forAll(aMesh.areaCentres(), facei)
    {
        dH[facei] =
            morphoAccFactor
            * (1/CsMax.value())
            * (deltaT / SfProj[facei])
            * dVsed[facei];
    }
}
else if (ABorder==1)
{
    scalar deltaT0 = runTime.time().deltaT0().value();  // old time step

    scalar coefT0 = deltaT * (1 + 0.5 * (deltaT / deltaT0));
    scalar coefT00 = - 0.5 * deltaT * deltaT / deltaT0;

    // Adams-Bashforth order 1
    forAll(aMesh.areaCentres(), facei)
    {
        dH[facei] =
            morphoAccFactor
            * (1/CsMax.value())
            * (coefT0 * dVsed[facei] + coefT00 * dVsed.oldTime()[facei])
            * (1 / SfProj[facei]);
    }
}

if (filterExner)
{
    labelList neighbourCount(aMesh.nFaces(), Zero);
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        ++neighbourCount[owner[edgei]];
        ++neighbourCount[neighbour[edgei]];
    }

    // Equal and opposite volume exchanges preserve sum(SfProj*dH).
    for (label pass = 0; pass < nFiltExner; ++pass)
    {
        scalarField volumeCorrection(aMesh.nFaces(), Zero);

        for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
        {
            const label own = owner[edgei];
            const label nei = neighbour[edgei];
            const scalar degree = max
            (
                scalar(neighbourCount[own]),
                scalar(neighbourCount[nei])
            );
            const scalar exchangedVolume =
                alphaFiltExner
               *min(SfProj[own], SfProj[nei])
               *(dH[nei] - dH[own])/(degree + VSMALL);

            volumeCorrection[own] += exchangedVolume;
            volumeCorrection[nei] -= exchangedVolume;
        }

        forAll(dH, facei)
        {
            dH[facei] += volumeCorrection[facei]/SfProj[facei];
        }
    }
}

// Morphodynamic CFL guard.  The transport equations remain conservative in
// their normal range; only a source spike large enough to move the bed by more
// than 1% of a 1 mm cell in one CFD step is clipped.
scalar maxRawBedChange = 0;
label limitedBedFaces = 0;
const scalar activeBedChangeLimit =
    min(maxBedChangePerStep, maxBedChangeRate*deltaT);
forAll(dH, facei)
{
    const scalar rawBedChange = dH[facei];
    maxRawBedChange = max(maxRawBedChange, mag(rawBedChange));
    dH[facei] =
        min
        (
            max(rawBedChange, -activeBedChangeLimit),
            activeBedChangeLimit
        );
    if (mag(dH[facei] - rawBedChange) > SMALL)
    {
        ++limitedBedFaces;
    }
}
reduce(maxRawBedChange, maxOp<scalar>());
reduce(limitedBedFaces, sumOp<label>());
if (limitedBedFaces > 0)
{
    Info<< "morphodynamic limiter: faces=" << limitedBedFaces
        << " rawMax=" << maxRawBedChange
        << " limitedMax=" << activeBedChangeLimit
        << " m (rate=" << maxBedChangeRate << " m/s)" << endl;
}

if (bed.bedMotion())
{
    #include "meshMove.H"
}
