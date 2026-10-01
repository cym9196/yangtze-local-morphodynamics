// 中文说明：把有限面积床面位移传递到体网格并执行动态网格移动。
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

pointVectorField& pointDisplacement = const_cast<pointVectorField&>
(
    mesh.objectRegistry::lookupObject<pointVectorField>
    (
        "pointDisplacement"
    )
);

vectorField& dispVals = dispValsPtr.ref();
for (const label patchID : bed.bedPatchesID())
{
    // A point field stores patch-point values in its internal point array;
    // patchInternalField() extracts those actual boundary point values.
    dispVals =
        pointDisplacement.boundaryField()[patchID].patchInternalField();
}


scalarField& dHpoints = dHpointsPtr.ref();
dHpoints = Zero;

//interpolate bed level increment from face center to vertices
bed.interpFaceToVertices(dH, dHpoints);

// zero displacements for faPatches
const labelList& fixedPointsLabels = fixedPatchPointLabelsPtr.ref();

for (const label i : fixedPointsLabels)
{
    dHpoints[i] = 0.;
}

// Candidate absolute vertical displacement of every finite-area bed point.
// Positive eta raises the bed because eg points in the gravity direction.
scalarField eta0(dHpoints.size(), Zero);
forAll(eta0, pointi)
{
    eta0[pointi] = -(eg & dispVals[pointi]) + dHpoints[pointi];
}

const uindirectPrimitivePatch& bedPatch = bed.aMesh().patch();
const labelList& bedMeshPoints = bedPatch.meshPoints();
const edgeList& bedEdges = bedPatch.edges();
const pointField& bedPoints = bedPatch.localPoints();

labelList fixedBedPoint(eta0.size(), Zero);
for (const label pointi : fixedPointsLabels)
{
    fixedBedPoint[pointi] = 1;
    eta0[pointi] = 0;
}

// Processor patches duplicate finite-area points.  Synchronising through the
// parent polyMesh makes the limiter decomposition-independent and keeps the
// moving bed watertight across MPI partitions.
auto syncBedPointAverage =
[
    &mesh,
    &bedMeshPoints
]
(scalarField& values)
{
    scalarField meshValues(mesh.nPoints(), Zero);
    labelList meshCounts(mesh.nPoints(), Zero);

    forAll(values, pointi)
    {
        const label meshPointi = bedMeshPoints[pointi];
        meshValues[meshPointi] += values[pointi];
        ++meshCounts[meshPointi];
    }

    syncTools::syncPointList
    (
        mesh,
        meshValues,
        plusEqOp<scalar>(),
        scalar(0)
    );
    syncTools::syncPointList
    (
        mesh,
        meshCounts,
        plusEqOp<label>(),
        label(0)
    );

    forAll(values, pointi)
    {
        const label meshPointi = bedMeshPoints[pointi];
        values[pointi] =
            meshValues[meshPointi]/max(meshCounts[meshPointi], label(1));
    }
};

syncBedPointAverage(eta0);
for (const label pointi : fixedPointsLabels)
{
    eta0[pointi] = 0;
}

scalar initialMaxSlope = 0;
forAll(bedEdges, edgei)
{
    const edge& e = bedEdges[edgei];
    const vector edgeVector = bedPoints[e.end()] - bedPoints[e.start()];
    const scalar horizontalLength =
        mag(edgeVector - eg*(eg & edgeVector));
    if (horizontalLength > VSMALL)
    {
        initialMaxSlope = max
        (
            initialMaxSlope,
            mag(eta0[e.end()] - eta0[e.start()])/horizontalLength
        );
    }
}
reduce(initialMaxSlope, maxOp<scalar>());

const label limiterPasses =
    initialMaxSlope > 1.2*maxBedSlope
  ? bedSurfaceRepairPasses
  : bedSurfaceLimiterPasses;

scalarField etaTarget(eta0);
for (label pass = 0; pass < limiterPasses; ++pass)
{
    scalarField neighbourSum(etaTarget.size(), Zero);
    labelList neighbourCount(etaTarget.size(), Zero);

    forAll(bedEdges, edgei)
    {
        const edge& e = bedEdges[edgei];
        neighbourSum[e.start()] += etaTarget[e.end()];
        neighbourSum[e.end()] += etaTarget[e.start()];
        ++neighbourCount[e.start()];
        ++neighbourCount[e.end()];
    }

    scalarField correction(etaTarget.size(), Zero);
    labelList correctionCount(etaTarget.size(), Zero);

    // Limit point-scale curvature relative to the immediate bed neighbours.
    forAll(etaTarget, pointi)
    {
        if (!fixedBedPoint[pointi] && neighbourCount[pointi] > 0)
        {
            const scalar neighbourAverage =
                neighbourSum[pointi]/neighbourCount[pointi];
            const scalar deviation = etaTarget[pointi] - neighbourAverage;
            const scalar limitedDeviation = min
            (
                max(deviation, -maxBedPointCurvature),
                maxBedPointCurvature
            );
            correction[pointi] +=
                neighbourAverage + limitedDeviation - etaTarget[pointi];
            ++correctionCount[pointi];
        }
    }

    // Enforce a slope below the 32 degree repose slope (tan(32)=0.625)
    // with numerical headroom.  Pair corrections are equal and opposite
    // unless one endpoint belongs to the fixed flume perimeter.
    forAll(bedEdges, edgei)
    {
        const edge& e = bedEdges[edgei];
        const label a = e.start();
        const label b = e.end();
        const vector edgeVector = bedPoints[b] - bedPoints[a];
        const scalar horizontalLength =
            mag(edgeVector - eg*(eg & edgeVector));
        const scalar allowedDifference = maxBedSlope*horizontalLength;
        const scalar difference = etaTarget[b] - etaTarget[a];
        const scalar excess = mag(difference) - allowedDifference;

        if (horizontalLength > VSMALL && excess > 0)
        {
            const scalar signDifference = difference >= 0 ? 1 : -1;

            if (!fixedBedPoint[a])
            {
                correction[a] += signDifference*excess
                    *(fixedBedPoint[b] ? 1 : 0.5);
                ++correctionCount[a];
            }
            if (!fixedBedPoint[b])
            {
                correction[b] -= signDifference*excess
                    *(fixedBedPoint[a] ? 1 : 0.5);
                ++correctionCount[b];
            }
        }
    }

    forAll(etaTarget, pointi)
    {
        if (correctionCount[pointi] > 0)
        {
            etaTarget[pointi] +=
                correction[pointi]/correctionCount[pointi];
        }
    }
}

syncBedPointAverage(etaTarget);
for (const label pointi : fixedPointsLabels)
{
    etaTarget[pointi] = 0;
}

// Repair an inherited rough restart gradually.  The 0.01 m/s default gives a
// mesh Courant number near 0.002 at dt=2e-4 s while removing a millimetre
// point defect in O(0.1 s), well before it can fold a bed face.
const scalar maxRegularisationCorrection =
    maxBedRegularisationRate*runTime.deltaTValue();
scalar maxAppliedRegularisation = 0;
scalarField eta(eta0);
forAll(eta, pointi)
{
    const scalar requestedCorrection = etaTarget[pointi] - eta0[pointi];
    const scalar appliedCorrection = min
    (
        max(requestedCorrection, -maxRegularisationCorrection),
        maxRegularisationCorrection
    );
    eta[pointi] += appliedCorrection;
    maxAppliedRegularisation = max
    (
        maxAppliedRegularisation,
        mag(appliedCorrection)
    );
}
reduce(maxAppliedRegularisation, maxOp<scalar>());

// Remove the global point-mean correction over movable points.  On this
// uniform Cartesian bed it is the discrete analogue of conserving sediment
// volume during the numerical repose/curvature projection.
scalar regularisationSum = 0;
label movablePointCount = 0;
forAll(eta, pointi)
{
    if (!fixedBedPoint[pointi])
    {
        regularisationSum += eta[pointi] - eta0[pointi];
        ++movablePointCount;
    }
}
reduce(regularisationSum, sumOp<scalar>());
reduce(movablePointCount, sumOp<label>());
const scalar meanRegularisation =
    regularisationSum/max(movablePointCount, label(1));
forAll(eta, pointi)
{
    if (!fixedBedPoint[pointi])
    {
        eta[pointi] -= meanRegularisation;
    }
}

syncBedPointAverage(eta);
for (const label pointi : fixedPointsLabels)
{
    eta[pointi] = 0;
}

scalar finalMaxSlope = 0;
scalarField finalNeighbourSum(eta.size(), Zero);
labelList finalNeighbourCount(eta.size(), Zero);
forAll(bedEdges, edgei)
{
    const edge& e = bedEdges[edgei];
    const vector edgeVector = bedPoints[e.end()] - bedPoints[e.start()];
    const scalar horizontalLength =
        mag(edgeVector - eg*(eg & edgeVector));
    if (horizontalLength > VSMALL)
    {
        finalMaxSlope = max
        (
            finalMaxSlope,
            mag(eta[e.end()] - eta[e.start()])/horizontalLength
        );
    }

    finalNeighbourSum[e.start()] += eta[e.end()];
    finalNeighbourSum[e.end()] += eta[e.start()];
    ++finalNeighbourCount[e.start()];
    ++finalNeighbourCount[e.end()];
}
reduce(finalMaxSlope, maxOp<scalar>());

scalar finalMaxPointCurvature = 0;
forAll(eta, pointi)
{
    if (!fixedBedPoint[pointi] && finalNeighbourCount[pointi] > 0)
    {
        finalMaxPointCurvature = max
        (
            finalMaxPointCurvature,
            mag
            (
                eta[pointi]
              - finalNeighbourSum[pointi]/finalNeighbourCount[pointi]
            )
        );
    }
}
reduce(finalMaxPointCurvature, maxOp<scalar>());

Info<< "bed-surface limiter: slope=" << initialMaxSlope
    << " -> " << finalMaxSlope
    << ", curvature=" << finalMaxPointCurvature << " m"
    << ", passes=" << limiterPasses
    << ", maxCorrection=" << maxAppliedRegularisation << " m"
    << ", meanCorrectionRemoved=" << meanRegularisation << " m" << endl;

// Preserve any non-vertical component and replace only the absolute vertical
// bed displacement with the bounded accumulated surface.
forAll(dispVals, pointi)
{
    dispVals[pointi] -= eg*(eta[pointi] + (eg & dispVals[pointi]));
}

// interpolate displacement on points, vertical orientation
//dispVals -= patchInterp.faceToPointInterpolate(dH) * eg;

for (const label patchID : bed.bedPatchesID())
{
    Info << "mesh vertices motion, max:" << max(mag(dispVals));
    Info << " m , min:" << min(mag(dispVals));
    Info << " m , mean:" << average(mag(dispVals)) << " m" << endl;
    pointDisplacement.boundaryFieldRef()[patchID] == dispVals;
}
