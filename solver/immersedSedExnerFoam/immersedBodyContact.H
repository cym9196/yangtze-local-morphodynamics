// 中文说明：计算铁块与实时冲淤床面的多点 Winkler–Coulomb 接触、摩擦和限位力。
/*---------------------------------------------------------------------------*\
  Parallel-safe multi-point Winkler/Coulomb contact with the Exner bed.
\*---------------------------------------------------------------------------*/

// Integrate the smoothed mask only over the water mesh.  Unlike a fixed full
// cuboid volume, this remains correct if the block later changes its
// submergence or orientation.
scalar currentDisplacedVolume = 0;
forAll(cellVolumes, cellI)
{
    currentDisplacedVolume += bodyMask[cellI]*cellVolumes[cellI];
}
reduce(currentDisplacedVolume, sumOp<scalar>());

contactForce = vector
(
    0,
    0,
    contactFluidDensity*currentDisplacedVolume*mag(g.value())
);
contactMoment = Zero;

const vector verticalNormal(0, 0, 1);
const vectorField& bedCentres = bed.aMesh().areaCentres();
const vectorField& bedAreaNormals = bed.aMesh().faceAreaNormals();

// Project only downward-facing surface patches onto the horizontal bed.  The
// summed projected area is 0.01 m2 upright and 0.02 m2 on either long side.
scalarField currentPointAreas(refContactPoints.size(), Zero);
scalar currentProjectedArea = 0;

forAll(refContactPoints, pointI)
{
    if (useDirectionalContactQuadrature)
    {
        const vector currentNormal =
            bodyMotion.orientation() & refContactNormals[pointI];
        currentPointAreas[pointI] =
            refContactPointAreas[pointI]
           *max(-(currentNormal & verticalNormal), scalar(0));
    }
    else
    {
        currentPointAreas[pointI] =
            contactArea/scalar(refContactPoints.size());
    }
    currentProjectedArea += currentPointAreas[pointI];
}

forAll(refContactPoints, pointI)
{
    const scalar pointArea = currentPointAreas[pointI];
    if (pointArea <= SMALL)
    {
        continue;
    }

    const scalar pointMass =
        bodyMotion.mass()*pointArea/(currentProjectedArea + VSMALL);
    const scalar pointStiffness = foundationStiffness*pointArea;
    const scalar pointDamping =
        2.0*contactDampingRatio
       *Foam::sqrt(doubleScalar(max(pointStiffness*pointMass, SMALL)));
    const scalar maxPointForce = bearingCapacity*pointArea;

    const point contactPoint = bodyMotion.transform(refContactPoints[pointI]);
    const vector contactVelocity = bodyMotion.velocity(contactPoint);

    scalar localDistanceSqr = GREAT;
    vector localBedCentre(Zero);
    vector localBedNormal(verticalNormal);

    forAll(bedCentres, bedFaceI)
    {
        const scalar dx = contactPoint.x() - bedCentres[bedFaceI].x();
        const scalar dy = contactPoint.y() - bedCentres[bedFaceI].y();
        const scalar distanceSqr = sqr(dx) + sqr(dy);

        if (distanceSqr < localDistanceSqr)
        {
            localDistanceSqr = distanceSqr;
            localBedCentre = bedCentres[bedFaceI];
            // Finite-area bed normals follow the volume-patch orientation;
            // enforce an upward contact normal for support forces.
            localBedNormal = -bedAreaNormals[bedFaceI];
            if ((localBedNormal & verticalNormal) < 0)
            {
                localBedNormal = -localBedNormal;
            }
        }
    }

    scalar globalDistanceSqr = localDistanceSqr;
    reduce(globalDistanceSqr, minOp<scalar>());

    const bool ownsNearest =
        mag(localDistanceSqr - globalDistanceSqr)
      < max(SMALL, 1.0e-10*max(globalDistanceSqr, scalar(1)));

    vector bedCentreContribution = ownsNearest ? localBedCentre : vector(Zero);
    vector bedNormalContribution = ownsNearest ? localBedNormal : vector(Zero);
    label nearestCount = ownsNearest ? 1 : 0;
    reduce(bedCentreContribution, sumOp<vector>());
    reduce(bedNormalContribution, sumOp<vector>());
    reduce(nearestCount, sumOp<label>());

    if
    (
        nearestCount == 0
     || globalDistanceSqr > sqr(contactSearchRadius)
    )
    {
        continue;
    }

    const vector nearestBedCentre =
        bedCentreContribution/scalar(nearestCount);
    vector contactNormal = bedNormalContribution/scalar(nearestCount);
    contactNormal /= mag(contactNormal) + VSMALL;

    // Signed distance to the local tangent plane is more accurate than a
    // vertical-only gap on the sloping sides of a scour hole.
    const scalar penetration =
        (nearestBedCentre - contactPoint) & contactNormal;

    if (penetration <= 0)
    {
        continue;
    }

    const scalar normalVelocity = contactVelocity & contactNormal;
    scalar normalForce =
        pointStiffness*penetration - pointDamping*normalVelocity;
    normalForce = min(max(normalForce, scalar(0)), maxPointForce);

    const vector tangentialVelocity =
        contactVelocity - normalVelocity*contactNormal;
    const scalar tangentialSpeed = mag(tangentialVelocity);
    vector tangentialForce(Zero);

    if (tangentialSpeed > SMALL)
    {
        const scalar mu =
            dynamicFriction
          + (staticFriction - dynamicFriction)
           *Foam::exp(-sqr(tangentialSpeed/(frictionVelocity + VSMALL)));

        tangentialForce =
            -mu*normalForce
            *Foam::tanh
            (
                doubleScalar(tangentialSpeed/(frictionVelocity + VSMALL))
            )
            *tangentialVelocity/(tangentialSpeed + VSMALL);
    }

    const vector pointForce = normalForce*contactNormal + tangentialForce;
    contactForce += pointForce;
    contactMoment +=
        (contactPoint - bodyMotion.centreOfMass()) ^ pointForce;
}

// Body-only retaining frame at the horizontal test-section limits.  The
// cuboid remains free to translate, roll and overturn inside the frame, but
// cannot leave the 0.05 x 0.05 m bed and subsequently fall without contact.
if (useBodyRetainer)
{
    const scalar retainerDamping =
        2.0*retainerDampingRatio
       *Foam::sqrt
        (
            doubleScalar(max(retainerStiffness*bodyMotion.mass(), SMALL))
        );

    for (direction axis = 0; axis < 3; ++axis)
    {
        if (retainerActiveAxes[axis] <= 0.5)
        {
            continue;
        }

        for (label side = 0; side < 2; ++side)
        {
            scalar extreme = side == 0 ? GREAT : -GREAT;
            point extremePoint(Zero);

            forAll(refContactPoints, pointI)
            {
                const point currentPoint =
                    bodyMotion.transform(refContactPoints[pointI]);
                if
                (
                    (side == 0 && currentPoint[axis] < extreme)
                 || (side == 1 && currentPoint[axis] > extreme)
                )
                {
                    extreme = currentPoint[axis];
                    extremePoint = currentPoint;
                }
            }

            vector wallNormal(Zero);
            wallNormal[axis] = side == 0 ? 1.0 : -1.0;
            const scalar penetration = side == 0
                ? retainerMin[axis] - extreme
                : extreme - retainerMax[axis];

            if (penetration <= 0)
            {
                continue;
            }

            const scalar normalVelocity =
                bodyMotion.velocity(extremePoint) & wallNormal;
            scalar normalForce =
                retainerStiffness*penetration
              - retainerDamping*normalVelocity;
            normalForce =
                min(max(normalForce, scalar(0)), retainerMaxForce);

            const vector frameForce = normalForce*wallNormal;
            contactForce += frameForce;
            contactMoment +=
                (extremePoint - bodyMotion.centreOfMass()) ^ frameForce;
        }
    }
}
